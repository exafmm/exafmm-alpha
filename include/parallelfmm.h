#ifndef parallelfmm_h
#define parallelfmm_h
#include "partition.h"

//! Handles all the communication of local essential trees
class ParallelFMM : public Partition {
private:
  std::vector<int>	sendBodyCnt;                            //!< Vector of body send counts
  std::vector<int>	sendBodyDsp;                            //!< Vector of body send displacements
  std::vector<int>	recvBodyCnt;                            //!< Vector of body recv counts
  std::vector<int>	recvBodyDsp;                            //!< Vector of body recv displacements
  std::vector<int>	sendBodyRanks;                          //!< Vector of ranks to send bodies to
  std::vector<int>	sendBodyCellCnt;                        //!< Vector of send counts for cells of bodies
  std::vector<C_iter>	sendBodyCells;                          //!< Vector of cell iterators for cells of bodies to send
  std::vector<int>	sendCellCnt;                            //!< Vector of cell send counts
  std::vector<int>	sendCellDsp;                            //!< Vector of cell send displacements
  std::vector<int>	recvCellCnt;                            //!< Vector of cell recv counts
  std::vector<int>	recvCellDsp;                            //!< Vector of cell recv displacements
  std::vector<vect>	xminAll;                                //!< Buffer for gathering XMIN
  std::vector<vect>	xmaxAll;                                //!< Buffer for gathering XMAX

  JBodies sendBodies;                                           //!< Send buffer for bodies
  JBodies recvBodies;                                           //!< Recv buffer for bodies
  JCells  sendCells;                                            //!< Send buffer for cells
  JCells  recvCells;                                            //!< Recv buffer for cells

private:
//! Gather bounds of other domain
  void gatherBounds() {
    xminAll.resize(MPISIZE);                                    // Resize buffer for gathering xmin
    xmaxAll.resize(MPISIZE);                                    // Resize buffer for gathering xmax
    sendBodyCnt.resize(MPISIZE);                                // Resize vector of body send counts
    sendBodyDsp.resize(MPISIZE);                                // Resize vector of body send displacements
    recvBodyCnt.resize(MPISIZE);                                // Resize vector of body recv counts
    recvBodyDsp.resize(MPISIZE);                                // Resize vector of body recv displacements
    sendCellCnt.resize(MPISIZE);                                // Resize vector of cell send counts
    sendCellDsp.resize(MPISIZE);                                // Resize vector of cell send displacements
    recvCellCnt.resize(MPISIZE);                                // Resize vector of cell recv counts
    recvCellDsp.resize(MPISIZE);                                // Resize vector of cell recv displacements
    MPI_Datatype MPI_TYPE = getType(XMIN[0]);                   // Get MPI data type
    MPI_Allgather(&XMIN[0],3,MPI_TYPE,                          // Gather XMIN
                  &xminAll[0],3,MPI_TYPE,MPI_COMM_WORLD);
    MPI_Allgather(&XMAX[0],3,MPI_TYPE,                          // Gather XMAX
                  &xmaxAll[0],3,MPI_TYPE,MPI_COMM_WORLD);
  }

//! Get neighbor ranks to send to
  void getSendRank(Cells &cells) {
    sendBodyRanks.clear();                                      // Clear send ranks
    sendBodyCellCnt.clear();                                    // Clear send counts
    sendBodyCells.clear();                                      // Clear send body cells
    int oldsize = 0;                                            // Per rank offset of the number of cells to send
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks
      int ic = 0;                                               //  Initialize neighbor dimension counter
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimensions
        if(xminAll[irank][d] < 2 * XMAX[d] - XMIN[d] &&         // If the two domains are touching or overlapping
           xmaxAll[irank][d] > 2 * XMIN[d] - XMAX[d]) {         // in all dimensions, they are neighboring domains
          ic++;                                                 //    Increment neighbor dimension counter
        }                                                       //   Endif for overlapping domains
      }                                                         //  End loop over dimensions
      ic = 3;
      if( ic == 3 && irank != MPIRANK ) {                       //  If ranks are neighbors in all dimensions
        for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {    //   Loop over cells
          if( C->NCHILD == 0 ) {                                //    If cell is a twig
            bool send = false;                                  //     Initialize logical for sending
            if( IMAGES == 0 ) {                                 //     If free boundary condition
              Xperiodic = 0;                                    //      Set periodic coordinate offset
              real R = getDistance(C,xminAll[irank],xmaxAll[irank]);//  Get distance to other domain
              send |= CLET * C->R > THETA * R - EPS2;           //      If the cell seems close enough for P2P
            } else {                                            //     If periodic boundary condition
              for( int ix=-IMAGEDIM[0]; ix<=IMAGEDIM[0]; ++ix ) {//     Loop over x periodic direction
                for( int iy=-IMAGEDIM[1]; iy<=IMAGEDIM[1]; ++iy ) {//    Loop over y periodic direction
                  for( int iz=-IMAGEDIM[2]; iz<=IMAGEDIM[2]; ++iz ) {//   Loop over z periodic direction
                    Xperiodic[0] = ix * 2 * R0;                 //         Coordinate offset for x periodic direction
                    Xperiodic[1] = iy * 2 * R0;                 //         Coordinate offset for y periodic direction
                    Xperiodic[2] = iz * 2 * R0;                 //         Coordinate offset for z periodic direction
                    real R = getDistance(C,xminAll[irank],xmaxAll[irank]);// Get distance to other domain
                    send |= CLET * C->R > THETA * R - EPS2;     //         If the cell seems close enough for P2P
                  }                                             //        End loop over z periodic direction
                }                                               //       End loop over y periodic direction
              }                                                 //      End loop over x periodic direction
            }                                                   //     Endif for periodic boundary condition
            if( send ) {                                        //     If the cell seems close enough for P2P
              sendBodyCells.push_back(C);                       //      Add cell iterator to scells
            }                                                   //     Endif for cell distance
          }                                                     //    Endif for twigs
        }                                                       //   End loop over cells
        sendBodyRanks.push_back(irank);                         //   Add current rank to sendBodyRanks
        sendBodyCellCnt.push_back(sendBodyCells.size()-oldsize);//   Add current cell count to sendBodyCellCnt
        oldsize = sendBodyCells.size();                         //   Set new offset for cell count
      }                                                         //  Endif for neighbor ranks
    }                                                           // End loop over ranks
  }

//! Get size of data to send
  void getSendCount(bool comm=true) {
    int ic = 0, ssize = 0;                                      // Initialize counter and offset for scells
    sendBodyCnt.assign(MPISIZE,0);                              // Initialize send count
    sendBodyDsp.assign(MPISIZE,0);                              // Initialize send displacement
    for( int i=0; i!=int(sendBodyRanks.size()); ++i ) {         // Loop over ranks to send to & recv from
      int irank = sendBodyRanks[i];                             //  Rank to send to & recv from
      for( int c=0; c!=sendBodyCellCnt[i]; ++c,++ic ) {         //  Loop over cells to send to that rank
        C_iter C = sendBodyCells[ic];                           //   Set cell iterator
        for( B_iter B=C->LEAF; B!=C->LEAF+C->NLEAF; ++B ) {     //   Loop over bodies in that cell
          JBody body;                                           //    Set compact body type for sending
          body.ICELL = B->ICELL;                                //    Set cell index of compact body type
          body.X     = B->X;                                    //    Set position of compact body type
          body.SRC   = B->SRC;                                  //    Set source values of compact body type
          body.dxdt  = B->dxdt;
          body.dxdt2 = B->dxdt2;
          body.dgdt  = B->dgdt;
          sendBodies.push_back(body);                           //    Push it into the send buffer
        }                                                       //   End loop over bodies
      }                                                         //  End loop over cells
      sendBodyCnt[irank] = sendBodies.size()-ssize;             //  Set send count of current rank
      sendBodyDsp[irank] = ssize;                               //  Set send displacement of current rank
      ssize += sendBodyCnt[irank];                              //  Increment offset for vector scells
    }                                                           // End loop over ranks
    if( comm ) {                                                // If communication is necessary
      MPI_Alltoall(&sendBodyCnt[0],1,MPI_INT,&recvBodyCnt[0],1,MPI_INT,MPI_COMM_WORLD);// Communicate the send counts
      int rsize = 0;                                            // Initialize total recv count
      for( int i=0; i!=MPISIZE; ++i ) {                         // Loop over ranks to recv from
        recvBodyDsp[i] = rsize;                                 //  Set recv displacements
        rsize += recvBodyCnt[i];                                //  Accumulate recv counts
      }                                                         // End loop over ranks to recv from
      recvBodies.resize(rsize);                                 // Resize recv buffer
    }
  }

//! Get disatnce to other domain
  real getDistance(C_iter C, vect xmin, vect xmax) {
    vect dist;                                                  // Distance vector
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      dist[d] = (C->X[d] + Xperiodic[d] > xmax[d])*             //  Calculate the distance between cell C and
                (C->X[d] + Xperiodic[d] - xmax[d])+             //  the nearest point in domain [xmin,xmax]^3
                (C->X[d] + Xperiodic[d] < xmin[d])*             //  Take the differnece from xmin or xmax
                (C->X[d] + Xperiodic[d] - xmin[d]);             //  or 0 if between xmin and xmax
    }                                                           // End loop over dimensions
    real R = std::sqrt(norm(dist));                             // Scalar distance
    return R;
  }

//! Determine which cells to send
  void getLET(C_iter C0, C_iter C, vect xmin, vect xmax) {
    int level = int(log(MPISIZE-1) / M_LN2 / 3) + 1;            // Level of local root cell
    if( MPISIZE == 1 ) level = 0;                               // Account for serial case
    for( int i=0; i!=C->NCHILD; i++ ) {                         // Loop over child cells
      C_iter CC = C0+C->CHILD[i];                               //  Iterator for child cell
      bool divide = false;                                      //  Initialize logical for dividing
      if( IMAGES == 0 ) {                                       //  If free boundary condition
        Xperiodic = 0;                                          //   Set periodic coordinate offset
        real R = getDistance(CC,xmin,xmax);                     //   Get distance to other domain
        divide |= CLET * CC->R > THETA * R - EPS2;              //   If the cell seems too close and not twig
      } else {                                                  //  If periodic boundary condition
        for( int ix=-IMAGEDIM[0]; ix<=IMAGEDIM[0]; ++ix ) {     //   Loop over x periodic direction
          for( int iy=-IMAGEDIM[1]; iy<=IMAGEDIM[1]; ++iy ) {   //    Loop over y periodic direction
            for( int iz=-IMAGEDIM[2]; iz<=IMAGEDIM[2]; ++iz ) { //     Loop over z periodic direction
              Xperiodic[0] = ix * 2 * R0;                       //      Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //      Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //      Coordinate offset for z periodic direction
              real R = getDistance(CC,xmin,xmax);               //      Get distance to other domain
              divide |= CLET * CC->R > THETA * R - EPS2;        //      If the cell seems too close and not twig
            }                                                   //     End loop over z periodic direction
          }                                                     //    End loop over y periodic direction
        }                                                       //   End loop over x periodic direction
      }                                                         //  Endif for periodic boundary condition
      divide |= R0 / (1 << level) + 1e-5 < CC->R;               //  If the cell is larger than the local root cell
      if( divide && CC->NCHILD != 0 ) {                         //  If the cell seems too close and not twig
        getLET(C0,CC,xmin,xmax);                                //   Traverse the tree further
      } else {                                                  //  If the cell is far or a twig
        assert( R0 / (1 << level) + 1e-5 > CC->R );             //   Can't send cells that are larger than local root
        JCell cell;                                             //   Set compact cell type for sending
        cell.ICELL = CC->ICELL;                                 //   Set index of compact cell type
        cell.M     = CC->M;                                     //   Set Multipoles of compact cell type
        sendCells.push_back(cell);                              //    Push cell into send buffer vector
      }                                                         //  Endif for interaction
    }                                                           // End loop over child cells
    if( C->ICELL == 0 && C->NCHILD == 0 ) {                     // If the root cell has no children
      JCell cell;                                               //  Set compact cell type for sending
      cell.ICELL = C->ICELL;                                    //  Set index of compact cell type
      cell.M     = C->M;                                        //  Set Multipoles of compact cell type
      sendCells.push_back(cell);                                //  Push cell into send buffer vector
    }                                                           // Endif for root cells children
  }

//! Turn recv bodies to twigs
  void rbodies2twigs(Bodies &bodies, Cells &twigs) {
    startTimer("Recv bodies  ");                                //  Start timer
    for( JB_iter JB=recvBodies.begin(); JB!=recvBodies.end(); ++JB ) {// Loop over recv bodies
      Body body;                                                //  Body structure
      body.ICELL = JB->ICELL;                                   //  Set index of body
      body.X     = JB->X;                                       //  Set position of body
      body.SRC   = JB->SRC;                                     //  Set source values of body
      body.dxdt  = JB->dxdt;
      body.dxdt2 = JB->dxdt2;
      body.dgdt  = JB->dgdt;
      bodies.push_back(body);                                   //  Push body into bodies vector
    }                                                           // End loop over recv bodies
    buffer.resize(bodies.size());                               // Resize sort buffer
    stopTimer("Recv bodies  ",printNow);                        //  Stop timer 
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs
  }

//! Turn cells to twigs
  void cells2twigs(Cells &cells, Cells &twigs, bool last) {
    while( !cells.empty() ) {                                   // While cell vector is not empty
      if( cells.back().NCHILD == 0 ) {                          //  If cell has no child
        if( cells.back().NLEAF == 0 || !last ) {                //   If cell has no leaf or is not last iteration
          cells.back().NLEAF = 0;                               //    Set number of leafs to 0
          twigs.push_back(cells.back());                        //    Push cell into twig vector
        }                                                       //   Endif for no leaf
      }                                                         //  Endif for no child
      cells.pop_back();                                         //  Pop last element from cell vector
    }                                                           // End while for cell vector
  }

//! Turn send buffer to twigs
  void send2twigs(Bodies &bodies, Cells &twigs, int offTwigs) {
    for( JC_iter JC=sendCells.begin(); JC!=sendCells.begin()+offTwigs; ++JC ) {// Loop over send buffer
      Cell cell;                                                //  Cell structure
      cell.ICELL = JC->ICELL;                                   //  Set index of cell
      cell.M     = JC->M;                                       //  Set multipole of cell
      cell.NLEAF = cell.NCHILD = 0;                             //  Set number of leafs and children
      cell.LEAF  = bodies.end();                                //  Set pointer to first leaf
      getCenter(cell);                                          //  Set center and radius
      twigs.push_back(cell);                                    //  Push cell into twig vector
    }                                                           // End loop over send buffer
    sendCells.clear();                                          // Clear send buffer
  }

//! Turn recv buffer to twigs
  void recv2twigs(Bodies &bodies, Cells &twigs) {
    for( JC_iter JC=recvCells.begin(); JC!=recvCells.end(); ++JC ) {// Loop over recv buffer
      Cell cell;                                                //  Cell structure
      cell.ICELL = JC->ICELL;                                   //  Set index of cell
      cell.M     = JC->M;                                       //  Set multipole of cell
      cell.NLEAF = cell.NCHILD = 0;                             //  Set number of leafs and children
      cell.LEAF  = bodies.end();                                //  Set pointer to first leaf
      getCenter(cell);                                          //  Set center and radius
      twigs.push_back(cell);                                    //  Push cell into twig vector
    }                                                           // End loop over recv buffer
  }

//! Zip two groups of twigs that overlap
  void zipTwigs(Cells &twigs, Cells &cells, Cells &sticks, bool last) {
    startTimer("Sort resize  ");                                // Start timer
    Cells cbuffer = twigs;                                      // Sort buffer for cells
    stopTimer("Sort resize  ",printNow);                        // Stop timer 
    sortCells(twigs,cbuffer);                                   // Sort twigs in ascending order
    startTimer("Ziptwigs     ");                                // Start timer
    bigint index = -1;                                          // Initialize index counter
    while( !twigs.empty() ) {                                   // While twig vector is not empty
      if( twigs.back().ICELL != index ) {                       //  If twig's index is different from previous
        cells.push_back(twigs.back());                          //   Push twig into cell vector
        index = twigs.back().ICELL;                             //   Update index counter
      } else if ( twigs.back().NLEAF == 0 || !last ) {          //  Elseif twig-twig collision
        cells.back().M += twigs.back().M;                       //   Accumulate the multipole
      } else if ( cells.back().NLEAF == 0 ) {                   //  Elseif twig-body collision
        coef M;                                                 //   Multipole for temporary storage
        M = cells.back().M;                                     //   Save multipoles from cells
        cells.back() = twigs.back();                            //   Copy twigs to cells
        cells.back().M = M;                                     //   Copy back multipoles to cells
        twigs.back().M = M - twigs.back().M;                    //   Take the difference of the two
        if( std::abs(twigs.back().M[0]/M[0]) > 1e-6 ) {         //   If the difference is non-zero
          sticks.push_back(twigs.back());                       //    Save this difference in the sticks vector
        }                                                       //   Endif for non-zero difference
      } else {                                                  //  Else body-body collision (don't do anything)
      }                                                         //  Endif for collision type
      twigs.pop_back();                                         //  Pop last element from twig vector
    }                                                           // End while for twig vector
    stopTimer("Ziptwigs     ",printNow);                        // Stop timer 
    sortCells(cells,cbuffer);                                   // Sort cells in ascending order
    startTimer("Ziptwigs     ");                                // Start timer
    twigs = cells;                                              // Copy cells to twigs
    cells.clear();                                              // Clear cells
    stopTimer("Ziptwigs     ",printNow);                        // Stop timer 
  }

//! Re-index bodies
  void reindexBodies(Bodies &bodies, Cells &twigs, Cells &cells ,Cells &sticks) {
    startTimer("Reindex      ");                                // Start timer
    while( !twigs.empty() ) {                                   // While twig vector is not empty
      if( twigs.back().NLEAF == 0 ) {                           //  If twig has no leafs
        cells.push_back(twigs.back());                          //   Push twig into cell vector
      }                                                         //  Endif for no leafs
      twigs.pop_back();                                         //  Pop last element from twig vector
    }                                                           // End while for twig vector
//    BottomUp::setIndex(bodies,-1,0,0,true);                     // Set index of bodies
    buffer.resize(bodies.size());                               // Resize sort buffer
    stopTimer("Reindex      ",printNow);                        // Stop timer 
//    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order
//    BottomUp::grow(bodies);                                     // Grow tree structure
    sortBodies(bodies,buffer,false);                              // Sort bodies in descending order
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs
    startTimer("Reindex      ");                                // Start timer
    for( C_iter C=twigs.begin(); C!=twigs.end(); ++C ) {        // Loop over cells
      if( sticks.size() > 0 ) {                                 //  If stick vector is not empty
        if( C->ICELL == sticks.back().ICELL ) {                 //   If twig's index is equal to stick's index
          C->M += sticks.back().M;                              //    Accumulate multipole
          sticks.pop_back();                                    //    Pop last element from stick vector
        }                                                       //   Endif for twig's index
      }                                                         //  Endif for stick vector
    }                                                           // End loop over cells
    cells.insert(cells.begin(),twigs.begin(),twigs.end());      // Add twigs to the end of cell vector
    cells.insert(cells.begin(),sticks.begin(),sticks.end());    // Add remaining sticks to the end of cell vector
    sticks.clear();                                             // Clear sticks
    Cells cbuffer = cells;                                      // Sort buffer for cells
    stopTimer("Reindex      ",printNow);                        // Stop timer 
    sortCells(cells,cbuffer);                                   // Sort cells in ascending order
    startTimer("Reindex      ");                                // Start timer
    twigs = cells;                                              // Copy cells to twigs
    cells.clear();                                              // Clear cells
    stopTimer("Reindex      ",printNow);                        // Stop timer 
  }

//! Validate number of send cells
  void checkNumCells(int l) {                                   // Only works with octsection
    int maxLevel = int(log(MPISIZE-1) / M_LN2 / 3) + 1;
    if( MPISIZE == 1 ) maxLevel = 0;
    int octant0 = -1;
    int numCells = 0;
    for( JC_iter JC=sendCells.begin(); JC!=sendCells.end(); ++JC ) {
      int level = getLevel(JC->ICELL);
      int index = JC->ICELL - ((1 << 3*level) - 1) / 7;
      int octant = index / (1 << 3 * (level - maxLevel));
      if( octant != octant0 ) {
        octant0 = octant;
        numCells++;
      }
    }
    int numCellsExpect = (1 << (3 * maxLevel - 1)) / (1 << l);  // Isn't true for far domains
    if( numCellsExpect != numCells && MPIRANK == 0) std::cout << numCells << " " << numCellsExpect << std::endl;
  }

//! Check total charge
  void checkSumMass(Cells &cells) {
    double localMass = 0;
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      if( C->NCHILD == 0 ) {
        localMass += C->M[0].real();
      }
    }
    double globalMass;
    MPI_Allreduce(&localMass,&globalMass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    print("localMass : ",0);
    print(localMass);
    print("globalMass : ",0);
    print(globalMass,0);
    print("\n",0);
  }

public:
//! Constructor
  ParallelFMM() : Partition() {}
//! Destructor
  ~ParallelFMM() {}

//! Set bodies to communicate
  void setCommBodies(Cells &cells) {
    startTimer("Gather bounds");                                // Start timer
    gatherBounds();                                             // Gather bounds of other domain
    stopTimer("Gather bounds",printNow);                        // Stop timer 
    startTimer("Get send rank");                                // Start timer
    getSendRank(cells);                                         // Get neighbor ranks to send to
    stopTimer("Get send rank",printNow);                        // Stop timer 
  }

//! Update bodies using the previous send count
  void updateBodies(bool comm=true) {
    startTimer("Get send cnt ");                                // Start timer
    getSendCount(comm);                                         // Get size of data to send
    stopTimer("Get send cnt ",printNow);                        // Stop timer 
    startTimer("Alltoall B   ");                                // Start timer
    int bytes = sizeof(sendBodies[0]);                          // Byte size of jbody structure
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      sendBodyCnt[i] *= bytes;                                  //  Multiply by bytes
      sendBodyDsp[i] *= bytes;                                  //  Multiply by bytes
      recvBodyCnt[i] *= bytes;                                  //  Multiply by bytes
      recvBodyDsp[i] *= bytes;                                  //  Multiply by bytes
    }                                                           // End loop over ranks
    MPI_Alltoallv(&sendBodies[0],&sendBodyCnt[0],&sendBodyDsp[0],MPI_BYTE,
                  &recvBodies[0],&recvBodyCnt[0],&recvBodyDsp[0],MPI_BYTE,MPI_COMM_WORLD);
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      sendBodyCnt[i] /= bytes;                                  //  Divide by bytes
      sendBodyDsp[i] /= bytes;                                  //  Divide by bytes
      recvBodyCnt[i] /= bytes;                                  //  Divide by bytes
      recvBodyDsp[i] /= bytes;                                  //  Divide by bytes
    }                                                           // End loop over ranks
    sendBodies.clear();                                         // Clear send buffer for bodies
    stopTimer("Alltoall B   ",printNow);                        // Stop timer 
  }

//! Communicate bodies in the local essential tree
  void commBodies(Cells &cells) {
    setCommBodies(cells);                                       // Set bodies to communicate
    updateBodies();                                             // Update bodies with alltoall
  }

//! Convert recvBodies to cells
  void bodies2cells(Bodies &bodies, Cells &cells) {
    Cells twigs,sticks;                                         // Twigs and sticks are special types of cells
    rbodies2twigs(bodies,twigs);                                // Put recv bodies into twig vector
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
  }

//! Communicate cells in the local essential tree
  void commCells(Bodies &bodies, Cells &cells) {
    vect xmin = 0, xmax = 0;                                    // Initialize domain boundaries
    Cells twigs,sticks;                                         // Twigs and sticks are special types of cells
    startTimer("Get LET      ");                                // Start timer
    int ssize = 0;                                              // Initialize offset for send cells
    sendCellCnt.assign(MPISIZE,0);                              // Initialize cell send count
    sendCellDsp.assign(MPISIZE,0);                              // Initialize cell send displacement
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks to send to
      getLET(cells.begin(),cells.end()-1,xminAll[irank],xmaxAll[irank]);//  Determine which cells to send
      sendCellCnt[irank] = sendCells.size()-ssize;              //  Set cell send count of current rank
      sendCellDsp[irank] = ssize;                               //  Set cell send displacement of current rank
      ssize += sendCellCnt[irank];                              //  Increment offset for vector send cells
    }                                                           // End loop over ranks
    stopTimer("Get LET      ",printNow);                        // Stop timer 
    startTimer("Alltoall C   ");                                // Start timer
    MPI_Alltoall(&sendCellCnt[0],1,MPI_INT,&recvCellCnt[0],1,MPI_INT,MPI_COMM_WORLD);// Communicate the send counts
    int rsize = 0;                                              // Initialize total recv count
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks to recv from
      recvCellDsp[irank] = rsize;                               //  Set recv displacements
      rsize += recvCellCnt[irank];                              //  Accumulate recv counts
    }                                                           // End loop over ranks to recv from
    recvCells.resize(rsize);                                    // Resize recv buffer
    int bytes = sizeof(sendCells[0]);                           // Byte size of jbody structure
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      sendCellCnt[i] *= bytes;                                  //  Multiply by bytes
      sendCellDsp[i] *= bytes;                                  //  Multiply by bytes
      recvCellCnt[i] *= bytes;                                  //  Multiply by bytes
      recvCellDsp[i] *= bytes;                                  //  Multiply by bytes
    }                                                           // End loop over ranks
    MPI_Alltoallv(&sendCells[0],&sendCellCnt[0],&sendCellDsp[0],MPI_BYTE,
                  &recvCells[0],&recvCellCnt[0],&recvCellDsp[0],MPI_BYTE,MPI_COMM_WORLD);
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      sendCellCnt[i] /= bytes;                                  //  Divide by bytes
      sendCellDsp[i] /= bytes;                                  //  Divide by bytes
      recvCellCnt[i] /= bytes;                                  //  Divide by bytes
      recvCellDsp[i] /= bytes;                                  //  Divide by bytes
    }                                                           // End loop over ranks
    stopTimer("Alltoall C   ",printNow);                        // Stop timer 
    rbodies2twigs(bodies,twigs);                                // Put recv bodies into twig vector
    startTimer("Cells2twigs  ");                                // Start timer
    cells2twigs(cells,twigs,true);                              // Put cells into twig vector
    stopTimer("Cells2twigs  ",printNow);                        // Stop timer 
    startTimer("Recv2twigs   ");                                // Start timer
    recv2twigs(bodies,twigs);                                   // Put recv buffer into twig vector
    stopTimer("Recv2twigs   ",printNow);                        // Stop timer 
    zipTwigs(twigs,cells,sticks,true);                          // Zip two groups of twigs that overlap
    reindexBodies(bodies,twigs,cells,sticks);                   // Re-index bodies
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
    sendCells.clear();                                          // Clear send buffer
    recvCells.clear();                                          // Clear recv buffer
#ifdef DEBUG
    print("M[0] @ root   : ",0);                                // Print identifier
    print((cells.end()-1)->M[0]);                               // Print monopole of root (should be 1 for test)
    print("bodies.size() : ",0);                                // Print identifier
    print(bodies.size());                                       // Print size of body vector
#endif
    sendCells.clear();                                          // Clear send buffer
    recvCells.clear();                                          // Clear recv buffer
  }

//! Remove cells that belong to current process
  void eraseLocalTree(Cells &cells) {
    int level = int(log(MPISIZE-1) / M_LN2 / 3) + 1;            // Level of process root cell
    if( MPISIZE == 1 ) level = 0;                               // Account for serial case
    int off = ((1 << 3 * level) - 1) / 7;                       // Levelwise offset of ICELL
    int size = (1 << 3 * level) / MPISIZE;                      // Number of cells to remove
    int begin = MPIRANK * size + off;                           // Begin index of cells to remove
    int end = (MPIRANK + 1) * size + off;                       // End index of cells to remove
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {        // Loop over cells
      if( begin <= C->ICELL && C->ICELL < end ) {               //  If cell is within the removal range
        C_iter CP = cells.begin()+C->PARENT;                    //   Iterator of parent cell
        int ic = 0;                                             //   Initialize child cell counter
        for( int c=0; c!=CP->NCHILD; ++c ) {                    //   Loop over child cells of parent
          C_iter CC = cells.begin()+CP->CHILD[c];               //    Iterator of child cell of parent
          if( CC->ICELL != C->ICELL ) {                         //    If child cell of parent is not current cell
            CP->CHILD[ic] = CP->CHILD[c];                       //     Copy pointer to child cell
            ic++;                                               //     Increment child cell counter
          }                                                     //    Endif for current child cell
        }                                                       //   End loop over child cells of parent
        CP->NCHILD--;                                           //   Decrement number of child cells of parent
      }                                                         //  Endif for removal range
    }                                                           // End loop over cells
  }
};

#endif
