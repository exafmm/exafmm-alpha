#ifndef evaluator_h
#define evaluator_h
#include "kernel.h"

//! Interface between tree and kernel
class Evaluator : public Kernel {
protected:
  C_iter      CI0;                                              //!< icells.begin()
  C_iter      CIB;                                              //!< icells begin per call
  C_iter      CIE;                                              //!< icells end per call
  C_iter      CJ0;                                              //!< jcells.begin()
  C_iter      CJB;                                              //!< jcells begin per call
  C_iter      CJE;                                              //!< jcells end per call
  Pairs       pairs;                                            //!< Stack of interacting cell pairs
  Lists       listM2L;                                          //!< M2L interaction list
  Lists       listM2P;                                          //!< M2P interaction list
  Lists       listP2P;                                          //!< P2P interaction list
  real        timeM2L;                                          //!< M2L execution time
  real        timeM2P;                                          //!< M2P execution time
  real        timeP2P;                                          //!< P2P execution time

  int         Iperiodic;                                        //!< Periodic image flag (using each bit for images)
  const int   Icenter;                                          //!< Periodic image flag at center
  const int   Iall;                                             //!< Periodic image flag for all neighbors
  Maps        flagM2L;                                          //!< Existance of periodic image for M2L
  Maps        flagM2P;                                          //!< Existance of periodic image for M2P
  Maps        flagP2P;                                          //!< Existance of periodic image for P2P

private:
//! Tree walk for treecode
  void treecode(C_iter Ci, C_iter Cj) {
    if( Ci->NCHILD == 0 && Cj->NCHILD == 0) {                   // If both cells are twigs
      if( Cj->NLEAF != 0 ) {                                    // If the twig has leafs
        testMACP2P(Ci,Cj);                                      //  Test multipole acceptance criteria for P2P kernel
      } else {                                                  // If the twig has no leafs
//#ifdef DEBUG
        std::cout << "Cj->ICELL=" << Cj->ICELL << " has no leaf. Doing M2P instead of P2P." << std::endl;
//#endif
        listM2P[Ci-CI0].push_back(Cj);                          // Push source cell into M2P interaction list
      }                                                         // Endif for twigs with leafs
    } else if ( Ci->NCHILD != 0 ) {                             // If target is not twig
      for( int i=0; i<Ci->NCHILD; i++ ) {                       //  Loop over child cells of target
        testMACM2P(CI0+Ci->CHILD[i],Cj);                        //   Test multipole acceptance criteria for M2P kernel
      }                                                         //  End loop over child cells of target
    } else {                                                    // If target is twig
      for( int i=0; i<Cj->NCHILD; i++ ) {                       //  Loop over child cells of source
        testMACM2P(Ci,CJ0+Cj->CHILD[i]);                        //   Test multipole acceptance criteria for M2P kernel
      }                                                         //  End loop over child cells of source
    }                                                           // Endif for type of interaction
  }

//! Tree walk for FMM
  void FMM(C_iter Ci, C_iter Cj) {
    if( Ci->NCHILD == 0 && Cj->NCHILD == 0 ) {                  // If both cells are twigs
      if( Cj->NLEAF != 0 ) {                                    // If the twig has leafs
        testMACP2P(Ci,Cj);                                      //  Test multipole acceptance criteria for P2P kernel
      } else {                                                  // If the twig has no leafs
//#ifdef DEBUG
        std::cout << "Cj->ICELL=" << Cj->ICELL << " has no leaf. Doing M2P instead of P2P." << std::endl;
//#endif
        listM2P[Ci-CI0].push_back(Cj);                          // Push source cell into M2P interaction list
      }                                                         // Endif for twigs with leafs
    } else if ( Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->R > Cj->R) ) {// If source is twig or target is larger
      for( int i=0; i<Ci->NCHILD; i++ ) {                       //  Loop over child cells of target
        testMACM2L(CI0+Ci->CHILD[i],Cj);                        //   Test multipole acceptance criteria for M2L kernel
      }                                                         //  End loop over child cells of target
    } else {                                                    // If target is twig or source is larger
      for( int i=0; i<Cj->NCHILD; i++ ) {                       //  Loop over child cells of source
        testMACM2L(Ci,CJ0+Cj->CHILD[i]);                        //   Test multipole acceptance criteria for M2L kernel
      }                                                         //  End loop over child cells of source
    }                                                           // Endif for type of interaction
  }

//! Tree walk for treecode-FMM hybrid
  void hybrid(C_iter Ci, C_iter Cj) {
    if( Ci->NCHILD == 0 && Cj->NCHILD == 0 ) {                  // If both cells are twigs
      if( Cj->NLEAF != 0 ) {                                    // If the twig has leafs
        testMACP2P(Ci,Cj);                                      //  Test multipole acceptance criteria for P2P kernel
      } else {                                                  // If the twig has no leafs
//#ifdef DEBUG
        std::cout << "Cj->ICELL=" << Cj->ICELL << " has no leaf. Doing M2P instead of P2P." << std::endl;
//#endif
        listM2P[Ci-CI0].push_back(Cj);                          // Push source cell into M2P interaction list
      }                                                         // Endif for twigs with leafs
    } else if ( Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->R > Cj->R) ) {// If source is twig or target is larger
      for( int i=0; i<Ci->NCHILD; i++ ) {                       //  Loop over child cells of target
        int Ni = (CI0+Ci->CHILD[i])->NLEAF;                     //   Number of target leafs
        int Nj = Cj->NLEAF;                                     //   Number of source leafs
        if( timeP2P*Nj < timeM2P && timeP2P*Ni*Nj < timeM2L ) { //   If P2P is fastest
          testMACP2P(CI0+Ci->CHILD[i],Cj);                      //    Test multipole acceptance criteria for P2P kernel
        } else if ( timeM2P < timeP2P*Nj && timeM2P*Ni < timeM2L ) {// If M2P is fastest
          testMACM2P(CI0+Ci->CHILD[i],Cj);                      //    Test multipole acceptance criteria for M2P kernel
        } else {                                                //   If M2L is fastest
          testMACM2L(CI0+Ci->CHILD[i],Cj);                      //    Test multipole acceptance criteria for M2L kernel
        }                                                       //   End if for fastest kernel
      }                                                         //  End loop over child cells of target
    } else {                                                    // If target is twig or source is larger
      for( int i=0; i<Cj->NCHILD; i++ ) {                       //  Loop over child cells of source
        int Ni = Ci->NLEAF;                                     //   Number of target leafs
        int Nj = (CJ0+Cj->CHILD[i])->NLEAF;                     //   Number of source leafs
        if( timeP2P*Nj < timeM2P && timeP2P*Ni*Nj < timeM2L ) { //   If P2P is fastest
          testMACP2P(Ci,CJ0+Cj->CHILD[i]);                      //    Test multipole acceptance criteria for P2P kernel
        } else if ( timeM2P < timeP2P*Nj && timeM2P*Ni < timeM2L ) {// If M2P is fastest
          testMACM2P(Ci,CJ0+Cj->CHILD[i]);                      //    Test multipole acceptance criteria for M2P kernel
        } else {                                                //   If M2L is fastest
          testMACM2L(Ci,CJ0+Cj->CHILD[i]);                      //    Test multipole acceptance criteria for M2L kernel
        }                                                       //   End if for fastest kernel
      }                                                         //  End loop over child cells of source
    }                                                           // Endif for type of interaction
  }

protected:
  void timeKernels();                                           //!< Time all kernels for auto-tuning

public:
//! Constructor
  Evaluator() : Icenter(8192), Iall(134217727) {}
//! Destructor
  ~Evaluator() {}

//! Set kernel name
  void setKernel(std::string name) {
    if( name == "Laplace" ) {                                   //  If Laplace
      kernelName = Laplace;                                     //   Set kernel name to Laplace
    } else if( name == "BiotSavart" ) {                         //  If BiotSavart
      kernelName = BiotSavart;                                  //   Set kernel name to BiotSavart
    } else if( name == "Stretching" ) {                         //  If Stretching
      kernelName = Stretching;                                  //   Set kernel name to Stretching
    } else if( name == "Gaussian" ) {                           //  If Gaussian
      kernelName = Gaussian;                                    //   Set kernel name to Gaussian
    } else if( name == "CoulombVdW" ) {                         //  If CoulombVdW
      kernelName = CoulombVdW;                                  //   Set kernel name to CoulombVdW
    }
  }

//! Traverse tree to get interaction list
  void traverse(Cells &cells, Cells &jcells, int method) {
    C_iter root = cells.end() - 1;                              // Iterator for root target cell
    C_iter jroot = jcells.end() - 1;                            // Iterator for root source cell
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      int numNeighbors = 1;                                     //  Number of neighbors
      for (int d=0; d<3; d++)                                   //  Loop over dimensions
        numNeighbors *= IMAGEDIM[d] == 0 ? 1 : 3;               //   Multiply by 3 for each dimension
      numNeighbors--;                                           //  Subtract self
      jroot = jcells.end() - 1 - numNeighbors * (IMAGES - 1);   //  The root is not at the end
    }                                                           // Endif for periodic boundary condition
    CI0 = cells.begin();                                        // Set begin iterator for target cells
    CJ0 = jcells.begin();                                       // Set begin iterator for source cells
    listM2L.resize(cells.size());                               // Resize M2L interaction list
    listM2P.resize(cells.size());                               // Resize M2P interaction list
    listP2P.resize(cells.size());                               // Resize P2P interaction list
    flagM2L.resize(cells.size());                               // Resize M2L periodic image flag
    flagM2P.resize(cells.size());                               // Resize M2P periodic image flag
    flagP2P.resize(cells.size());                               // Resize P2P periodic image flag
    if( IMAGES == 0 ) {                                         // If free boundary condition
      Iperiodic = Icenter;                                      //  Set periodic image flag to center
      Xperiodic = 0;                                            //  Set periodic coordinate offset
      Pair pair(root,jroot);                                    //  Form pair of root cells
      pairs.push(pair);                                         //  Push pair to stack
      while( !pairs.empty() ) {                                 //  While interaction stack is not empty
        pair = pairs.top();                                     //   Get interaction pair from top of stack
        pairs.pop();                                            //   Pop interaction stack
        switch (method) {                                       //   Swtich between methods
        case 0 : treecode(pair.first,pair.second); break;       //    0 : treecode
        case 1 : FMM(pair.first,pair.second);      break;       //    1 : FMM
        case 2 : hybrid(pair.first,pair.second);   break;       //    2 : hybrid
        }                                                       //   End switch between methods
      }                                                         //  End while loop for interaction stack
    } else {                                                    // If periodic boundary condition
      int I = 0;                                                //  Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //    Loop over z periodic direction
	    if( !ix|IMAGEDIM[0] && !iy|IMAGEDIM[1] && !iz|IMAGEDIM[2] ) { // If match periodic dimension
              Iperiodic = 1 << I;                               //      Set periodic image flag
              Xperiodic[0] = ix * 2 * R0;                       //      Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //      Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //      Coordinate offset for z periodic direction
              Pair pair(root,jroot);                            //      Form pair of root cells
              pairs.push(pair);                                 //      Push pair to stack
              while( !pairs.empty() ) {                         //      While interaction stack is not empty
                pair = pairs.top();                             //       Get interaction pair from top of stack
                pairs.pop();                                    //       Pop interaction stack
                switch (method) {                               //       Swtich between methods
                case 0 : treecode(pair.first,pair.second); break;//      0 : treecode
                case 1 : FMM(pair.first,pair.second);      break;//      1 : FMM
                case 2 : hybrid(pair.first,pair.second);   break;//      2 : hybrid
                }                                               //       End switch between methods
              }                                                 //      End while loop for interaction stack
	    }                                                   //     Endif for periodic dimension
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {   //  Loop over target cells
        listM2L[Ci-CI0].sort();                                 //   Sort interaction list
        listM2L[Ci-CI0].unique();                               //   Eliminate duplicate periodic entries
        listM2P[Ci-CI0].sort();                                 //   Sort interaction list
        listM2P[Ci-CI0].unique();                               //   Eliminate duplicate periodic entries
        listP2P[Ci-CI0].sort();                                 //   Sort interaction list
        listP2P[Ci-CI0].unique();                               //   Eliminate duplicate periodic entries
      }                                                         //  End loop over target cells
    }                                                           // Endif for periodic boundary condition
  }

//! Upward phase for periodic cells
  void upwardPeriodic(Cells &jcells) {
    Cells pccells, pjcells;                                     // Periodic jcells for M2L/M2P & M2M
    pccells.push_back(jcells.back());                           // Root cell is first periodic cell
    for( int level=0; level<IMAGES-1; ++level ) {               // Loop over sublevels of tree
      C_iter C = pccells.end() - 1;                             //  Set previous periodic cell as source
      for( int ix=-IMAGEDIM[0]; ix<=IMAGEDIM[0]; ++ix ) {       //  Loop over x periodic direction
        for( int iy=-IMAGEDIM[1]; iy<=IMAGEDIM[1]; ++iy ) {     //   Loop over y periodic direction
          for( int iz=-IMAGEDIM[2]; iz<=IMAGEDIM[2]; ++iz ) {   //    Loop over z periodic direction
            Cell cell;                                          //     New periodic jcell for M2M
            cell.X[0] = C->X[0] + ix * 2 * C->R;                //     Set new x coordinate for periodic image
            cell.X[1] = C->X[1] + iy * 2 * C->R;                //     Set new y cooridnate for periodic image
            cell.X[2] = C->X[2] + iz * 2 * C->R;                //     Set new z coordinate for periodic image
            cell.M = C->M;                                      //     Copy multipoles to new periodic image
            pjcells.push_back(cell);                            //     Push cell into periodic jcell vector
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      Cell cell;                                                //  New periodic cell at next sublevel
      cell.X = C->X;                                            //  This is the center cell
      cell.R = 3 * C->R;                                        //  The cell size increase three times
      pccells.push_back(cell);                                  //  Push cell into periodic cell vector
      CI = pccells.end() - 1;                                   //  Set current cell as target for M2M
      while( !pjcells.empty() ) {                               //  While there are periodic jcells remaining
        CJ = pjcells.end() - 1;                                 //   Set current jcell as source for M2M
        selectM2M_CPU();                                        //   Select M2M_CPU kernel
        pjcells.pop_back();                                     //   Pop last element from periodic jcell vector
      }                                                         //  End while for remaining periodic jcells
      for( int ix=-IMAGEDIM[0]; ix<=IMAGEDIM[0]; ++ix ) {       //  Loop over x periodic direction
        for( int iy=-IMAGEDIM[1]; iy<=IMAGEDIM[1]; ++iy ) {     //   Loop over y periodic direction
          for( int iz=-IMAGEDIM[2]; iz<=IMAGEDIM[2]; ++iz ) {   //    Loop over z periodic direction
            if( ix != 0 || iy != 0 || iz != 0 ) {               //     If periodic cell is not at center
              cell.X[0]  = CI->X[0] + ix * 2 * CI->R;           //      Set new x coordinate for periodic image
              cell.X[1]  = CI->X[1] + iy * 2 * CI->R;           //      Set new y cooridnate for periodic image
              cell.X[2]  = CI->X[2] + iz * 2 * CI->R;           //      Set new z coordinate for periodic image
              cell.M     = CI->M;                               //      Copy multipoles to new periodic image
              cell.NLEAF = cell.NCHILD = 0;                     //      Initialize NLEAF & NCHILD
              jcells.push_back(cell);                           //      Push cell into periodic jcell vector
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
    }                                                           // End loop over sublevels of tree
  }

//! Initialize GPU
  void initialize() {
    if( kernelName == Laplace ) {                               // If Laplace kernel
      LaplaceInit();                                            //  Initialize GPU
    } else if ( kernelName == BiotSavart ) {                    // If Biot Savart kernel
      BiotSavartInit();                                         //  Initialize GPU
    } else if ( kernelName == Stretching ) {                    // If Stretching kernel
      StretchingInit();                                         //  Initialize GPU
    } else if ( kernelName == Gaussian ) {                      // If Gaussian kernel
      GaussianInit();                                           //  Initialize GPU
    } else if ( kernelName == CoulombVdW ) {                    // If CoulombVdW kernel
      CoulombVdWInit();                                         //  Initialize GPU
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type in initialize" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }

//! Finalize GPU
  void finalize() {
    if( kernelName == Laplace ) {                               // If Laplace kernel
      LaplaceFinal();                                           //  Finalize GPU
    } else if ( kernelName == BiotSavart ) {                    // If Biot Savart kernel
      BiotSavartFinal();                                        //  Finalize GPU
    } else if ( kernelName == Stretching ) {                    // If Stretching kernel
      StretchingFinal();                                        //  Finalize GPU
    } else if ( kernelName == Gaussian ) {                      // If Gaussian kernel
      GaussianFinal();                                          //  Finalize GPU
    } else if ( kernelName == CoulombVdW ) {                    // If CoulombVdW kernel
      CoulombVdWFinal();                                        //  Finalize GPU
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type in finalize" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }

  void setSourceBody();                                         //!< Set source buffer for bodies
  void setSourceCell(bool isM);                                 //!< Set source buffer for cells
  void setTargetBody(Lists lists, Maps flags);                  //!< Set target buffer for bodies
  void setTargetCell(Lists lists, Maps flags);                  //!< Set target buffer for cells
  void getTargetBody(Lists &lists);                             //!< Get body values from target buffer
  void getTargetCell(Lists &lists, bool isM);                   //!< Get cell values from target buffer
  void clearBuffers();                                          //!< Clear GPU buffers

  void testMACP2P(C_iter Ci, C_iter Cj);                        //!< Test multipole acceptance criteria for P2P kernel
  void testMACM2L(C_iter Ci, C_iter Cj);                        //!< Test multipole acceptance criteria for M2L kernel
  void testMACM2P(C_iter Ci, C_iter Cj);                        //!< Test multipole acceptance criteria for M2P kernel
  void traversePeriodic(Cells &cells, Cells &jcells, int method);//!< Traverse tree for periodic cells
  void evalP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU=false);//!< Evaluate P2P kernel (all pairs)
  void evalP2M(Cells &twigs);                                   //!< Evaluate P2M kernel
  void evalM2M(Cells &cells);                                   //!< Evaluate M2M kernel
  void evalM2L(Cells &cells, bool kernel=false);                //!< Evaluate M2L kernel
  void evalM2P(Cells &cells, bool kernel=false);                //!< Evaluate M2P kernel
  void evalP2P(Cells &cells, bool kernel=false);                //!< Evaluate P2P kernel (near field)
  void evalL2L(Cells &cells);                                   //!< Evaluate L2L kernel
  void evalL2P(Cells &cells);                                   //!< Evaluate L2P kernel
};
#if cpu
#include "../kernel/cpuEvaluator.cxx"
#elif gpu
#include "../kernel/gpuEvaluator.cxx"
#endif

#endif
