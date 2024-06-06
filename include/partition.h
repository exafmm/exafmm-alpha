#ifndef partition_h
#define partition_h
#include "mympi.h"
#include "serialfmm.h"

//! Handles all the partitioning of domains
class Partition : public MyMPI, public SerialFMM {
private:
  int numCells1D;                                               //!< Number of cells in one dimension (leaf level)

protected:
  vect XMIN;                                                    //!< Minimum position vector of bodies
  vect XMAX;                                                    //!< Maximum position vector of bodies

public:
//! Constructor
  Partition() : SerialFMM() {}
//! Destructor
  ~Partition() {}

//! Set bounds of domain to be partitioned
  void setGlobDomain(Bodies &bodies, vect x0=0, real r0=M_PI) {
    numCells1D = 1 << getMaxLevel(bodies);                      // Set initial number of bodies
    B_iter B = bodies.begin();                                  // Reset body iterator
    XMIN = XMAX = B->X;                                         // Initialize xmin,xmax
    MPI_Datatype MPI_TYPE = getType(XMIN[0]);                   // Get MPI data type
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->X[d] < XMIN[d]) XMIN[d] = B->X[d];           //   Determine xmin
        else if(B->X[d] > XMAX[d]) XMAX[d] = B->X[d];           //   Determine xmax
      }                                                         //  End loop over each dimension
    }                                                           // End loop over bodies
    vect X;                                                     // Recv buffer
    MPI_Allreduce(XMAX,X,3,MPI_TYPE,MPI_MAX,MPI_COMM_WORLD);    // Reduce global maximum
    XMAX = X;                                                   // Get data from buffer
    MPI_Allreduce(XMIN,X,3,MPI_TYPE,MPI_MIN,MPI_COMM_WORLD);    // Reduce global minimum
    XMIN = X;                                                   // Get data from buffer
    for( int d=0; d!=3; ++d ) {                                 // Loop over each dimension
      X0[d] = (XMAX[d] + XMIN[d]) / 2;                          //  Calculate center of domain
      X0[d] = int(X0[d]+.5);                                    //  Shift center to nearest integer
      R0 = std::max(XMAX[d] - X0[d], R0);                       //  Calculate max distance from center
      R0 = std::max(X0[d] - XMIN[d], R0);                       //  Calculate max distance from center
    }                                                           // End loop over each dimension
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      if( X0[0]-R0 < x0[0]-r0 || x0[0]+r0 < X0[0]+R0            //  Check for outliers in x direction
       || X0[1]-R0 < x0[1]-r0 || x0[1]+r0 < X0[1]+R0            //  Check for outliers in y direction
       || X0[2]-R0 < x0[2]-r0 || x0[2]+r0 < X0[2]+R0 ) {        //  Check for outliers in z direction
        std::cout << "Error: Particles located outside periodic domain @ rank "
                  << MPIRANK << std::endl;                      //   Print error message
      }                                                         //  End if for outlier checking
      X0 = x0;                                                  //  Center is [0, 0, 0]
      R0 = r0;                                                  //  Radius is M_PI
    } else {                                                    // If not periodic boundary condition
      R0 += 1e-5;                                               //  Add some leeway to root radius
    }                                                           // Endif for periodic boundary condition
    XMAX = X0 + R0;                                             // Reposition global maximum
    XMIN = X0 - R0;                                             // Reposition global minimum
  }

//! Send bodies to next rank (round robin)
  void shiftBodies(Bodies &bodies) {
    int newSize;                                                // New number of bodies
    int oldSize = bodies.size();                                // Current number of bodies
    const int bytes = sizeof(bodies[0]);                        // Byte size of body structure
    const int isend = (MPIRANK + 1          ) % MPISIZE;        // Send to next rank (wrap around)
    const int irecv = (MPIRANK - 1 + MPISIZE) % MPISIZE;        // Receive from previous rank (wrap around)
    MPI_Request sreq,rreq;                                      // Send, recv request handles

    MPI_Isend(&oldSize,1,MPI_INT,irecv,0,MPI_COMM_WORLD,&sreq); // Send current number of bodies
    MPI_Irecv(&newSize,1,MPI_INT,isend,0,MPI_COMM_WORLD,&rreq); // Receive new number of bodies
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for recv to complete

    buffer.resize(newSize);                                     // Resize buffer to new number of bodies
    MPI_Isend(&bodies[0],oldSize*bytes,MPI_BYTE,irecv,          // Send bodies to next rank
              1,MPI_COMM_WORLD,&sreq);
    MPI_Irecv(&buffer[0],newSize*bytes,MPI_BYTE,isend,          // Receive bodies from previous rank
              1,MPI_COMM_WORLD,&rreq);
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for recv to complete
    bodies = buffer;                                            // Copy bodies from buffer
  }

//! Partition by recursive octsection
  void octsection(Bodies &bodies) {
    startTimer("Partition    ");                                // Start timer
    int byte = sizeof(bodies[0]);                               // Byte size of body structure
    int level = int(log(MPISIZE-1) / M_LN2 / 3) + 1;            // Level of local root cell
    if( MPISIZE == 1 ) level = 0;                               // For serial execution local root cell is root cell
    BottomUp::setIndex(bodies,level);                           // Set index of bodies for that level
    buffer.resize(bodies.size());                               // Resize sort buffer
    stopTimer("Partition    ");                                 // Stop timer 
    sortBodies(bodies,buffer);                                  // Sort bodies in ascending order
    startTimer("Partition    ");                                // Start timer
    int *scnt = new int [MPISIZE];                              // Send count
    int *sdsp = new int [MPISIZE];                              // Send displacement
    int *rcnt = new int [MPISIZE];                              // Recv count
    int *rdsp = new int [MPISIZE];                              // Recv displacement
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      scnt[i] = 0;                                              //  Initialize send counts
    }                                                           // End loop over ranks
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int index = B->ICELL - ((1 << 3*level) - 1) / 7;          //  Get levelwise index
      int irank = index / (int(pow(8,level)) / MPISIZE);        //  Get rank which the cell belongs to
      scnt[irank]++;                                            //  Fill send count bucket
    }                                                           // End loop over bodies
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD); // Communicate send count to get recv count
    sdsp[0] = rdsp[0] = 0;                                      // Initialize send/recv displacements
    for( int irank=0; irank!=MPISIZE-1; ++irank ) {             // Loop over ranks
      sdsp[irank+1] = sdsp[irank] + scnt[irank];                //  Set send displacement based on send count
      rdsp[irank+1] = rdsp[irank] + rcnt[irank];                //  Set recv displacement based on recv count
    }                                                           // End loop over ranks
    buffer.resize(rdsp[MPISIZE-1]+rcnt[MPISIZE-1]);             // Resize recv buffer
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks
      scnt[irank] *= byte;                                      //  Multiply send count by byte size of data
      sdsp[irank] *= byte;                                      //  Multiply send displacement by byte size of data
      rcnt[irank] *= byte;                                      //  Multiply recv count by byte size of data
      rdsp[irank] *= byte;                                      //  Multiply recv displacement by byte size of data
    }                                                           // End loop over ranks
    MPI_Alltoallv(&bodies[0],scnt,sdsp,MPI_BYTE,&buffer[0],rcnt,rdsp,MPI_BYTE,MPI_COMM_WORLD);// Communicat bodies
    bodies = buffer;                                            // Copy recv buffer to bodies
    delete[] scnt;                                              // Delete send count
    delete[] sdsp;                                              // Delete send displacement
    delete[] rcnt;                                              // Delete recv count
    delete[] rdsp;                                              // Delete recv displacement
    stopTimer("Partition    ",printNow);                        // Stop timer 
  }

//! Send bodies back to where they came from
  void unpartition(Bodies &bodies) {
    startTimer("Unpartition  ");                                // Start timer
    int byte = sizeof(bodies[0]);                               // Byte size of body structure
    int *scnt = new int [MPISIZE];                              // Send count
    int *sdsp = new int [MPISIZE];                              // Send displacement
    int *rcnt = new int [MPISIZE];                              // Recv count
    int *rdsp = new int [MPISIZE];                              // Recv displacement
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->ICELL = B->IPROC;                                      //  Copy process rank to cell index for sorting
    }                                                           // End loop over bodies
    buffer.resize(bodies.size());                               // Resize sort buffer
    stopTimer("Unpartition  ");                                 // Stop timer 
    sortBodies(bodies,buffer);                                  // Sort bodies in ascending order
    startTimer("Unpartition  ");                                // Start timer
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      scnt[i] = 0;                                              //  Initialize send counts
    }                                                           // End loop over ranks
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int irank = B->IPROC;                                     //  Get rank which the body belongs to
      scnt[irank]++;                                            //  Fill send count bucket
    }                                                           // End loop over bodies
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD); // Communicate send count to get recv count
    sdsp[0] = rdsp[0] = 0;                                      // Initialize send/recv displacements
    for( int irank=0; irank!=MPISIZE-1; ++irank ) {             // Loop over ranks
      sdsp[irank+1] = sdsp[irank] + scnt[irank];                //  Set send displacement based on send count
      rdsp[irank+1] = rdsp[irank] + rcnt[irank];                //  Set recv displacement based on recv count
    }                                                           // End loop over ranks
    buffer.resize(rdsp[MPISIZE-1]+rcnt[MPISIZE-1]);             // Resize recv buffer
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks
      scnt[irank] *= byte;                                      //  Multiply send count by byte size of data
      sdsp[irank] *= byte;                                      //  Multiply send displacement by byte size of data
      rcnt[irank] *= byte;                                      //  Multiply recv count by byte size of data
      rdsp[irank] *= byte;                                      //  Multiply recv displacement by byte size of data
    }                                                           // End loop over ranks
    MPI_Alltoallv(&bodies[0],scnt,sdsp,MPI_BYTE,&buffer[0],rcnt,rdsp,MPI_BYTE,MPI_COMM_WORLD);// Communicat bodies
    bodies = buffer;                                            // Copy recv buffer to bodies
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      assert(B->IPROC == MPIRANK);                              //  Make sure bodies are in the right rank
    }                                                           // End loop over bodies
    delete[] scnt;                                              // Delete send count
    delete[] sdsp;                                              // Delete send displacement
    delete[] rcnt;                                              // Delete recv count
    delete[] rdsp;                                              // Delete recv displacement
    stopTimer("Unpartition  ",printNow);                        // Stop timer 
  }
};

#endif
