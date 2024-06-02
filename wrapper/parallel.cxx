#include "dataset.h"
#include "parallelfmm.h"

int main() {
  const int numBodies = 10000;
  const int numTarget = 100;
  Bodies bodies(numBodies);
  Cells cells;
  ParallelFMM FMM;
  FMM.setKernel("BiotSavart");
  FMM.initialize();
  FMM.IMAGES = 1;
  FMM.THETA = 1 / sqrtf(4);
  bool printNow = (MPIRANK == 0);

  FMM.startTimer("Set bodies   ");
  Dataset Data;
  Data.setKernel("BiotSavart");
  Data.random(bodies,MPIRANK+1);
  FMM.stopTimer("Set bodies   ",printNow);

  FMM.startTimer("Set domain   ");
  FMM.setGlobDomain(bodies);
  FMM.stopTimer("Set domain   ",printNow);

  FMM.startTimer("Partition    ");
  FMM.octsection(bodies);
  FMM.stopTimer("Partition    ",printNow);

  FMM.startTimer("Build tree   ");
  FMM.bottomup(bodies,cells);
  FMM.stopTimer("Build tree   ",printNow);

  FMM.startTimer("Set LET      ");
  FMM.setCommBodies(cells);
  FMM.stopTimer("Set LET      ",printNow);

  FMM.startTimer("Send LET     ");
  Bodies jbodies;
  Cells jcells = cells;
  if( MPISIZE != 1 ) {
#pragma omp parallel sections num_threads(2)
    {
#pragma omp section
      {
        FMM.downward(cells,jcells,1,false);
      }
#pragma omp section
      {
        FMM.updateBodies();
      }
    }
    jbodies = bodies;
    jcells = cells;
    FMM.commCells(jbodies,jcells);
    FMM.eraseLocalTree(jcells);
  }
  FMM.stopTimer("Send LET     ",printNow);

  FMM.startTimer("Traversal    ");
  FMM.downward(cells,jcells,1);
  FMM.stopTimer("Traversal    ",printNow);

  FMM.startTimer("Direct sum   ");
  jbodies = bodies;
  Data.sampleBodies(bodies,numTarget);
  Bodies bodies2 = bodies;
  Data.initTarget(bodies2);
  for (int i=0; i<MPISIZE; i++) {
    FMM.shiftBodies(jbodies);
    FMM.evalP2P(bodies2,jbodies);
    if(printNow) std::cout << "Direct loop   : " << i+1 << "/" << MPISIZE << std::endl;
  }
  FMM.stopTimer("Direct sum   ",printNow);

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
  Data.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  MPI_Datatype MPI_TYPE = FMM.getType(diff1);
  MPI_Reduce(&diff1,&diff3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm1,&norm3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&diff2,&diff4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm2,&norm4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  if(printNow) Data.printError(diff3,norm3,diff4,norm4);
  FMM.finalize();
}
