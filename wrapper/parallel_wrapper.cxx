#include "parallelfmm.h"

ParallelFMM * FMM;

extern "C" void fmm_init_(int & images) {
  FMM = new ParallelFMM;
  FMM->setKernel("BiotSavart");
  FMM->initialize();
  FMM->IMAGES = images;
  FMM->THETA = 1 / sqrtf(4);
}

extern "C" void fmm_finalize_() {
  FMM->finalize();
  delete FMM;
}

extern "C" void fmm_partition_(int & numBodies, double * x, double * g, double * u) {
  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(numBodies);
  B_iter B=bodies.begin();
  for (int b=0; b<numBodies; b++,B++) {
    B->IBODY = b;
    B->IPROC = MPIRANK;
    for (int d=0; d<3; d++) {
      B->X[d] = x[3*b+d];
      B->SRC[d] = g[3*b+d];
    }
    B->SRC[3] = 1e-6;
  }
  FMM->stopTimer("Set bodies   ",printNow);

  FMM->startTimer("Set domain   ");
  FMM->setGlobDomain(bodies);
  FMM->stopTimer("Set domain   ",printNow);

  FMM->startTimer("Partition    ");
  FMM->octsection(bodies);
  FMM->stopTimer("Partition    ",printNow);

  numBodies = bodies.size();
  B=bodies.begin();
  for(int b=0; b<numBodies; b++,B++) {
    for (int d=0; d<3; d++) {
      x[3*b+d] = B->X[d];
      g[3*b+d] = B->SRC[d];
    }
  }
}

extern "C" void fmm_biot_savart_(int & numBodies, double * x, double * g, double * u) {
  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(numBodies);
  B_iter B=bodies.begin();
  for (int b=0; b<numBodies; b++,B++) {
    B->IBODY = b;
    for (int d=0; d<3; d++) {
      B->X[d] = x[3*b+d];
      B->SRC[d] = g[3*b+d];
    }
    B->SRC[3] = 1e-6;
    B->TRG = 0;
  }
  FMM->stopTimer("Set bodies   ",printNow);

  FMM->startTimer("Build tree   ");
  Cells cells;
  FMM->bottomup(bodies,cells);
  FMM->stopTimer("Build tree   ",printNow);

  FMM->startTimer("Set LET      ");
  FMM->setCommBodies(cells);
  FMM->stopTimer("Set LET      ",printNow);

  FMM->startTimer("Send LET     ");
  Bodies jbodies;
  Cells jcells = cells;
  if( MPISIZE != 1 ) {
#pragma omp parallel sections num_threads(2)
    {
#pragma omp section
      {
        FMM->downward(cells,jcells,1,false);
      }
#pragma omp section
      {
        FMM->updateBodies();
      }
    }
    jbodies = bodies;
    jcells = cells;
    FMM->commCells(jbodies,jcells);
    FMM->eraseLocalTree(jcells);
  }
  FMM->stopTimer("Send LET     ",printNow);

  FMM->startTimer("Traversal    ");
  FMM->downward(cells,jcells,1);
  FMM->stopTimer("Traversal    ",printNow);

  FMM->startTimer("Get bodies   ");
  B=bodies.begin();
  for(int b=0; b<numBodies; b++,B++) {
    for (int d=0; d<3; d++) {
      u[3*B->IBODY+d] = B->TRG[d];
    }
  }
  FMM->stopTimer("Get bodies   ",printNow);
}

extern "C" void direct_biot_savart_(int & numBodies, double * x, double * g, double * u) {
  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(numBodies);
  B_iter B=bodies.begin();
  for (int b=0; b<numBodies; b++,B++) {
    B->IBODY = b;
    B->IPROC = 0;
    for (int d=0; d<3; d++) {
      B->X[d] = x[3*b+d];
      B->SRC[d] = g[3*b+d];
    }
    B->SRC[3] = 1e-6;
    B->TRG = 0;
  }
  FMM->stopTimer("Set bodies   ",printNow);

  FMM->startTimer("Direct sum   ");
  Bodies jbodies = bodies;
  for( int i=0; i!=MPISIZE; ++i ) {
    FMM->shiftBodies(jbodies);
    FMM->evalP2P(bodies,jbodies);
    if(printNow) std::cout << "Direct loop   : " << i+1 << "/" << MPISIZE << std::endl;
  }
  FMM->stopTimer("Direct sum   ",printNow);

  FMM->startTimer("Get bodies   ");
  B=bodies.begin();
  for(int b=0; b<numBodies; b++,B++) {
    for (int d=0; d<3; d++) {
      u[3*B->IBODY+d] = B->TRG[d];
    }
  }
  FMM->stopTimer("Get bodies   ",printNow);
}
