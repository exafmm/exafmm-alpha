#include "types.h"
#include "serialfmm.h"

SerialFMM * FMM;

extern "C" void fmm_init_() {
  FMM = new SerialFMM;
  FMM->setKernel("BiotSavart");
  FMM->initialize();
  FMM->IMAGES = 0;
  FMM->THETA = 1 / sqrtf(4);
}

extern "C" void fmm_finalize_() {
  FMM->finalize();
  delete FMM;
}

extern "C" void fmm_biot_savart_(int & numBodies, double * x, double * g, double * u) {
  bool printNow = true;
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

  FMM->startTimer("Set domain   ");
  FMM->setDomain(bodies);
  FMM->stopTimer("Set domain   ",printNow);

  FMM->startTimer("Build tree   ");
  Cells cells;
  FMM->bottomup(bodies,cells);
  FMM->stopTimer("Build tree   ",printNow);

  FMM->startTimer("Traversal    ");
  Cells jcells = cells;
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
  bool printNow = true;
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
  FMM->evalP2P(bodies,jbodies);
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
