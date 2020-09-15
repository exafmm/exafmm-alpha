#include "types.h"
#include "serialfmm.h"

SerialFMM * FMM;

extern "C" void fmm_init_(int & images, double * x0, double & r0) {
  FMM = new SerialFMM;
  FMM->setKernel("BiotSavart");
  FMM->initialize();
  FMM->IMAGES = images;
  for (int d=0; d<3; d++) FMM->X0[d] = x0[d];
  FMM->R0 = r0;
  FMM->THETA = 1 / sqrtf(4);
}

extern "C" void fmm_finalize_() {
  FMM->finalize();
  delete FMM;
}

extern "C" void fmm_biot_savart_(int & ni, double * xi, double * ui, int & nj, double * xj, double * gj) {
//  bool printNow = true;
  bool printNow = false;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(ni);
  B_iter B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
    B->IBODY = b;
    B->IPROC = 0;
    for (int d=0; d<3; d++) {
      B->X[d] = xi[3*b+d];
    }
    B->TRG = 0;
  }
  Bodies jbodies(nj);
  B=jbodies.begin();
  for (int b=0; b<nj; b++,B++) {
    B->IBODY = b;
    B->IPROC = 0;
    for (int d=0; d<3; d++) {
      B->X[d] = xj[3*b+d];
      B->SRC[d] = gj[3*b+d];
    }
    B->SRC[3] = 1e-6;
  }
  FMM->stopTimer("Set bodies   ",printNow);

  FMM->startTimer("Build tree   ");
  Cells cells, jcells;
  FMM->setDomain(bodies);
  FMM->bottomup(bodies,cells);
  FMM->setDomain(jbodies);
  FMM->bottomup(jbodies,jcells);
  FMM->stopTimer("Build tree   ",printNow);

  FMM->startTimer("Traversal    ");
  FMM->downward(cells,jcells,1);
  FMM->stopTimer("Traversal    ",printNow);

  FMM->startTimer("Get bodies   ");
  B=bodies.begin();
  for(int b=0; b<ni; b++,B++) {
    for (int d=0; d<3; d++) {
      ui[3*B->IBODY+d] = B->TRG[d];
    }
  }
  FMM->stopTimer("Get bodies   ",printNow);
}

extern "C" void direct_biot_savart_(int & ni, double * xi, double * ui, int & nj, double * xj, double * gj) {
//  bool printNow = true;
  bool printNow = false;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(ni);
  B_iter B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
    B->IBODY = b;
    B->IPROC = 0;
    for (int d=0; d<3; d++) {
      B->X[d] = xi[3*b+d];
    }
    B->TRG = 0;
  }
  Bodies jbodies(nj);
  B=jbodies.begin();
  for (int b=0; b<nj; b++,B++) {
    B->IBODY = b;
    B->IPROC = 0;
    for (int d=0; d<3; d++) {
      B->X[d] = xj[3*b+d];
      B->SRC[d] = gj[3*b+d];
    }
    B->SRC[3] = 1e-6;
  }
  FMM->stopTimer("Set bodies   ",printNow);

  FMM->startTimer("Direct sum   ");
  FMM->evalP2P(bodies,jbodies);
  FMM->stopTimer("Direct sum   ",printNow);

  FMM->startTimer("Get bodies   ");
  B=bodies.begin();
  for(int b=0; b<ni; b++,B++) {
    for (int d=0; d<3; d++) {
      ui[3*B->IBODY+d] = B->TRG[d];
    }
  }
  FMM->stopTimer("Get bodies   ",printNow);
}
