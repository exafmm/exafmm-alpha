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

extern "C" void fmm_partition_(int & ni, double * xi, int & nj, double * xj, double * gj) {
  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(ni);
  B_iter B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
    B->IBODY = b;
    B->IPROC = 0;
    for (int d=0; d<3; d++) {
      B->X[d] = xi[3*b+d];
    }
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
  }
  FMM->stopTimer("Set bodies   ",printNow);

  FMM->startTimer("Set domain   ");
  bodies.insert(bodies.end(),jbodies.begin(),jbodies.end());
  FMM->setGlobDomain(bodies);
  bodies.resize(ni);
  FMM->stopTimer("Set domain   ",printNow);

  FMM->startTimer("Partition    ");
  FMM->octsection(bodies);
  FMM->octsection(jbodies);
  FMM->stopTimer("Partition    ",printNow);

  ni = bodies.size();
  nj = jbodies.size();
  assert(ni*nj != 0);
  B=bodies.begin();
  for(int b=0; b<ni; b++,B++) {
    for (int d=0; d<3; d++) {
      xi[3*b+d] = B->X[d];
    }
  }
  B=jbodies.begin();
  for(int b=0; b<nj; b++,B++) {
    for (int d=0; d<3; d++) {
      xj[3*b+d] = B->X[d];
      gj[3*b+d] = B->SRC[d];
    }
  }
}

extern "C" void fmm_biot_savart_(int & ni, double * xi, double * ui, int & nj, double * xj, double * gj) {
  bool printNow = MPIRANK == 0;
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
  FMM->bottomup(bodies,cells);
  FMM->bottomup(jbodies,jcells);
  FMM->stopTimer("Build tree   ",printNow);

  FMM->startTimer("Set LET      ");
  Bodies jbodies2 = jbodies; // Remove this
  Cells jcells2 = jcells; // Remove this
  FMM->setCommBodies(jcells);
  FMM->stopTimer("Set LET      ",printNow);

  FMM->startTimer("Send LET     ");
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
    jbodies = jbodies2;
    jcells = jcells2;
    FMM->commCells(jbodies,jcells);
    FMM->eraseLocalTree(jcells);
  }
  FMM->stopTimer("Send LET     ",printNow);

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
  bool printNow = MPIRANK == 0;
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
  for( int i=0; i!=MPISIZE; ++i ) {
    FMM->shiftBodies(jbodies);
    FMM->evalP2P(bodies,jbodies);
    if(printNow) std::cout << "Direct loop   : " << i+1 << "/" << MPISIZE << std::endl;
  }
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
