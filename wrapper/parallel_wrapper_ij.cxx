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

// partition of original version
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

  assert(bodies.size() < 2*ni);
  assert(jbodies.size() < 2*nj);
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


// 1 package version of partition, biot, unpartion with the same order output as VFM using nioff for giving global index
// --- begin part_fmm_biot_unpart_order ---
extern "C" void part_fmm_biot_unpart_order_(int & ni, int & nioff, double * xi, double * ui, int & nj, double * xj, double * gj) {
  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(ni);
  B_iter B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
//    B->IBODY = b;             // org
    B->IBODY = b+nioff*MPIRANK; // get global index for VFM order
    B->IPROC = 0;
//    std::cout << "b1: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
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

  FMM->startTimer("Set domain in");
  bodies.insert(bodies.end(),jbodies.begin(),jbodies.end());
  FMM->stopTimer("Set domain in",printNow);
  FMM->startTimer("Set domain G ");
  FMM->setGlobDomain(bodies);
  FMM->stopTimer("Set domain G " ,printNow);
  FMM->startTimer("Set domain re");
  bodies.resize(ni);
  FMM->stopTimer("Set domain re",printNow);

  FMM->startTimer("Partition    ");
  FMM->octsection(bodies);
  FMM->octsection(jbodies);
  FMM->stopTimer("Partition    ",printNow);

  FMM->startTimer("Get bodies   ");
  ni = bodies.size();
  nj = jbodies.size();
  assert(ni*nj != 0);
  B=bodies.begin();
  for(int b=0; b<ni; b++,B++) {
//    std::cout << "b2: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
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
  FMM->stopTimer("Get bodies   ",printNow);


//extern "C" void fmm_biot_savart_(int & ni, double * xi, double * ui, int & nj, double * xj, double * gj) {
//  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
//  Bodies bodies(ni);
//  B_iter B=bodies.begin();
  int ibody[ni];           // declare int array of ibody to save global index
  B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
    ibody[b] =  B->IBODY; // added to save B->IBODY of global index
    B->IBODY = b;         // B->IBODY becomes descending order after FMM
    B->IPROC = 0;
    for (int d=0; d<3; d++) {
      B->X[d] = xi[3*b+d];
    }
    B->TRG = 0;
  }
//  Bodies jbodies(nj);
  B=jbodies.begin();
  for (int b=0; b<nj; b++,B++) {
//    B->IBODY = b;
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
      ui[3*B->IBODY+d] = B->TRG[d];  // B->IBODY becomes descending order after FMM
    }
  }
  FMM->stopTimer("Get bodies   ",printNow);
//}


//extern "C" void fmm_unpartition_disorder_xi_ui_ibody_(int & ni, int * ibody, double * xi, double * ui, int & nj, double * xj, double * gj) {
//  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
//  Bodies bodies(ni);
//  B_iter B=bodies.begin();
  B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
//    B->IBODY = b;      // org
    B->IBODY = ibody[b]; // to give global index
//    std::cout << "b3: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
    B->IPROC = 0;
    for (int d=0; d<3; d++) {
      B->X[d] = xi[3*b+d];
      B->TRG[d] = ui[3*b+d]; // added for ui
    }
  }
//!//  Bodies jbodies(nj);
//!  B=jbodies.begin();
//!  for (int b=0; b<nj; b++,B++) {
//!//    B->IBODY = b;
//!    B->IPROC = 0;
//!    for (int d=0; d<3; d++) {
//!      B->X[d] = xj[3*b+d];
//!      B->SRC[d] = gj[3*b+d];
//!    }
//!  }
  FMM->stopTimer("Set bodies   ",printNow);

  FMM->startTimer("Set domain   ");
//!  bodies.insert(bodies.end(),jbodies.begin(),jbodies.end());
  FMM->setGlobDomain(bodies);
//!  bodies.resize(ni);
  FMM->stopTimer("Set domain   ",printNow);

  FMM->startTimer("UnPartition  ");
  FMM->unpartition(bodies);
//!  FMM->unpartition(jbodies);
  FMM->stopTimer("UnPartition  ",printNow);

  FMM->startTimer("Get bodies   ");
  ni = bodies.size();
//!  nj = jbodies.size();
//  assert(ni*nj != 0); // comment out for error in mpirank > 0
  B=bodies.begin();
  for(int b=0; b<ni; b++,B++) {
//    std::cout << "b4: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
    for (int d=0; d<3; d++) {
      xi[3*B->IBODY+d] = B->X[d];   // arrange global index order by B->IBODY
      ui[3*B->IBODY+d] = B->TRG[d]; // added for ui, arrange global index order by B->IBODY
    }
  }
//!  B=jbodies.begin();
//!  for(int b=0; b<nj; b++,B++) {
//!    for (int d=0; d<3; d++) {
//!      xj[3*b+d] = B->X[d];
//!      gj[3*b+d] = B->SRC[d];
//!//      xj[3*B->IBODY+d] = B->X[d];
//!//      gj[3*B->IBODY+d] = B->SRC[d];
//!    }
//!  }
  FMM->stopTimer("Get bodies   ",printNow);

}
// --- end part_fmm_biot_unpart_order ---



// 1 package version of partition, biot, unpartion with the same order output as VFM using nioff for giving global index
// --- begin part_direct_biot_unpart_order ---
extern "C" void part_direct_biot_unpart_order_(int & ni, int & nioff, double * xi, double * ui, int & nj, double * xj, double * gj) {
  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(ni);
  B_iter B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
//    B->IBODY = b;             // org
    B->IBODY = b+nioff*MPIRANK; // get global index for VFM order
    B->IPROC = 0;
//    std::cout << "b1: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
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
//    std::cout << "b2: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
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


//extern "C" void direct_biot_savart_(int & ni, double * xi, double * ui, int & nj, double * xj, double * gj) {
//  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
//  Bodies bodies(ni);
//  B_iter B=bodies.begin();
  int ibody[ni];           // declare int array of ibody to save global index
  B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
    ibody[b] =  B->IBODY; // added to save B->IBODY of global index
    B->IBODY = b;
    B->IPROC = 0;
    for (int d=0; d<3; d++) {
      B->X[d] = xi[3*b+d];
    }
    B->TRG = 0;
  }
//  Bodies jbodies(nj);
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
//}


//extern "C" void fmm_unpartition_disorder_xi_ui_ibody_(int & ni, int * ibody, double * xi, double * ui, int & nj, double * xj, double * gj) {
//  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
//  Bodies bodies(ni);
//  B_iter B=bodies.begin();
  B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
//    B->IBODY = b;      // org
    B->IBODY = ibody[b]; // to give global index
//    std::cout << "b3: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
    B->IPROC = 0;
    for (int d=0; d<3; d++) {
      B->X[d] = xi[3*b+d];
      B->TRG[d] = ui[3*b+d]; // added for ui
    }
  }
//!//  Bodies jbodies(nj);
//!  B=jbodies.begin();
//!  for (int b=0; b<nj; b++,B++) {
//!//    B->IBODY = b;
//!    B->IPROC = 0;
//!    for (int d=0; d<3; d++) {
//!      B->X[d] = xj[3*b+d];
//!      B->SRC[d] = gj[3*b+d];
//!    }
//!  }
  FMM->stopTimer("Set bodies   ",printNow);

  FMM->startTimer("Set domain   ");
//!  bodies.insert(bodies.end(),jbodies.begin(),jbodies.end());
  FMM->setGlobDomain(bodies);
//!  bodies.resize(ni);
  FMM->stopTimer("Set domain   ",printNow);

  FMM->startTimer("UnPartition    ");
  FMM->unpartition(bodies);
//!  FMM->unpartition(jbodies);
  FMM->stopTimer("UnPartition    ",printNow);

  ni = bodies.size();
//!  nj = jbodies.size();
//  assert(ni*nj != 0); // comment out for error in mpirank > 0
  B=bodies.begin();
  for(int b=0; b<ni; b++,B++) {
//    std::cout << "b4: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
    for (int d=0; d<3; d++) {
      xi[3*B->IBODY+d] = B->X[d];   // arrange global index order by B->IBODY
      ui[3*B->IBODY+d] = B->TRG[d]; // added for ui, arrange global index order by B->IBODY
    }
  }
//!  B=jbodies.begin();
//!  for(int b=0; b<nj; b++,B++) {
//!    for (int d=0; d<3; d++) {
//!      xj[3*b+d] = B->X[d];
//!      gj[3*b+d] = B->SRC[d];
//!//      xj[3*B->IBODY+d] = B->X[d];
//!//      gj[3*B->IBODY+d] = B->SRC[d];
//!    }
//!  }


}
// --- end part_direct_biot_unpart_order ---



// partition using ibody with nioff for giving global index
extern "C" void fmm_partition_ibody_(int & ni, int & nioff, int * ibody, double * xi, int & nj, double * xj, double * gj) {
  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(ni);
  B_iter B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
//    B->IBODY = b;
    B->IBODY = b+nioff*MPIRANK;
    B->IPROC = 0;
//    std::cout << "b1: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
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
    ibody[b] = B->IBODY;
//    std::cout << "b2: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
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

// unpartition with disorder, not unsorted
extern "C" void fmm_unpartition_disorder_(int & ni, double * xi, int & nj, double * xj, double * gj) {
  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(ni);
  B_iter B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
    B->IBODY = b;
//    std::cout << "b3: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
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

  FMM->startTimer("UnPartition    ");
  FMM->unpartition(bodies);
  FMM->unpartition(jbodies);
  FMM->stopTimer("UnPartition    ",printNow);

  ni = bodies.size();
  nj = jbodies.size();
//  assert(ni*nj != 0); // comment out for error in mpirank > 0
  B=bodies.begin();
  for(int b=0; b<ni; b++,B++) {
//    std::cout << "b4: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
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

// unpartition for xi, ui with disorder, not unsorted
extern "C" void fmm_unpartition_disorder_xi_ui_(int & ni, double * xi, double * ui, int & nj, double * xj, double * gj) {
  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(ni);
  B_iter B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
    B->IBODY = b;
//    std::cout << "b3: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
    B->IPROC = 0;
    for (int d=0; d<3; d++) {
      B->X[d] = xi[3*b+d];
      B->TRG[d] = ui[3*b+d]; // added for ui
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

  FMM->startTimer("UnPartition    ");
  FMM->unpartition(bodies);
  FMM->unpartition(jbodies);
  FMM->stopTimer("UnPartition    ",printNow);

  ni = bodies.size();
  nj = jbodies.size();
//  assert(ni*nj != 0); // comment out for error in mpirank > 0
  B=bodies.begin();
  for(int b=0; b<ni; b++,B++) {
//    std::cout << "b4: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
    for (int d=0; d<3; d++) {
      xi[3*b+d] = B->X[d];
      ui[3*b+d] = B->TRG[d]; // added for ui
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

// unpartition for xi, ui using ibody with disorder, not unsorted
extern "C" void fmm_unpartition_disorder_xi_ui_ibody_(int & ni, int * ibody, double * xi, double * ui, int & nj, double * xj, double * gj) {
  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(ni);
  B_iter B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
//    B->IBODY = b;
    B->IBODY = ibody[b];
//    std::cout << "b3: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
    B->IPROC = 0;
    for (int d=0; d<3; d++) {
      B->X[d] = xi[3*b+d];
      B->TRG[d] = ui[3*b+d]; // added for ui
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

  FMM->startTimer("UnPartition    ");
  FMM->unpartition(bodies);
  FMM->unpartition(jbodies);
  FMM->stopTimer("UnPartition    ",printNow);

  ni = bodies.size();
  nj = jbodies.size();
//  assert(ni*nj != 0); // comment out for error in mpirank > 0
  B=bodies.begin();
  for(int b=0; b<ni; b++,B++) {
//    std::cout << "b4: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
    for (int d=0; d<3; d++) {
//      xi[3*b+d] = B->X[d];
      xi[3*B->IBODY+d] = B->X[d];
      ui[3*b+d] = B->TRG[d]; // added for ui
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

// unpartition for ud with disorder, not unsorted
extern "C" void fmm_unpartition_disorder_ud_(int & ni, double * ud) {
  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(ni);
  B_iter B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
    B->IBODY = b;
    B->IPROC = 0;
    for (int d=0; d<3; d++) {
//      B->X[d] = xi[3*b+d];
      B->TRG[d] = ud[3*b+d]; // added for ud
    }
  }
  FMM->stopTimer("Set bodies   ",printNow);

  FMM->startTimer("Set domain   ");
  FMM->setGlobDomain(bodies);
  FMM->stopTimer("Set domain   ",printNow);

  FMM->startTimer("UnPartition    ");
  FMM->unpartition(bodies);
  FMM->stopTimer("UnPartition    ",printNow);

  ni = bodies.size();
//  assert(ni*nj != 0); // comment out for error in mpirank > 0
  B=bodies.begin();
  for(int b=0; b<ni; b++,B++) {
    for (int d=0; d<3; d++) {
//      xi[3*b+d] = B->X[d];
      ud[3*b+d] = B->TRG[d]; // added for ud
    }
  }
}

// get global index info after partition using partition and unpartition
extern "C" void fmm_get_imap_(int & ni, int & nioff, double * xi, int * imap) {
  bool printNow = MPIRANK == 0;
  FMM->startTimer("Set bodies   ");
  Bodies bodies(ni);
  B_iter B=bodies.begin();
  for (int b=0; b<ni; b++,B++) {
//    B->IBODY = b;
    B->IBODY = b+nioff*MPIRANK;
    B->IPROC = 0;
    for (int d=0; d<3; d++) {  // begin get global index
      B->X[d] = xi[3*b+d];     // get global index
    }                          // end get global index
  }
  FMM->stopTimer("Set bodies   ",printNow);

  FMM->startTimer("Set domain   ");
  FMM->setGlobDomain(bodies);
  FMM->stopTimer("Set domain   ",printNow);

  FMM->startTimer("Partition    ");
  FMM->octsection(bodies);
  FMM->stopTimer("Partition    ",printNow);

  ni = bodies.size();
  B=bodies.begin();
//  for(int b=0; b<ni; b++,B++) {
//    std::cout << "b: " << B->IBODY << "\n"; // write i_global : b=map[i_global]
//  }

  FMM->startTimer("Set bodies   ");
  for (int b=0; b<ni; b++,B++) {
//    B->IBODY = b;
    B->IPROC = 0;
  }
  FMM->stopTimer("Set bodies   ",printNow);

  FMM->startTimer("Set domain   ");
  FMM->setGlobDomain(bodies);
  FMM->stopTimer("Set domain   ",printNow);

  FMM->startTimer("UnPartition    ");
  FMM->unpartition(bodies);
  FMM->stopTimer("UnPartition    ",printNow);

  std::map<int, int> map;       // set map array

  ni = bodies.size();
  B=bodies.begin();
  for(int b=0; b<ni; b++,B++) {
    for (int d=0; d<3; d++) {  // begin give global index after unpartition
      xi[3*b+d] = B->X[d];     // give global index after unpartition
    }                          // end give global index after unpartition

//    std::cout << "bu: " << B->IBODY << "\n"; // write i_global : b=map[i_global]

    int i_global = B->IBODY;   // get global index for map array
    map[i_global] = b+1;       // get map array
  }
//  for (auto itr = map.begin(); itr != map.end(); ++itr) {     // begin write i_global : b=map[i_global]
//      std::cout << itr->first << ": " << itr->second << "\n"; // end write i_global : b=map[i_global]
//  }
  for(int i_global=0; i_global<ni; i_global++) { // begin get imap array for fortran
    imap[i_global] = map[i_global];            // get imap array for fortran
  }                                              // end get imap array for fortran
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
