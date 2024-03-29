void Evaluator::testMACP2P(C_iter Ci, C_iter Cj) {              // Test multipole acceptance criteria for P2P kernel
  BI0 = Ci->LEAF;                                               // Set target bodies begin iterator
  BIN = Ci->LEAF + Ci->NLEAF;                                   // Set target bodies end iterator
  BJ0 = Cj->LEAF;                                               // Set source bodies begin iterator
  BJN = Cj->LEAF + Cj->NLEAF;                                   // Set source bodies end iterator
  vect Xp = Xperiodic;                                          // Define coordinate offset vector
  selectP2P_CPU(Xp);                                            // Select P2P_CPU kernel
  NP2P++;                                                       // Count P2P kernel execution
}

void Evaluator::testMACM2L(C_iter Ci, C_iter Cj) {              // Test multipole acceptance criteria for M2L kernel
  vect dist = Ci->X - Cj->X - Xperiodic;                        // Distance vector between cells
  real R = std::sqrt(norm(dist));                               // Distance between cells
  if( Ci->R + Cj->R > THETA*R ) {                               // If cell is too large
    Pair pair(Ci,Cj);                                           //  Form pair of interacting cells
    pairs.push(pair);                                           //  Push interacting pair into stack
  } else {                                                      // If cell is small enough
    CI = Ci;                                                    //  Set global target iterator
    CJ = Cj;                                                    //  Set global source iterator
    selectM2L();                                                //  Select M2L kernel
    NM2L++;                                                     //  Count M2L kernel execution
  }                                                             // Endif for interaction
}

void Evaluator::testMACM2P(C_iter Ci, C_iter Cj) {              // Test multipole acceptance criteria for M2P kernel
  vect dist = Ci->X - Cj->X - Xperiodic;                        // Distance vector between cells
  real R = std::sqrt(norm(dist));                               // Distance between cells
  if( Ci->NCHILD != 0 || Ci->R + Cj->R > THETA*R ) {            // If target is not twig or cell is too large
    Pair pair(Ci,Cj);                                           //  Form pair of interacting cells
    pairs.push(pair);                                           //  Push interacting pair into stack
  } else {                                                      // If target is twig and cell is small enough
    CI = Ci;                                                    //  Set global target iterator
    CJ = Cj;                                                    //  Set global source iterator
    selectM2P();                                                //  Select M2P kernel
    NM2P++;                                                     //  Count M2P kernel execution
  }                                                             // Endif for interaction
}

void Evaluator::timeKernels() {                                 // Time all kernels for auto-tuning
  Bodies ibodies(1000), jbodies(1000);                          // Artificial bodies
  for( B_iter Bi=ibodies.begin(),Bj=jbodies.begin(); Bi!=ibodies.end(); ++Bi, ++Bj ) {// Loop over artificial bodies
    Bi->X = 0;                                                  //  Set coordinates of target body
    Bj->X = 1;                                                  //  Set coordinates of source body
  }                                                             // End loop over artificial bodies
  Cells cells;                                                  // Artificial cells
  cells.resize(2);                                              // Two artificial cells
  C_iter Ci = cells.begin(), Cj = cells.begin()+1;              // Artificial target & source cell
  Ci->X = 0;                                                    // Set coordinates of target cell
  Ci->NLEAF = 10;                                               // Number of leafs in target cell
  Ci->LEAF = ibodies.begin();                                   // Leaf iterator in target cell
  Cj->X = 1;                                                    // Set coordinates of source cell
  Cj->NLEAF = 1000;                                             // Number of leafs in source cell
  Cj->LEAF = jbodies.begin();                                   // Leaf iterator in source cell
  BI0 = Ci->LEAF;                                               // Set target bodies begin iterator
  BIN = Ci->LEAF + Ci->NLEAF;                                   // Set target bodies end iterator
  BJ0 = Cj->LEAF;                                               // Set source bodies begin iterator
  BJN = Cj->LEAF + Cj->NLEAF;                                   // Set source bodies end iterator
  vect Xp = 0;                                                  // Define coordinate offset vector
  startTimer("P2P kernel   ");                                  // Start timer
  for( int i=0; i!=1; ++i ) selectP2P_CPU(Xp);                  // Select P2P_CPU kernel
  timeP2P = stopTimer("P2P kernel   ") / 10000;                 // Stop timer
  CI = Ci;                                                      // Set global target cell iterator
  CJ = Cj;                                                      // Set global source cell iterator
  startTimer("M2L kernel   ");                                  // Start timer
  for( int i=0; i!=1000; ++i ) selectM2L();                     // Select M2L kernel
  timeM2L = stopTimer("M2L kernel   ") / 1000;                  // Stop timer
  startTimer("M2P kernel   ");                                  // Start timer
  for( int i=0; i!=100; ++i ) selectM2P();                      // Select M2P kernel
  timeM2P = stopTimer("M2P kernel   ") / 1000;                  // Stop timer
}

void Evaluator::traversePeriodic(Cells &cells, Cells &jcells, int method) {// Traverse tree for periodic cells
  Xperiodic = 0;                                                // Set periodic coordinate offset
  C_iter Cj = jcells.end()-1;                                   // Initialize iterator for periodic source cell
  for( int level=0; level<IMAGES-1; ++level ) {                 // Loop over sublevels of tree
    for( int I=0; I!=26; ++I, --Cj ) {                          //  Loop over periodic images (exclude center)
      switch (method) {                                         //   Switch between method
      case 0 :                                                  //   0 : treecode
        for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) { //   Loop over cells
          if( Ci->NCHILD == 0 ) {                               //     If cell is twig
            CI = Ci;                                            //      Set global target iterator
            CJ = Cj;                                            //      Set global source iterator
            selectM2P();                                        //      Select M2P kernel
          }                                                     //     Endif for twig
        }                                                       //    End loop over cells
        break;                                                  //    Terminate this case
      case 1 :                                                  //   1 : FMM
      case 2 :                                                  //   2 : hybrid
        CI = cells.end() - 1;                                   //    Set root cell as target iterator
        CJ = Cj;                                                //    Set global source iterator
        selectM2L();                                            //    Select M2P kernel
        break;                                                  //    Terminate this case
      }                                                         //   End switch between methods
    }                                                           //  End loop over x periodic direction
  }                                                             // End loop over sublevels of tree
}

void Evaluator::evalP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU) {// Evaluate P2P
  BI0 = ibodies.begin();                                        // Set target bodies begin iterator
  BIN = ibodies.end();                                          // Set target bodies end iterator
  BJ0 = jbodies.begin();                                        // Set source bodies begin iterator
  BJN = jbodies.end();                                          // Set source bodies end iterator
  int range = (pow(3,IMAGES) - 1) / 2;                          // Compute periodic range
#pragma omp parallel for collapse(3)
  for( int ix=-range; ix<=range; ++ix ) {                       // Loop over x periodic direction
    for( int iy=-range; iy<=range; ++iy ) {                     //  Loop over y periodic direction
      for( int iz=-range; iz<=range; ++iz ) {                   //   Loop over z periodic direction
	vect Xp;                                                //    Define coordinate offset vector
        Xp[0] = ix * 2 * R0;                                    //    Coordinate offset for x periodic direction
        Xp[1] = iy * 2 * R0;                                    //    Coordinate offset for y periodic direction
        Xp[2] = iz * 2 * R0;                                    //    Coordinate offset for z periodic direction
        selectP2P_CPU(Xp);                                      //    Select P2P_CPU kernel
      }                                                         //   End loop over x periodic direction
    }                                                           //  End loop over y periodic direction
  }                                                             // End loop over z periodic direction
}

void Evaluator::evalP2M(Cells &cells) {                         // Evaluate P2M
  startTimer("evalP2M      ");                                  // Start timer
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over cells
    CI->M = CI->L = 0;                                          //  Initialize multipole & local coefficients
    if( CI->NCHILD == 0 ) {                                     //  If cell is a twig
      selectP2M();                                              //   Select P2M kernel
    }                                                           //  Endif for twig
  }                                                             // End loop over cells
  stopTimer("evalP2M      ");                                   // Stop timer
}

void Evaluator::evalM2M(Cells &cells) {                         // Evaluate M2M
  startTimer("evalM2M      ");                                  // Start timer
  CJ0 = cells.begin();                                          // Set begin iterator
  for( CJ=cells.begin(); CJ!=cells.end()-1; ++CJ ) {            // Loop over cells bottomup (except root cell)
    CI = CJ0 + CJ->PARENT;                                      //  Set target cell iterator
    selectM2M_CPU();                                            //  Select M2M_CPU kernel
  }                                                             // End loop over cells
  stopTimer("evalM2M      ");                                   // Stop timer
}

void Evaluator::evalM2L(Cells &cells, bool kernel) {            // Evaluate M2L
  if( kernel ) {                                                // If this is a kernel unit test
    startTimer("evalM2L      ");                                // Start timer
    CI0 = cells.begin();                                        // Set begin iterator
    for( CI=cells.begin(); CI!=cells.end(); ++CI ) {            // Loop over cells
      while( !listM2L[CI-CI0].empty() ) {                       //  While M2L interaction list is not empty
        CJ = listM2L[CI-CI0].back();                            //   Set source cell iterator
        Iperiodic = flagM2L[CI-CI0][CJ];                        //   Set periodic image flag
        int I = 0;                                              //   Initialize index of periodic image
        for( int ix=-1; ix<=1; ++ix ) {                         //   Loop over x periodic direction
          for( int iy=-1; iy<=1; ++iy ) {                       //    Loop over y periodic direction
            for( int iz=-1; iz<=1; ++iz, ++I ) {                //     Loop over z periodic direction
              if( Iperiodic & (1 << I) ) {                      //      If periodic flag is on
                Xperiodic[0] = ix * 2 * R0;                     //       Coordinate offset for x periodic direction
                Xperiodic[1] = iy * 2 * R0;                     //       Coordinate offset for y periodic direction
                Xperiodic[2] = iz * 2 * R0;                     //       Coordinate offset for z periodic direction
                selectM2L();                                    //       Select M2L kernel
              }                                                 //      Endif for periodic flag
            }                                                   //     End loop over x periodic direction
          }                                                     //    End loop over y periodic direction
        }                                                       //   End loop over z periodic direction
        listM2L[CI-CI0].pop_back();                             //   Pop last element from M2L interaction list
      }                                                         //  End while for M2L interaction list
    }                                                           // End loop over cells topdown
    listM2L.clear();                                            // Clear interaction lists
    flagM2L.clear();                                            // Clear periodic image flags
    stopTimer("evalM2L      ");                                 // Stop timer
  }                                                             // End if for kernel unit test
}

void Evaluator::evalM2P(Cells &cells, bool kernel) {            // Evaluate M2P
  if( kernel ) {                                                // If this is a kernel unit test
    startTimer("evalM2P      ");                                // Start timer
    CI0 = cells.begin();                                        // Set begin iterator
    for( CI=cells.begin(); CI!=cells.end(); ++CI ) {            // Loop over cells
      while( !listM2P[CI-CI0].empty() ) {                       //  While M2P interaction list is not empty
        CJ = listM2P[CI-CI0].back();                            //   Set source cell iterator
        Iperiodic = flagM2P[CI-CI0][CJ];                        //   Set periodic image flag
        int I = 0;                                              //   Initialize index of periodic image
        for( int ix=-1; ix<=1; ++ix ) {                         //   Loop over x periodic direction
          for( int iy=-1; iy<=1; ++iy ) {                       //    Loop over y periodic direction
            for( int iz=-1; iz<=1; ++iz, ++I ) {                //     Loop over z periodic direction
              if( Iperiodic & (1 << I) ) {                      //      If periodic flag is on
                Xperiodic[0] = ix * 2 * R0;                     //       Coordinate offset for x periodic direction
                Xperiodic[1] = iy * 2 * R0;                     //       Coordinate offset for y periodic direction
                Xperiodic[2] = iz * 2 * R0;                     //       Coordinate offset for z periodic direction
                selectM2P();                                    //       Select M2P kernel
              }                                                 //      Endif for periodic flag
            }                                                   //     End loop over x periodic direction
          }                                                     //    End loop over y periodic direction
        }                                                       //   End loop over z periodic direction
        listM2P[CI-CI0].pop_back();                             //   Pop last element from M2P interaction list
      }                                                         //  End while for M2P interaction list
    }                                                           // End loop over cells topdown
    listM2P.clear();                                            // Clear interaction lists
    flagM2P.clear();                                            // Clear periodic image flags
    stopTimer("evalM2P      ");                                 // Stop timer
  }                                                             // End if for kernel unit test
}

void Evaluator::evalP2P(Cells &cells, bool kernel) {            // Evaluate P2P
  if( kernel ) {                                                // If this is a kernel unit test
    startTimer("evalP2P      ");                                // Start timer
    CI0 = cells.begin();                                        // Set begin iterator
    for( CI=cells.begin(); CI!=cells.end(); ++CI ) {            // Loop over cells
      BI0 = CI->LEAF;                                           //  Set target bodies begin iterator
      BIN = CI->LEAF + CI->NLEAF;                               //  Set target bodies end iterator
      while( !listP2P[CI-CI0].empty() ) {                       //  While M2P interaction list is not empty
        CJ = listP2P[CI-CI0].back();                            //   Set source cell iterator
        BJ0 = CJ->LEAF;                                         //   Set source bodies begin iterator
        BJN = CJ->LEAF + CJ->NLEAF;                             //   Set source bodies end iterator
        Iperiodic = flagP2P[CI-CI0][CJ];                        //   Set periodic image flag
        int I = 0;                                              //   Initialize index of periodic image
        for( int ix=-1; ix<=1; ++ix ) {                         //   Loop over x periodic direction
          for( int iy=-1; iy<=1; ++iy ) {                       //    Loop over y periodic direction
            for( int iz=-1; iz<=1; ++iz, ++I ) {                //     Loop over z periodic direction
              if( Iperiodic & (1 << I) ) {                      //      If periodic flag is on
		vect Xp;                                        //       Define coordinate offset vector
                Xp[0] = ix * 2 * R0;                            //       Coordinate offset for x periodic direction
                Xp[1] = iy * 2 * R0;                            //       Coordinate offset for y periodic direction
                Xp[2] = iz * 2 * R0;                            //       Coordinate offset for z periodic direction
                selectP2P_CPU(Xp);                              //       Select P2P_CPU kernel
              }                                                 //      Endif for periodic flag
            }                                                   //     End loop over x periodic direction
          }                                                     //    End loop over y periodic direction
        }                                                       //   End loop over z periodic direction
        listP2P[CI-CI0].pop_back();                             //   Pop last element from M2P interaction list
      }                                                         //  End while for M2P interaction list
    }                                                           // End loop over cells topdown
    listP2P.clear();                                            // Clear interaction lists
    flagP2P.clear();                                            // Clear periodic image flags
    stopTimer("evalP2P      ");                                 // Stop timer
  }                                                             // End if for kernel unit test
}

void Evaluator::evalL2L(Cells &cells) {                         // Evaluate L2L
  startTimer("evalL2L      ");                                  // Start timer
  CI0 = cells.begin();                                          // Set begin iterator
  for( CI=cells.end()-2; CI!=cells.begin()-1; --CI ) {          // Loop over cells topdown (except root cell)
    CJ = CI0 + CI->PARENT;                                      //  Set source cell iterator
    selectL2L();                                                //  Select L2L kernel
  }                                                             // End loop over cells topdown
  stopTimer("evalL2L      ");                                   // Stop timer
}

void Evaluator::evalL2P(Cells &cells) {                         // Evaluate L2P
  startTimer("evalL2P      ");                                  // Start timer
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over cells
    if( CI->NCHILD == 0 ) {                                     //  If cell is a twig
      selectL2P();                                              //   Select L2P kernel
    }                                                           //  Endif for twig
  }                                                             // End loop over cells topdown
  stopTimer("evalL2P      ");                                   // Stop timer
}
