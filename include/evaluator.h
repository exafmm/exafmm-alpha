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

public:
//! Constructor
  Evaluator() {}
//! Destructor
  ~Evaluator() {}

//! Add single list for kernel unit test
  void addM2L(C_iter Cj) {
    listM2L.resize(1);                                          // Resize vector of M2L interation lists
    listM2L[0].push_back(Cj);                                   // Push single cell into list
  }

//! Add single list for kernel unit test
  void addM2P(C_iter Cj) {
    listM2P.resize(1);                                          // Resize vector of M2P interation lists
    listM2P[0].push_back(Cj);                                   // Push single cell into list
  }

//! Traverse tree to get interaction list
  void traverse(Cells &cells, Cells &jcells, int method) {
    C_iter root = cells.end() - 1;                              // Iterator for root target cell
    C_iter jroot = jcells.end() - 1;                            // Iterator for root source cell
    CI0 = cells.begin();                                        // Set begin iterator for target cells
    CJ0 = jcells.begin();                                       // Set begin iterator for source cells
    listM2L.resize(cells.size());                               // Resize M2L interaction list
    listM2P.resize(cells.size());                               // Resize M2P interaction list
    listP2P.resize(cells.size());                               // Resize P2P interaction list
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
  }

//! Initialize GPU
  void initialize() {
    LaplaceInit();
  }

//! Finalize GPU
  void finalize() {
    LaplaceFinal();
  }

  void setSourceBody() {                               // Set source buffer for bodies
    startTimer("Set sourceB  ");                                  // Start timer
    for( M_iter M=sourceSize.begin(); M!=sourceSize.end(); ++M ) {// Loop over source map
      CJ = M->first;                                              //  Set source cell
      sourceBegin[CJ] = sourceHost.size() / 7;                    //  Key : iterator, Value : offset of source leafs
      for( B_iter B=CJ->LEAF; B!=CJ->LEAF+CJ->NLEAF; ++B ) {      //  Loop over leafs in source cell
        sourceHost.push_back(B->X[0]);                            //   Copy x position to GPU buffer
        sourceHost.push_back(B->X[1]);                            //   Copy y position to GPU buffer
        sourceHost.push_back(B->X[2]);                            //   Copy z position to GPU buffer
        sourceHost.push_back(B->SRC[0]);                          //   Copy 1st source value to GPU buffer
        sourceHost.push_back(B->SRC[1]);                          //   Copy 2nd source value to GPU buffer
        sourceHost.push_back(B->SRC[2]);                          //   Copy 3rd source value to GPU buffer
        sourceHost.push_back(B->SRC[3]);                          //   Copy 4th source value to GPU buffer
      }                                                           //  End loop over leafs
    }                                                             // End loop over source map
    stopTimer("Set sourceB  ");                                   // Stop timer
  }

  void setSourceCell(bool isM=true) {                  // Set source buffer for cells
    startTimer("Set sourceC  ");                                  // Start timer
    for( M_iter M=sourceSize.begin(); M!=sourceSize.end(); ++M ) {// Loop over source map
      CJ = M->first;                                              //  Set source cell
      sourceBegin[CJ] = sourceHost.size();                        //  Key : iterator, Value : offset of sources
      sourceHost.push_back(CJ->X[0]);                             //  Copy x position to GPU buffer
      sourceHost.push_back(CJ->X[1]);                             //  Copy y position to GPU buffer
      sourceHost.push_back(CJ->X[2]);                             //  Copy z position to GPU buffer
      if( isM ) {                                                 //  If source is M
        for( int i=0; i!=NCOEF; ++i ) {                           //   Loop over coefs in source cell
          sourceHost.push_back((CJ->M[i]).real());                //    Copy real multipole to GPU buffer
          sourceHost.push_back((CJ->M[i]).imag());                //    Copy imaginary multipole to GPU buffer
        }                                                         //   End loop over coefs
      } else {                                                    //  If source is L
        for( int i=0; i!=NCOEF; ++i ) {                           //   Loop over coefs in source cell
          sourceHost.push_back((CJ->L[i]).real());                //    Copy real multipole to GPU buffer
          sourceHost.push_back((CJ->L[i]).imag());                //    Copy imaginary multipole to GPU buffer
        }                                                         //   End loop over coefs
      }                                                           //  Endif for source type
    }                                                             // End loop over source map
    stopTimer("Set sourceC  ");                                   // Stop timer
  }

  void setTargetBody(Lists lists) {        // Set target buffer for bodies
    startTimer("Set targetB  ");                                  // Start timer
    int key = 0;                                                  // Initialize key to range of coefs in source cells
    for( CI=CIB; CI!=CIE; ++CI ) {                                // Loop over target cells
      if( !lists[CI-CI0].empty() ) {                              //  If the interation list is not empty
        BI0 = CI->LEAF;                                           //   Set target bodies begin iterator
        BIN = CI->LEAF + CI->NLEAF;                               //   Set target bodies end iterator
        int blocks = (BIN - BI0 - 1) / THREADS + 1;               //   Number of thread blocks needed for this target cell
        for( int i=0; i!=blocks; ++i ) {                          //   Loop over thread blocks
          keysHost.push_back(key);                                //    Save key to range of leafs in source cells
        }                                                         //   End loop over thread blocks
        key += 3*lists[CI-CI0].size()+1;                          //   Increment key counter
        rangeHost.push_back(lists[CI-CI0].size());                //   Save size of interaction list
        for( L_iter L=lists[CI-CI0].begin(); L!=lists[CI-CI0].end(); ++L ) {//  Loop over interaction list
          CJ = *L;                                                //   Set source cell
          rangeHost.push_back(sourceBegin[CJ]);                   //    Set begin index of coefs in source cell
          rangeHost.push_back(sourceSize[CJ]);                    //    Set number of coefs in source cell
          rangeHost.push_back(0);                 //    Set periodic image flag of source cell
        }                                                         //   End loop over interaction list
        targetBegin[CI] = targetHost.size() / 6;                  //   Key : iterator, Value : offset of target leafs
        for( B_iter B=BI0; B!=BIN; ++B ) {                        //   Loop over leafs in target cell
          targetHost.push_back(B->X[0]);                          //    Copy x position to GPU buffer
          targetHost.push_back(B->X[1]);                          //    Copy y position to GPU buffer
          targetHost.push_back(B->X[2]);                          //    Copy z position to GPU buffer
          targetHost.push_back(B->SRC[0]);                        //    Copy 1st target value to GPU buffer
          targetHost.push_back(B->SRC[1]);                        //    Copy 2nd target value to GPU buffer
          targetHost.push_back(B->SRC[2]);                        //    Copy 3rd target value to GPU buffer
        }                                                         //   End loop over leafs
        int numPad = blocks * THREADS - (BIN - BI0);              //   Number of elements to pad in target GPU buffer
        for( int i=0; i!=numPad; ++i ) {                          //   Loop over elements to pad
          targetHost.push_back(0);                                //    Pad x position in GPU buffer
          targetHost.push_back(0);                                //    Pad y position in GPU buffer
          targetHost.push_back(0);                                //    Pad z position in GPU buffer
          targetHost.push_back(0);                                //    Pad 1st target value to GPU buffer
          targetHost.push_back(0);                                //    Pad 2nd target value to GPU buffer
          targetHost.push_back(0);                                //    Pad 3rd target value to GPU buffer
        }                                                         //   End loop over elements to pad
      }                                                           //  End if for empty interation list
    }                                                             // End loop over target cells
    stopTimer("Set targetB  ");                                   // Stop timer
  }

  void setTargetCell(Lists lists) {        // Set target buffer for cells
    startTimer("Set targetC  ");                                  // Start timer
    int key = 0;                                                  // Initialize key to range of coefs in target cells
    for( CI=CIB; CI!=CIE; ++CI ) {                                // Loop over target cells
      if( !lists[CI-CI0].empty() ) {                              //  If the interation list is not empty
        keysHost.push_back(key);                                  //   Save key to range of coefs in target cells
        key += 3*lists[CI-CI0].size()+1;                          //   Increment key counter
        rangeHost.push_back(lists[CI-CI0].size());                //   Save size of interaction list
        for( L_iter L=lists[CI-CI0].begin(); L!=lists[CI-CI0].end(); ++L ) {//  Loop over interaction list
          CJ = *L;                                                //    Set target cell
          int begin = sourceBegin[CJ];                            //    Get begin index of coefs in source cell
          int size = sourceSize[CJ];                              //    Get number of coefs in source cell
          rangeHost.push_back(begin);                             //    Set begin index of coefs in source cell
          rangeHost.push_back(size);                              //    Set number of coefs in source cell
          rangeHost.push_back(0);                 //    Set periodic image flag of source cell
        }                                                         //   End loop over interaction list
        targetBegin[CI] = targetHost.size();                      //   Key : iterator, Value : offset of target coefs
        targetHost.push_back(CI->X[0]);                           //   Copy x position to GPU buffer
        targetHost.push_back(CI->X[1]);                           //   Copy y position to GPU buffer
        targetHost.push_back(CI->X[2]);                           //   Copy z position to GPU buffer
        for( int i=0; i!=2*NCOEF; ++i ) {                         //   Loop over coefs in target cell
          targetHost.push_back(0);                                //    Pad GPU buffer
        }                                                         //   End loop over coefs
        int numPad = 2 * THREADS * NCOEF / NTERM - 2 * NCOEF - 3; //   Number of elements to pad in target GPU buffer
        assert(numPad >= 0);                                      //   THREADS must be large enough
        for( int i=0; i!=numPad; ++i ) {                          //   Loop over elements to pad
          targetHost.push_back(0);                                //    Pad GPU buffer
        }                                                         //   End loop over elements to pad
      }                                                           //  End if for empty interation list
    }                                                             // End loop over target cells
    stopTimer("Set targetC  ");                                   // Stop timer
  }

  void getTargetBody(Lists &lists) {                   // Get body values from target buffer
    startTimer("Get targetB  ");                                  // Start timer
    for( CI=CIB; CI!=CIE; ++CI ) {                                // Loop over target cells
      if( !lists[CI-CI0].empty() ) {                              //  If the interation list is not empty
        BI0 = CI->LEAF;                                           //   Set target bodies begin iterator
        BIN = CI->LEAF + CI->NLEAF;                               //   Set target bodies end iterator
        int begin = targetBegin[CI];                              //   Offset of target leafs
        for( B_iter B=BI0; B!=BIN; ++B ) {                        //    Loop over target bodies
          B->TRG[0] += targetHost[6*(begin+B-BI0)+0];             //     Copy 1st target value from GPU buffer
          B->TRG[1] += targetHost[6*(begin+B-BI0)+1];             //     Copy 2nd target value from GPU buffer
          B->TRG[2] += targetHost[6*(begin+B-BI0)+2];             //     Copy 3rd target value from GPU buffer
          B->TRG[3] += targetHost[6*(begin+B-BI0)+3];             //     Copy 4th target value from GPU buffer
        }                                                         //    End loop over target bodies
        lists[CI-CI0].clear();                                    //   Clear interaction list
      }                                                           //  End if for empty interation list
    }                                                             // End loop over target cells
    stopTimer("Get targetB  ");                                   // Stop timer
  }

  void getTargetCell(Lists &lists, bool isM=true) {    // Get body values from target buffer
    startTimer("Get targetC  ");                                  // Start timer
    for( CI=CIB; CI!=CIE; ++CI ) {                                // Loop over target cells
      if( !lists[CI-CI0].empty() ) {                              //  If the interation list is not empty
        int begin = targetBegin[CI];                              //   Offset of target coefs
        if( isM ) {                                               //   If target is M
          for( int i=0; i!=NCOEF; ++i ) {                         //    Loop over coefs in target cell
            CI->M[i] += complex(targetHost[begin+2*i+0],targetHost[begin+2*i+1]); // Copy real target values from GPU buffer
          }                                                       //    End loop over coefs
        } else {                                                  //   If target is L
          for( int i=0; i!=NCOEF; ++i ) {                         //    Loop over coefs in target cell
            CI->L[i] += complex(targetHost[begin+2*i+0],targetHost[begin+2*i+1]); // Copy real target values from GPU buffer
          }                                                       //    End loop over coefs
        }                                                         //   Endif for target type
        lists[CI-CI0].clear();                                    //   Clear interaction list
      }                                                           //  End if for empty interation list
    }                                                             // End loop over target cells
    stopTimer("Get targetC  ");                                   // Stop timer
  }

  void clearBuffers() {                                // Clear GPU buffers
    startTimer("Clear buffer ");                                  // Start timer
    keysHost.clear();                                             // Clear keys vector
    rangeHost.clear();                                            // Clear range vector
    constHost.clear();                                            // Clear const vector
    sourceHost.clear();                                           // Clear source vector
    targetHost.clear();                                           // Clear target vector
    sourceBegin.clear();                                          // Clear map for offset of source cells
    sourceSize.clear();                                           // Clear map for size of source cells
    targetBegin.clear();                                          // Clear map for offset of target cells
    stopTimer("Clear buffer ");                                   // Stop timer
  }

  void testMACP2P(C_iter Ci, C_iter Cj) {              // Test multipole acceptance criteria for P2P kernel
    listP2P[Ci-CI0].push_back(Cj);                                // Push source cell into P2P interaction list
    NP2P++;                                                       // Count P2P kernel execution
  }

  void testMACM2L(C_iter Ci, C_iter Cj) {              // Test multipole acceptance criteria for M2L kernel
    vect dist = Ci->X - Cj->X;                        // Distance vector between cells
    real R = std::sqrt(norm(dist));                               // Distance between cells
    if( Ci->R + Cj->R > THETA*R ) {                               // If cell is too large
      Pair pair(Ci,Cj);                                           //  Form pair of interacting cells
      pairs.push(pair);                                           //  Push interacting pair into stack
    } else {                                                      // If cell is small enough
      listM2L[Ci-CI0].push_back(Cj);                              //  Push source cell into M2L interaction list
      NM2L++;                                                     //  Count M2L kernel execution
    }                                                             // Endif for interaction
  }

  void testMACM2P(C_iter Ci, C_iter Cj) {              // Test multipole acceptance criteria for M2P kernel
    vect dist = Ci->X - Cj->X;                        // Distance vector between cells
    real R = std::sqrt(norm(dist));                               // Distance between cells
    if( Ci->NCHILD != 0 || Ci->R + Cj->R > THETA*R ) {            // If target is not twig or cell is too large
      Pair pair(Ci,Cj);                                           //  Form pair of interacting cells
      pairs.push(pair);                                           //  Push interacting pair into stack
    } else {                                                      // If target is twig and cell is small enough
      listM2P[Ci-CI0].push_back(Cj);                              //  Push source cell into M2P interaction list
      NM2P++;                                                     //  Count M2P kernel execution
    }                                                             // Endif for interaction
  }

  void timeKernels() {                                 // Time all kernels for auto-tuning
    Bodies ibodies(NCRIT), jbodies(NCRIT);                        // Artificial bodies
    for( B_iter Bi=ibodies.begin(),Bj=jbodies.begin(); Bi!=ibodies.end(); ++Bi, ++Bj ) {// Loop over artificial bodies
      Bi->X = 0;                                                  //  Set coordinates of target body
      Bj->X = 1;                                                  //  Set coordinates of source body
    }                                                             // End loop over artificial bodies
    Cells icells, jcells;                                         // Artificial cells
    icells.resize(100);                                           // 100 artificial target cells
    jcells.resize(100);                                           // 100 artificial source cells
    CI0 = icells.begin();                                         // Set global begin iterator for source
    for( C_iter Ci=icells.begin(); Ci!=icells.end(); ++Ci ) {     // Loop over target cells
      Ci->X = 0;                                                  //  Set coordinates of target cell
      Ci->NLEAF = NCRIT;                                          //  Number of leafs in target cell
      Ci->LEAF = ibodies.begin();                                 //  Leaf iterator in target cell
    }                                                             // End loop over target cells
    for( C_iter Cj=jcells.begin(); Cj!=jcells.end(); ++Cj ) {     // Loop over source cells
      Cj->X = 1;                                                  //  Set coordinates of source cell
      Cj->NLEAF = NCRIT;                                          //  Number of leafs in source cell
      Cj->LEAF = jbodies.begin();                                 //  Leaf iterator in source cell
    }                                                             // End loop over source cells
    listM2L.resize(icells.size());                                // Resize M2L interaction list
    listM2P.resize(icells.size());                                // Resize M2P interaction list
    listP2P.resize(icells.size());                                // Resize P2P interaction list
    for( C_iter Ci=icells.begin(); Ci!=icells.end(); ++Ci ) {     // Loop over target cells
      for( C_iter Cj=jcells.begin(); Cj!=jcells.end(); ++Cj ) {   //  Loop over source cells
        listP2P[Ci-CI0].push_back(Cj);                            //   Push source cell into P2P interaction list
        listM2P[Ci-CI0].push_back(Cj);                            //   Push source cell into P2P interaction list
        listM2L[Ci-CI0].push_back(Cj);                            //   Push source cell into P2P interaction list
      }                                                           //  End loop over source cells
    }                                                             // End loop over target cells
    startTimer("P2P kernel   ");                                  // Start timer
    evalP2P(icells);                                              // Evaluate P2P kernel
    timeP2P = stopTimer("P2P kernel   ") / NCRIT / NCRIT;         // Stop timer
    startTimer("M2L kernel   ");                                  // Start timer
    evalM2L(icells);                                              // Evaluate M2L kernel
    timeM2L = stopTimer("M2L kernel   ");                         // Stop timer
    startTimer("M2P kernel   ");                                  // Start timer
    evalM2P(icells);                                              // Evaluate M2P kernel
    timeM2P = stopTimer("M2P kernel   ") / NCRIT;                 // Stop timer
  }

  void evalP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU=false) {// Evaluate P2P
    int numIcall = int(ibodies.size()-1)/MAXBODY+1;               // Number of icall loops
    int numJcall = int(jbodies.size()-1)/MAXBODY+1;               // Number of jcall loops
    int ioffset = 0;                                              // Initialzie offset for icall loops
    for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
      BI0 = ibodies.begin()+ioffset;                              //  Set target bodies begin iterator
      BIN = ibodies.begin()+std::min(ioffset+MAXBODY,int(ibodies.size()));// Set target bodies end iterator
      int joffset = 0;                                            //  Initialize offset for jcall loops
      for( int jcall=0; jcall!=numJcall; ++jcall ) {              //  Loop over jcall
        BJ0 = jbodies.begin()+joffset;                            //  Set source bodies begin iterator
        BJN = jbodies.begin()+std::min(joffset+MAXBODY,int(jbodies.size()));// Set source bodies end iterator
        constHost.push_back(2*R0);                              //   Copy domain size to GPU buffer
        for( B_iter B=BJ0; B!=BJN; ++B ) {                      //   Loop over source bodies
          sourceHost.push_back(B->X[0]);                        //   Copy x position to GPU buffer
          sourceHost.push_back(B->X[1]);                        //   Copy y position to GPU buffer
          sourceHost.push_back(B->X[2]);                        //   Copy z position to GPU buffer
          sourceHost.push_back(B->SRC[0]);                      //   Copy 1st source value to GPU buffer
          sourceHost.push_back(B->SRC[1]);                      //   Copy 2nd source value to GPU buffer
          sourceHost.push_back(B->SRC[2]);                      //   Copy 3rd source value to GPU buffer
          sourceHost.push_back(B->SRC[3]);                      //   Copy 4th source value to GPU buffer
        }                                                       //   End loop over source bodies
        int key = 0;                                            //   Initialize key to range of leafs in source cells
        int blocks = (BIN - BI0 - 1) / THREADS + 1;             //   Number of thread blocks needed for this target cell
        for( int i=0; i!=blocks; ++i ) {                        //   Loop over thread blocks
          keysHost.push_back(key);                              //    Save key to range of leafs in source cells
        }                                                       //   End loop over thread blocks
        rangeHost.push_back(1);                                 //   Save size of interaction list
        rangeHost.push_back(0);                                 //   Set begin index of leafs
        rangeHost.push_back(BJN-BJ0);                           //   Set number of leafs
        rangeHost.push_back(0);                           //   Set periodic image flag
        for( B_iter B=BI0; B!=BIN; ++B ) {                      //   Loop over target bodies
          targetHost.push_back(B->X[0]);                        //    Copy x position to GPU buffer
          targetHost.push_back(B->X[1]);                        //    Copy y position to GPU buffer
          targetHost.push_back(B->X[2]);                        //    Copy z position to GPU buffer
          targetHost.push_back(B->SRC[0]);                      //    Copy 1st target value to GPU buffer
          targetHost.push_back(B->SRC[1]);                      //    Copy 2nd target value to GPU buffer
          targetHost.push_back(B->SRC[2]);                      //    Copy 3rd target value to GPU buffer
        }                                                       //   End loop over target bodies
        int numPad = blocks * THREADS - (BIN - BI0);            //   Number of elements to pad in target GPU buffer
        for( int i=0; i!=numPad; ++i ) {                        //   Loop over elements to pad
          targetHost.push_back(0);                              //    Pad x position in GPU buffer
          targetHost.push_back(0);                              //    Pad y position in GPU buffer
          targetHost.push_back(0);                              //    Pad z position in GPU buffer
          targetHost.push_back(0);                              //    Pad 1st target value to GPU buffer
          targetHost.push_back(0);                              //    Pad 2nd target value to GPU buffer
          targetHost.push_back(0);                              //    Pad 3rd target value to GPU buffer
        }                                                       //   End loop over elements to pad
        selectP2P();                                            //   Select P2P kernel
        for( B_iter B=BI0; B!=BIN; ++B ) {                      //    Loop over target bodies
          B->TRG[0] += targetHost[6*(B-BI0)+0];                 //     Copy 1st target value from GPU buffer
          B->TRG[1] += targetHost[6*(B-BI0)+1];                 //     Copy 2nd target value from GPU buffer
          B->TRG[2] += targetHost[6*(B-BI0)+2];                 //     Copy 3rd target value from GPU buffer
          B->TRG[3] += targetHost[6*(B-BI0)+3];                 //     Copy 4th target value from GPU buffer
        }                                                       //   Endif for Gaussian kernel
        keysHost.clear();                                       //   Clear keys vector
        rangeHost.clear();                                      //   Clear range vector
        constHost.clear();                                      //   Clear const vector
        targetHost.clear();                                     //   Clear target vector
        sourceHost.clear();                                     //   Clear source vector
        joffset += MAXBODY;                                       //  Increment jcall offset
      }                                                           // End loop over jcall
      ioffset += MAXBODY;                                         // Increment icall offset
    }                                                             // End loop over icall
  }

  void evalP2M(Cells &cells) {                         // Evaluate P2M
    CI0 = cells.begin();                                          // Set begin iterator for target
    const int numCell = MAXCELL/NCRIT/7;                          // Number of cells per icall
    int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
    int ioffset = 0;                                              // Initialzie offset for icall loops
    for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
      CIB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
      CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
      constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
      startTimer("Get list     ");                                //  Start timer
      Lists listP2M(cells.size());                                //  Define P2M interation list vector
      for( CI=CIB; CI!=CIE; ++CI ) {                              //  Loop over target cells
        CI->M = CI->L = 0;                                        //   Initialize multipole & local coefficients
        if( CI->NCHILD == 0 ) {                                   //   If cell is a twig
          listP2M[CI-CI0].push_back(CI);                          //    Push source cell into P2M interaction list
          sourceSize[CI] = CI->NLEAF;                             //    Key : iterator, Value : number of leafs
        }                                                         //   End loop over cells topdown
      }                                                           //  End loop over source map
      stopTimer("Get list     ");                                 //  Stop timer
      setSourceBody();                                            //  Set source buffer for bodies
      setTargetCell(listP2M);                             //  Set target buffer for cells
      selectP2M();                                                //  Select P2M kernel
      getTargetCell(listP2M);                                     //  Get body values from target buffer
      clearBuffers();                                             //  Clear GPU buffers
      ioffset += numCell;                                         //  Increment ioffset
    }                                                             // End loop over icall
  }

  void evalM2M(Cells &cells) {                         // Evaluate M2M
    CI0 = cells.begin();                                          // Set begin iterator for target
    const int numCell = MAXCELL/NCOEF/2;                          // Number of cells per icall
    int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
    int level = getLevel(CI0->ICELL);                             // Level of twig
    while( level != -1 ) {                                        // While level of target is not past root level
      int ioffset = 0;                                            //  Initialzie offset for icall loops
      for( int icall=0; icall!=numIcall; ++icall ) {              //  Loop over icall
        CIB = cells.begin()+ioffset;                              //   Set begin iterator for target per call
        CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
        constHost.push_back(2*R0);                                //   Copy domain size to GPU buffer
        startTimer("Get list     ");                              //   Start timer
        Lists listM2M(cells.size());                              //   Define M2M interation list vector
        for( CI=CIB; CI!=CIE; ++CI ) {                            //   Loop over cells bottomup (except root cell)
          if( getLevel(CI->ICELL) == level ) {                    //    If target cell is at current level
            for( int i=0; i<CI->NCHILD; ++i ) {                   //     Loop over child cells
              CJ = CI0 + CI->CHILD[i];                            //      Set iterator for source cell
              listM2M[CI-CI0].push_back(CJ);                      //      Push source cell into M2M interaction list
              sourceSize[CJ] = 2 * NCOEF;                         //      Key : iterator, Value : number of coefs
            }                                                     //     End loop over child cells
          }                                                       //    Endif for current level
        }                                                         //   End loop over cells
        stopTimer("Get list     ");                               //   Stop timer
        setSourceCell();                                          //   Set source buffer for cells
        setTargetCell(listM2M);                           //   Set target buffer for cells
        selectM2M();                                              //   Select M2M kernel
        getTargetCell(listM2M);                                   //   Get body values from target buffer
        clearBuffers();                                           //   Clear GPU buffers
        ioffset += numCell;                                       //   Increment ioffset
      }                                                           //  End loop over icall
      level--;                                                    //  Decrement level
    }                                                             // End while loop over levels
  }

  void evalM2L(Cells &cells, bool kernel=false) {            // Evaluate M2L
    CI0 = cells.begin();                                          // Set begin iterator
    const int numCell = MAXCELL/NCOEF/2;                          // Number of cells per icall
    int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
    int ioffset = 0 * kernel;                                     // Initialzie offset for icall loops
    for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
      CIB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
      CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
      constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
      startTimer("Get list     ");                                //  Start timer
      for( CI=CIB; CI!=CIE; ++CI ) {                              //  Loop over target cells
        for( L_iter L=listM2L[CI-CI0].begin(); L!=listM2L[CI-CI0].end(); ++L ) {//  Loop over interaction list
          CJ = *L;                                                //    Set source cell
          sourceSize[CJ] = 2 * NCOEF;                             //    Key : iterator, Value : number of coefs
        }                                                         //   End loop over interaction list
      }                                                           //  End loop over target cells
      stopTimer("Get list     ");                                 //  Stop timer
      setSourceCell();                                            //  Set source buffer for cells
      setTargetCell(listM2L);                             //  Set target buffer for cells
      selectM2L();                                                //  Select M2L kernel
      getTargetCell(listM2L,false);                               //  Get body values from target buffer
      clearBuffers();                                             //  Clear GPU buffers
      ioffset += numCell;                                         //  Increment ioffset
    }                                                             // End loop over icall
    listM2L.clear();                                              // Clear interaction lists
  }

  void evalM2P(Cells &cells, bool kernel=false) {            // Evaluate M2P
    CI0 = cells.begin();                                          // Set begin iterator for target
    const int numCell = MAXCELL/NCRIT/7;                          // Number of cells per icall
    int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
    int ioffset = 0 * kernel;                                     // Initialzie offset for icall loops
    for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
      CIB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
      CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
      constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
      startTimer("Get list     ");                                //  Start timer
      for( CI=CIB; CI!=CIE; ++CI ) {                              //  Loop over target cells
        for( L_iter L=listM2P[CI-CI0].begin(); L!=listM2P[CI-CI0].end(); ++L ) {//  Loop over interaction list
          CJ = *L;                                                //    Set source cell
          sourceSize[CJ] = 2 * NCOEF;                             //    Key : iterator, Value : number of coefs
        }                                                         //   End loop over interaction list
      }                                                           //  End loop over target cells
      stopTimer("Get list     ");                                 //  Stop timer
      setSourceCell();                                            //  Set source buffer for cells
      setTargetBody(listM2P);                             //  Set target buffer for bodies
      selectM2P();                                                //  Select M2P kernel
      getTargetBody(listM2P);                                     //  Get body values from target buffer
      clearBuffers();                                             //  Clear GPU buffers
      ioffset += numCell;                                         //  Increment ioffset
    }                                                             // End loop over icall
    listM2P.clear();                                              // Clear interaction lists
  }

  void evalP2P(Cells &cells, bool kernel=false) {            // Evaluate P2P
    CI0 = cells.begin();                                          // Set begin iterator
    const int numCell = MAXCELL/NCRIT/7;                          // Number of cells per icall
    int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
    int ioffset = 0 * kernel;                                     // Initialzie offset for icall loops
    for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
      CIB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
      CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
      constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
      startTimer("Get list     ");                                //  Start timer
      for( CI=CIB; CI!=CIE; ++CI ) {                              //  Loop over target cells
        for( L_iter L=listP2P[CI-CI0].begin(); L!=listP2P[CI-CI0].end(); ++L ) {//  Loop over interaction list
          CJ = *L;                                                //    Set source cell
          sourceSize[CJ] = CJ->NLEAF;                             //    Key : iterator, Value : number of leafs
        }                                                         //   End loop over interaction list
      }                                                           //  End loop over target cells
      stopTimer("Get list     ");                                 //  Stop timer
      setSourceBody();                                            //  Set source buffer for bodies
      setTargetBody(listP2P);                             //  Set target buffer for bodies
      selectP2P();                                                //  Select P2P kernel
      getTargetBody(listP2P);                                     //  Get body values from target buffer
      clearBuffers();                                             //  Clear GPU buffers
      ioffset += numCell;                                         //  Increment ioffset
    }                                                             // End loop over icall
    listP2P.clear();                                              // Clear interaction lists
  }

  void evalL2L(Cells &cells) {                         // Evaluate L2L
    CI0 = cells.begin();                                          // Set begin iterator
    const int numCell = MAXCELL/NCOEF/2;                          // Number of cells per icall
    int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
    int maxLevel = getLevel(CI0->ICELL);                          // Level of twig
    int level = 1;                                                // Start level from 1
    while( level != maxLevel+1 ) {                                // While level of source is not root level
      int ioffset = 0;                                            //  Initialzie offset for icall loops
      for( int icall=0; icall!=numIcall; ++icall ) {              //  Loop over icall
        CIB = cells.begin()+ioffset;                              //   Set begin iterator for target per call
        CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
        constHost.push_back(2*R0);                                //   Copy domain size to GPU buffer
        startTimer("Get list     ");                              //   Start timer
        Lists listL2L(cells.size());                              //   Define L2L interation list vector
        for( CI=CIE-2; CI!=CIB-1; --CI ) {                        //   Loop over cells topdown (except root cell)
          if( getLevel(CI->ICELL) == level ) {                    //    If target cell is at current level
            CJ = CI0 + CI->PARENT;                                //     Set source cell iterator
            listL2L[CI-CI0].push_back(CJ);                        //     Push source cell into L2L interaction list
            if( sourceSize[CJ] == 0 ) {                           //     If the source cell has not been stored yet
              sourceSize[CJ] = 2 * NCOEF;                         //      Key : iterator, Value : number of coefs
            }                                                     //     Endif for current level
          }                                                       //    Endif for stored source cell
        }                                                         //   End loop over cells topdown
        stopTimer("Get list     ");                               //   Stop timer
        setSourceCell(false);                                     //   Set source buffer for cells
        setTargetCell(listL2L);                           //   Set target buffer for cells
        selectL2L();                                              //   Select L2L kernel
        getTargetCell(listL2L,false);                             //   Get body values from target buffer
        clearBuffers();                                           //   Clear GPU buffers
        ioffset += numCell;                                       //   Increment ioffset
      }                                                           //  End loop over icall
      level++;                                                    //  Increment level
    }                                                             // End while loop over levels
  }

  void evalL2P(Cells &cells) {                         // Evaluate L2P
    CI0 = cells.begin();                                          // Set begin iterator
    const int numCell = MAXCELL/NCRIT/7;                          // Number of cells per icall
    int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
    int ioffset = 0;                                              // Initialzie offset for icall loops
    for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
      CIB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
      CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
      constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
      startTimer("Get list     ");                                //  Start timer
      Lists listL2P(cells.size());                                //  Define L2P interation list vector
      for( CI=CIB; CI!=CIE; ++CI ) {                              //  Loop over cells
        if( CI->NCHILD == 0 ) {                                   //   If cell is a twig evaluate L2P kernel
          listL2P[CI-CI0].push_back(CI);                          //    Push source cell into L2P interaction list
          sourceSize[CI] = 2 * NCOEF;                             //    Key : iterator, Value : number of coefs
        }                                                         //   Endif for twig cells
      }                                                           //  End loop over cells topdown
      stopTimer("Get list     ");                                 //  Stop timer
      setSourceCell(false);                                       //  Set source buffer for cells
      setTargetBody(listL2P);                             //  Set target buffer for bodies
      selectL2P();                                                //  Select L2P kernel
      getTargetBody(listL2P);                                     //  Get body values from target buffer
      clearBuffers();                                             //  Clear GPU buffers
      ioffset += numCell;                                         //  Increment ioffset
    }                                                             // End loop over icall
  }
};

#endif
