#ifndef dataset_h
#define dataset_h
#include "types.h"

//! Contains all the different datasets
class Dataset {
private:
  long filePosition;                                            //!< Position of file stream

public:
//! Constructor
  Dataset() : filePosition(0) {}
//! Destructor
  ~Dataset() {}

//! Initialize source values
  void initSource(Bodies &bodies) {
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->SRC = 0;                                               //  Clear previous source values
      B->SRC[0] = 1. / bodies.size() / MPISIZE;                 //   Initialize mass/charge
    }                                                           // End loop over bodies
  }

//! Initialize target values
  void initTarget(Bodies &bodies, bool IeqJ=true) {
    srand(1);                                                   // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->TRG = 0 * IeqJ;                                        //  Clear previous target values (IeqJ is dummy)
      B->TRG[0] = -B->SRC[0] / std::sqrt(EPS2) * IeqJ;          //   Initialize potential (0 if I != J)
    }                                                           // End loop over bodies
  }

//! Random distribution in [-1,1]^3 cube
  void random(Bodies &bodies, int seed=1, int numSplit=1) {
    srand(seed);                                                // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit) ) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand(seed);                                            //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;   //   Initialize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

//! Random distribution on r = 1 sphere
  void sphere(Bodies &bodies, int seed=1, int numSplit=1) {
    srand(seed);                                                // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit) ) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand(seed);                                            //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] = rand() / (1. + RAND_MAX) * 2 - 1;             //   Initialize positions
      }                                                         //  End loop over dimension
      real r = std::sqrt(norm(B->X));                           //  Distance from center
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] /= r * 1.1;                                     //   Normalize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

//! Uniform distribution on [-1,1]^3 lattice (for debugging)
  void lattice(Bodies &bodies) {
    int level = int(log(bodies.size()*MPISIZE+1.)/M_LN2/3);     // Level of tree
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int d = 0, l = 0;                                         //  Initialize dimension and level
      int index = MPIRANK * bodies.size() + (B-bodies.begin()); //  Set index of body iterator
      vec<3,int> nx = 0;                                        //  Initialize 3-D cell index
      while( index != 0 ) {                                     //  Deinterleave bits while index is nonzero
        nx[d] += (index % 2) * (1 << l);                        //   Add deinterleaved bit to 3-D cell index
        index >>= 1;                                            //   Right shift the bits
        d = (d+1) % 3;                                          //   Increment dimension
        if( d == 0 ) l++;                                       //   If dimension is 0 again, increment level
      }                                                         //  End while loop for deinterleaving bits
      for( d=0; d!=3; ++d ) {                                   //  Loop over dimensions
        B->X[d] = -1 + (2 * nx[d] + 1.) / (1 << level);         //   Calculate cell center from 3-D cell index
      }                                                         //  End loop over dimensions
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

  void sampleBodies(Bodies &bodies, int numTargets) {
    int n = bodies.size();
    int p = n / numTargets;
    assert(p > 0);
    for (int i=0; i<numTargets; i++) {
      assert(i * p < n);
      bodies[i] = bodies[i*p];
    }
    bodies.resize(numTargets);
  }

//! Evaluate relative L2 norm error
  void evalError(Bodies &bodies, Bodies &bodies2,
                 real &diff1, real &norm1, real &diff2, real &norm2) {
    B_iter B2 = bodies2.begin();                              //  Set iterator for bodies2
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
      std::cout << B->ICELL << " " << B->TRG[0] << " " << B2->TRG[0] << std::endl;// Compare every element
#endif
      diff1 += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);// Difference of potential
      norm1 += B2->TRG[0] * B2->TRG[0];                       //  Value of potential
      diff2 += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]);// Difference of x acceleration
      diff2 += (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]);// Difference of y acceleration
      diff2 += (B->TRG[3] - B2->TRG[3]) * (B->TRG[3] - B2->TRG[3]);// Difference of z acceleration
      norm2 += B2->TRG[1] * B2->TRG[1];                       //  Value of x acceleration
      norm2 += B2->TRG[2] * B2->TRG[2];                       //  Value of y acceleration
      norm2 += B2->TRG[3] * B2->TRG[3];                       //  Value of z acceleration
    }                                                         //  End loop over bodies & bodies2
  }

//! Print relative L2 norm error
  void printError(real diff1, real norm1, real diff2, real norm2) {
    std::cout << "Error (pot)   : " << std::sqrt(diff1/norm1) << std::endl;
    std::cout << "Error (acc)   : " << std::sqrt(diff2/norm2) << std::endl;
  }
};

#endif
