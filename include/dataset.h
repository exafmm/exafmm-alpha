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

  void initBodies(Bodies &bodies) {
    srand(0);
    int b = 0;
    while (b < int(bodies.size())) {
#if 1
      bodies[b].X[0] = drand48();
      bodies[b].X[1] = drand48();
      bodies[b].X[2] = drand48();
      b++;
#else
      double X1 = drand48();
      double X2 = drand48();
      double X3 = drand48();
      double R = 1.0 / sqrt( (pow(X1, -2.0 / 3.0) - 1.0) );
      if (R < 100.0) {
        double Z = (1.0 - 2.0 * X2) * R;
        double X = sqrt(R * R - Z * Z) * cos(2.0 * M_PI * X3);
        double Y = sqrt(R * R - Z * Z) * sin(2.0 * M_PI * X3);
        double scale = 3.0 * M_PI / 16.0;
        X *= scale; Y *= scale; Z *= scale;
        bodies[b].X[0] = X;
        bodies[b].X[1] = Y;
        bodies[b].X[2] = Z;
        b++;
      }
#endif
    }
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
