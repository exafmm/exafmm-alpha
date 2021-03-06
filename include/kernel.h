#ifndef kernel_h
#define kernel_h
#define KERNEL
#include "sort.h"
#undef KERNEL
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

const int  P2 = P * P;                                          //!< P^2
const int  P4 = P2 * P2;                                        //!< P^4
const real EPS = 1e-6;                                          //!< Single precision epsilon

//! Unified CPU/GPU kernel class
class Kernel : public Sort {
protected:
  B_iter      BI0;                                              //!< Target bodies begin iterator
  B_iter      BIN;                                              //!< Target bodies end iterator
  B_iter      BJ0;                                              //!< Source bodies begin iterator
  B_iter      BJN;                                              //!< Source bodies end iterator
  C_iter      CI;                                               //!< Target cell iterator
  C_iter      CJ;                                               //!< Source cell iterator
  vect        X0;                                               //!< Center of root cell
  real        R0;                                               //!< Radius of root cell
  vect        Xperiodic;                                        //!< Coordinate offset of periodic image
  KernelName  kernelName;                                       //!< Name of kernel

  int                  ATOMS;                                   //!< Number of atom types in Van der Waals
  std::vector<real>    RSCALE;                                  //!< Scaling parameter for Van der Waals
  std::vector<real>    GSCALE;                                  //!< Scaling parameter for Van der Waals

  std::vector<int>     keysHost;                                //!< Offsets for rangeHost
  std::vector<int>     rangeHost;                               //!< Offsets for sourceHost
  std::vector<gpureal> constHost;                               //!< Constants on host
  std::vector<gpureal> sourceHost;                              //!< Sources on host
  std::vector<gpureal> targetHost;                              //!< Targets on host
  Map                  sourceBegin;                             //!< Define map for offset of source cells
  Map                  sourceSize;                              //!< Define map for size of source cells
  Map                  targetBegin;                             //!< Define map for offset of target cells

  double *factorial;                                            //!< Factorial
  double *prefactor;                                            //!< \f$ \sqrt{ \frac{(n - |m|)!}{(n + |m|)!} } \f$
  double *Anm;                                                  //!< \f$ (-1)^n / \sqrt{ \frac{(n + m)!}{(n - m)!} } \f$
  complex *Ynm;                                                 //!< \f$ r^n Y_n^m \f$
  complex *YnmTheta;                                            //!< \f$ \theta \f$ derivative of \f$ r^n Y_n^m \f$
  complex *Cnm;                                                 //!< M2L translation matrix \f$ C_{jn}^{km} \f$
public:
  int IMAGES;                                                   //!< Number of periodic image sublevels
  real THETA;                                                   //!< Box opening criteria
  real NP2P;                                                    //!< Number of P2P kernel call
  real NM2P;                                                    //!< Number of M2P kernel call
  real NM2L;                                                    //!< Number of M2L kernel call

private:
//! Get r,theta,phi from x,y,z
  void cart2sph(real& r, real& theta, real& phi, vect dist) {
    r = std::sqrt(norm(dist))+EPS;                              // r = sqrt(x^2 + y^2 + z^2) + eps
    theta = std::acos(dist[2] / r);                             // theta = acos(z / r)
    if( std::abs(dist[0]) + std::abs(dist[1]) < EPS ) {         // If |x| < eps & |y| < eps
      phi = 0;                                                  //  phi can be anything so we set it to 0
    } else if( std::abs(dist[0]) < EPS ) {                      // If |x| < eps
      phi = dist[1] / std::abs(dist[1]) * M_PI * 0.5;           //  phi = sign(y) * pi / 2
    } else if( dist[0] > 0 ) {                                  // If x > 0
      phi = std::atan(dist[1] / dist[0]);                       //  phi = atan(y / x)
    } else {                                                    // If x < 0
      phi = std::atan(dist[1] / dist[0]) + M_PI;                //  phi = atan(y / x) + pi
    }                                                           // End if for x,y cases
  }

//! Spherical to cartesian coordinates
  template<typename T>
  void sph2cart(real r, real theta, real phi, T spherical, T &cartesian) {
    cartesian[0] = sin(theta) * cos(phi) * spherical[0]         // x component (not x itself)
                 + cos(theta) * cos(phi) / r * spherical[1]
                 - sin(phi) / r / sin(theta) * spherical[2];
    cartesian[1] = sin(theta) * sin(phi) * spherical[0]         // y component (not y itself)
                 + cos(theta) * sin(phi) / r * spherical[1]
                 + cos(phi) / r / sin(theta) * spherical[2];
    cartesian[2] = cos(theta) * spherical[0]                    // z component (not z itself)
                 - sin(theta) / r * spherical[1];
  }

//! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
  void evalMultipole(real rho, real alpha, real beta) {
    const complex I(0.,1.);                                     // Imaginary unit
    double x = std::cos(alpha);                                 // x = cos(alpha)
    double y = std::sin(alpha);                                 // y = sin(alpha)
    double fact = 1;                                            // Initialize 2 * m + 1
    double pn = 1;                                              // Initialize Legendre polynomial Pn
    double rhom = 1;                                            // Initialize rho^m
    for( int m=0; m!=P; ++m ) {                                 // Loop over m in Ynm
      complex eim = std::exp(I * double(m * beta));             //  exp(i * m * beta)
      double p = pn;                                            //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * prefactor[npn] * eim;               //  rho^m * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      double p1 = p;                                            //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;// theta derivative of r^n * Ynm
      rhom *= rho;                                              //  rho^m
      double rhon = rhom;                                       //  rho^n
      for( int n=m+1; n!=P; ++n ) {                             //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * prefactor[npm] * eim;             //   rho^n * Ynm
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        double p2 = p1;                                         //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;// theta derivative
        rhon *= rho;                                            //   Update rho^n
      }                                                         //  End loop over n in Ynm
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
    }                                                           // End loop over m in Ynm
  }

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
  void evalLocal(real rho, real alpha, real beta) {
    const complex I(0.,1.);                                     // Imaginary unit
    double x = std::cos(alpha);                                 // x = cos(alpha)
    double y = std::sin(alpha);                                 // y = sin(alpha)
    double fact = 1;                                            // Initialize 2 * m + 1
    double pn = 1;                                              // Initialize Legendre polynomial Pn
    double rhom = 1.0 / rho;                                    // Initialize rho^(-m-1)
    for( int m=0; m!=2*P; ++m ) {                               // Loop over m in Ynm
      complex eim = std::exp(I * double(m * beta));             //  exp(i * m * beta)
      double p = pn;                                            //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * prefactor[npn] * eim;               //  rho^(-m-1) * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      double p1 = p;                                            //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;// theta derivative of r^n * Ynm
      rhom /= rho;                                              //  rho^(-m-1)
      double rhon = rhom;                                       //  rho^(-n-1)
      for( int n=m+1; n!=2*P; ++n ) {                           //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * prefactor[npm] * eim;             //   rho^n * Ynm for m > 0
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        double p2 = p1;                                         //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;// theta derivative
        rhon /= rho;                                            //   rho^(-n-1)
      }                                                         //  End loop over n in Ynm
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
    }                                                           // End loop over m in Ynm
  }

protected:
//! Get level from cell index
  int getLevel(bigint index) {
    int level = -1;                                             // Initialize level counter
    while( index >= 0 ) {                                       // While cell index is non-negative
      level++;                                                  //  Increment level
      index -= 1 << 3*level;                                    //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

public:
//! Constructor
  Kernel() : X0(0), R0(0), NP2P(0), NM2P(0), NM2L(0) {}
//! Destructor
  ~Kernel() {}

//! Get center of root cell
  vect getX0() {return X0;}
//! Get radius of root cell
  real getR0() {return R0;}

//! Set center and size of root cell
  void setDomain(Bodies &bodies, vect x0=0, real r0=M_PI) {
    vect xmin,xmax;                                             // Min,Max of domain
    B_iter B = bodies.begin();                                  // Reset body iterator
    xmin = xmax = B->X;                                         // Initialize xmin,xmax
    X0 = 0;                                                     // Initialize center and size of root cell
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->X[d] < xmin[d]) xmin[d] = B->X[d];           //   Determine xmin
        else if(B->X[d] > xmax[d]) xmax[d] = B->X[d];           //   Determine xmax
      }                                                         //  End loop over each dimension
      X0 += B->X;                                               //  Sum positions
    }                                                           // End loop over bodies
    X0 /= bodies.size();                                        // Calculate average position
    for( int d=0; d!=3; ++d ) {                                 // Loop over each dimension
      X0[d] = int(X0[d]+.5);                                    //  Shift center to nearest integer
      R0 = std::max(xmax[d] - X0[d], R0);                       //  Calculate max distance from center
      R0 = std::max(X0[d] - xmin[d], R0);                       //  Calculate max distance from center
    }                                                           // End loop over each dimension
    R0 += 1e-5;                                                 // Add some leeway to root radius
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      if( X0[0]-R0 < x0[0]-r0 || x0[0]+r0 < X0[0]+R0            //  Check for outliers in x direction
       || X0[1]-R0 < x0[1]-r0 || x0[1]+r0 < X0[1]+R0            //  Check for outliers in y direction
       || X0[2]-R0 < x0[2]-r0 || x0[2]+r0 < X0[2]+R0 ) {        //  Check for outliers in z direction
        std::cout << "Error: Particles located outside periodic domain" << std::endl;// Print error message
      }                                                         //  End if for outlier checking
      X0 = x0;                                                  //  Center is [0, 0, 0]
      R0 = r0;                                                  //  Radius is r0
    }                                                           // Endif for periodic boundary condition
  }

//! Set scaling paramters in Van der Waals
  void setVanDerWaals(int atoms, double *rscale, double *gscale) {
    ATOMS = atoms;                                              // Set number of atom types
    RSCALE.resize(ATOMS*ATOMS);                                 // Resize rscale vector
    GSCALE.resize(ATOMS*ATOMS);                                 // Resize gscale vector
    for( int i=0; i!=ATOMS*ATOMS; ++i ) {                       // Loop over scale vector
      RSCALE[i] = rscale[i];                                    //  Set rscale vector
      GSCALE[i] = gscale[i];                                    //  Set gscale vector
    }                                                           // End loop over scale vector
  }

//! Precalculate M2L translation matrix
  void preCalculation() {
    const complex I(0.,1.);                                     // Imaginary unit
    factorial = new double  [P];                                // Factorial
    prefactor = new double  [4*P2];                             // sqrt( (n - |m|)! / (n + |m|)! )
    Anm       = new double  [4*P2];                             // (-1)^n / sqrt( (n + m)! / (n - m)! )
    Ynm       = new complex [4*P2];                             // r^n * Ynm
    YnmTheta  = new complex [4*P2];                             // theta derivative of r^n * Ynm
    Cnm       = new complex [P4];                               // M2L translation matrix Cjknm

    factorial[0] = 1;                                           // Initialize factorial
    for( int n=1; n!=P; ++n ) {                                 // Loop to P
      factorial[n] = factorial[n-1] * n;                        //  n!
    }                                                           // End loop to P

    for( int n=0; n!=2*P; ++n ) {                               // Loop over n in Anm
      for( int m=-n; m<=n; ++m ) {                              //  Loop over m in Anm
        int nm = n*n+n+m;                                       //   Index of Anm
        int nabsm = abs(m);                                     //   |m|
        double fnmm = 1.0;                                      //   Initialize (n - m)!
        for( int i=1; i<=n-m; ++i ) fnmm *= i;                  //   (n - m)!
        double fnpm = 1.0;                                      //   Initialize (n + m)!
        for( int i=1; i<=n+m; ++i ) fnpm *= i;                  //   (n + m)!
        double fnma = 1.0;                                      //   Initialize (n - |m|)!
        for( int i=1; i<=n-nabsm; ++i ) fnma *= i;              //   (n - |m|)!
        double fnpa = 1.0;                                      //   Initialize (n + |m|)!
        for( int i=1; i<=n+nabsm; ++i ) fnpa *= i;              //   (n + |m|)!
        prefactor[nm] = std::sqrt(fnma/fnpa);                   //   sqrt( (n - |m|)! / (n + |m|)! )
        Anm[nm] = ODDEVEN(n)/std::sqrt(fnmm*fnpm);              //   (-1)^n / sqrt( (n + m)! / (n - m)! )
      }                                                         //  End loop over m in Anm
    }                                                           // End loop over n in Anm

    for( int j=0, jk=0, jknm=0; j!=P; ++j ) {                   // Loop over j in Cjknm
      for( int k=-j; k<=j; ++k, ++jk ){                         //  Loop over k in Cjknm
        for( int n=0, nm=0; n!=P; ++n ) {                       //   Loop over n in Cjknm
          for( int m=-n; m<=n; ++m, ++nm, ++jknm ) {            //    Loop over m in Cjknm
            const int jnkm = (j+n)*(j+n)+j+n+m-k;               //     Index C_{j+n}^{m-k}
            Cnm[jknm] = std::pow(I,double(abs(k-m)-abs(k)-abs(m)))*(ODDEVEN(j)*Anm[nm]*Anm[jk]/Anm[jnkm]);// Cjknm
          }                                                     //    End loop over m in Cjknm
        }                                                       //   End loop over n in Cjknm
      }                                                         //  End loop over in k in Cjknm
    }                                                           // End loop over in j in Cjknm
  }

//! Free temporary allocations
  void postCalculation() {
    delete[] factorial;                                         // Free factorial
    delete[] prefactor;                                         // Free sqrt( (n - |m|)! / (n + |m|)! )
    delete[] Anm;                                               // Free (-1)^n / sqrt( (n + m)! / (n - m)! )
    delete[] Ynm;                                               // Free r^n * Ynm
    delete[] YnmTheta;                                          // Free theta derivative of r^n * Ynm
    delete[] Cnm;                                               // Free M2L translation matrix Cjknm
  }

  void LaplaceInit();                                           //!< Initialize Laplace kernels
  void LaplaceP2M();                                            //!< Evaluate Laplace P2M kernel
  void LaplaceM2M();                                            //!< Evaluate Laplace M2M kernel
  void LaplaceM2M_CPU();                                        //!< Evaluate Laplace M2M kernel on CPU
  void LaplaceM2L();                                            //!< Evaluate Laplace M2L kernel
  void LaplaceM2P();                                            //!< Evaluate Laplace M2P kernel
  void LaplaceP2P();                                            //!< Evaluate Laplace P2P kernel
  void LaplaceP2P_CPU();                                        //!< Evaluate Laplace P2P kernel on CPU
  void LaplaceL2L();                                            //!< Evaluate Laplace L2L kernel
  void LaplaceL2P();                                            //!< Evaluate Laplace L2P kernel
  void LaplaceFinal();                                          //!< Finalize Lapalce kernels

  void BiotSavartInit();                                        //!< Initialize Biot-Savart kernels
  void BiotSavartP2M();                                         //!< Evaluate Biot-Savart P2M kernel
  void BiotSavartM2M();                                         //!< Evalaute Biot-Savart M2M kernel
  void BiotSavartM2M_CPU();                                     //!< Evaluate Biot-Savart M2M kernel on CPU
  void BiotSavartM2L();                                         //!< Evaluate Biot-Savart M2L kernel
  void BiotSavartM2P();                                         //!< Evaluate Biot-Savart M2P kernel
  void BiotSavartP2P();                                         //!< Evaluate Biot-Savart P2P kernel
  void BiotSavartP2P_CPU();                                     //!< Evaluate Biot-Savart P2P kernel on CPU
  void BiotSavartL2L();                                         //!< Evaluate Biot-Savart L2L kernel
  void BiotSavartL2P();                                         //!< Evaluate Biot-Savart L2P kernel
  void BiotSavartFinal();                                       //!< Finalize Biot-Savart kernels

  void StretchingInit();                                        //!< Initialize Stretching kernels
  void StretchingP2M();                                         //!< Evaluate Stretching P2M kernel
  void StretchingM2M();                                         //!< Evaluate Stretching M2M kernel
  void StretchingM2M_CPU();                                     //!< Evaluate Stretching M2M kernel on CPU
  void StretchingM2L();                                         //!< Evaluate Stretching M2L kernel
  void StretchingM2P();                                         //!< Evaluate Stretching M2P kernel
  void StretchingP2P();                                         //!< Evaluate Stretching P2P kernel
  void StretchingP2P_CPU();                                     //!< Evaluate Stretching P2P kernel on CPU
  void StretchingL2L();                                         //!< Evaluate Stretching L2L kernel
  void StretchingL2P();                                         //!< Evaluate Stretching L2P kernel
  void StretchingFinal();                                       //!< Finalize Stretching kernels

  void GaussianInit();                                          //!< Initialize Gaussian kernels
  void GaussianP2M();                                           //!< Dummy
  void GaussianM2M();                                           //!< Dummy
  void GaussianM2M_CPU();                                       //!< Dummy
  void GaussianM2L();                                           //!< Dummy
  void GaussianM2P();                                           //!< Dummy
  void GaussianP2P();                                           //!< Evaluate Gaussian P2P kernel
  void GaussianP2P_CPU();                                       //!< Evaluate Gaussian P2P kernel on CPU
  void GaussianL2L();                                           //!< Dummy
  void GaussianL2P();                                           //!< Dummy
  void GaussianFinal();                                         //!< Finalize Gaussian kernels

  void CoulombVdWInit();                                        //!< Initialize CoulombVdW kernels
  void CoulombVdWP2M();                                         //!< Evaluate CoulombVdW P2M kernel
  void CoulombVdWM2M();                                         //!< Evaluate CoulombVdW M2M kernel
  void CoulombVdWM2M_CPU();                                     //!< Evaluate CoulombVdW M2M kernel on CPU
  void CoulombVdWM2L();                                         //!< Evaluate CoulombVdW M2L kernel
  void CoulombVdWM2P();                                         //!< Evaluate CoulombVdW M2P kernel
  void CoulombVdWP2P();                                         //!< Evaluate CoulombVdW P2P kernel
  void CoulombVdWP2P_CPU();                                     //!< Evaluate CoulombVdW P2P kernel on CPU
  void CoulombVdWL2L();                                         //!< Evaluate CoulombVdW L2L kernel
  void CoulombVdWL2P();                                         //!< Evaluate CoulombVdW L2P kernel
  void CoulombVdWFinal();                                       //!< Finalize CoulombVdW kernels

//! Select P2P kernel
  void selectP2P() {
    if( kernelName == Laplace ) {                               // If Laplace kernel
      LaplaceP2P();                                             //  Evaluate P2P kernel
    } else if ( kernelName == BiotSavart ) {                    // If Biot Savart kernel
      BiotSavartP2P();                                          //  Evaluate P2P kernel
    } else if ( kernelName == Stretching ) {                    // If Stretching kernel
      StretchingP2P();                                          //  Evaluate P2P kernel
    } else if ( kernelName == Gaussian ) {                      // If Gaussian kernel
      GaussianP2P();                                            //  Evaluate P2P kernel
    } else if ( kernelName == CoulombVdW ) {                    // If CoulombVdW kernel
      CoulombVdWP2P();                                          //  Evaluate P2P kernel
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }

//! Select P2P_CPU kernel
  void selectP2P_CPU() {
    if( kernelName == Laplace ) {                               // If Laplace kernel
      LaplaceP2P_CPU();                                         //  Evaluate P2P_CPU kernel
    } else if ( kernelName == BiotSavart ) {                    // If Biot Savart kernel
      BiotSavartP2P_CPU();                                      //  Evaluate P2P_CPU kernel
    } else if ( kernelName == Stretching ) {                    // If Stretching kernel
      StretchingP2P_CPU();                                      //  Evaluate P2P_CPU kernel
    } else if ( kernelName == Gaussian ) {                      // If Gaussian kernel
      GaussianP2P_CPU();                                        //  Evaluate P2P_CPU kernel
    } else if ( kernelName == CoulombVdW ) {                    // If CoulombVdW kernel
      CoulombVdWP2P_CPU();                                      //  Evaluate P2P_CPU kernel
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }

//! Select P2M kernel
  void selectP2M() {
    if( kernelName == Laplace ) {                               // If Laplace kernel
      LaplaceP2M();                                             //  Evaluate P2M kernel
    } else if ( kernelName == BiotSavart ) {                    // If Biot Savart kernel
      BiotSavartP2M();                                          //  Evaluate P2M kernel
    } else if ( kernelName == Stretching ) {                    // If Stretching kernel
      StretchingP2M();                                          //  Evaluate P2M kernel
    } else if ( kernelName == Gaussian ) {                      // If Gaussian kernel
      GaussianP2M();                                            //  Evaluate P2M kernel
    } else if ( kernelName == CoulombVdW ) {                    // If CoulombVdW kernel
      CoulombVdWP2M();                                          //  Evaluate P2M kernel
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }

//! Select M2M kernel
  void selectM2M() {
    if( kernelName == Laplace ) {                               // If Laplace kernel
      LaplaceM2M();                                             //  Evaluate M2M kernel
    } else if ( kernelName == BiotSavart ) {                    // If Biot Savart kernel
      BiotSavartM2M();                                          //  Evaluate M2M kernel
    } else if ( kernelName == Stretching ) {                    // If Stretching kernel
      StretchingM2M();                                          //  Evaluate M2M kernel
    } else if ( kernelName == Gaussian ) {                      // If Gaussian kernel
      GaussianM2M();                                            //  Evaluate M2M kernel
    } else if ( kernelName == CoulombVdW ) {                    // If CoulombVdW kernel
      CoulombVdWM2M();                                          //  Evaluate M2M kernel
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }

//! Select M2M_CPU kernel
  void selectM2M_CPU() {
    if( kernelName == Laplace ) {                               // If Laplace kernel
      LaplaceM2M_CPU();                                         //  Evaluate M2M_CPU kernel
    } else if ( kernelName == BiotSavart ) {                    // If Biot Savart kernel
      BiotSavartM2M_CPU();                                      //  Evaluate M2M_CPU kernel
    } else if ( kernelName == Stretching ) {                    // If Stretching kernel
      StretchingM2M_CPU();                                      //  Evaluate M2M_CPU kernel
    } else if ( kernelName == Gaussian ) {                      // If Gaussian kernel
      GaussianM2M_CPU();                                        //  Evaluate M2M_CPU kernel
    } else if ( kernelName == CoulombVdW ) {                    // If CoulombVdW kernel
      CoulombVdWM2M_CPU();                                      //  Evaluate M2M_CPU kernel
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }

//! Select M2L kernel
  void selectM2L() {
    if( kernelName == Laplace ) {                               // If Laplace kernel
      LaplaceM2L();                                             //  Evaluate M2L kernel
    } else if ( kernelName == BiotSavart ) {                    // If Biot Savart kernel
      BiotSavartM2L();                                          //  Evaluate M2L kernel
    } else if ( kernelName == Stretching ) {                    // If Stretching kernel
      StretchingM2L();                                          //  Evaluate M2L kernel
    } else if ( kernelName == Gaussian ) {                      // If Gaussian kernel
      GaussianM2L();                                            //  Evaluate M2L kernel
    } else if ( kernelName == CoulombVdW ) {                    // If CoulombVdW kernel
      CoulombVdWM2L();                                          //  Evaluate M2L kernel
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }

//! Select M2P kernel
  void selectM2P() {
    if( kernelName == Laplace ) {                               // If Laplace kernel
      LaplaceM2P();                                             //  Evaluate M2P kernel
    } else if ( kernelName == BiotSavart ) {                    // If Biot Savart kernel
      BiotSavartM2P();                                          //  Evaluate M2P kernel
    } else if ( kernelName == Stretching ) {                    // If Stretching kernel
      StretchingM2P();                                          //  Evaluate M2P kernel
    } else if ( kernelName == Gaussian ) {                      // If Gaussian kernel
      GaussianM2P();                                            //  Evaluate M2P kernel
    } else if ( kernelName == CoulombVdW ) {                    // If CoulombVdW kernel
      CoulombVdWM2P();                                          //  Evaluate M2P kernel
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }

//! Select L2L kernel
  void selectL2L() {
    if( kernelName == Laplace ) {                               // If Laplace kernel
      LaplaceL2L();                                             //  Evaluate L2L kernel
    } else if ( kernelName == BiotSavart ) {                    // If Biot Savart kernel
      BiotSavartL2L();                                          //  Evaluate L2L kernel
    } else if ( kernelName == Stretching ) {                    // If Stretching kernel
      StretchingL2L();                                          //  Evaluate L2L kernel
    } else if ( kernelName == Gaussian ) {                      // If Gaussian kernel
      GaussianL2L();                                            //  Evaluate L2L kernel
    } else if ( kernelName == CoulombVdW ) {                    // If CoulombVdW kernel
      CoulombVdWL2L();                                          //  Evaluate L2L kernel
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }

//! Select L2P kernel
  void selectL2P() {
    if( kernelName == Laplace ) {                               // If Laplace kernel
      LaplaceL2P();                                             //  Evaluate L2P kernel
    } else if ( kernelName == BiotSavart ) {                    // If Biot Savart kernel
      BiotSavartL2P();                                          //  Evaluate L2P kernel
    } else if ( kernelName == Stretching ) {                    // If Stretching kernel
      StretchingL2P();                                          //  Evaluate L2P kernel
    } else if ( kernelName == Gaussian ) {                      // If Gaussian kernel
      GaussianL2P();                                            //  Evaluate L2P kernel
    } else if ( kernelName == CoulombVdW ) {                    // If CoulombVdW kernel
      CoulombVdWL2P();                                          //  Evaluate L2P kernel
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }
};

#endif
