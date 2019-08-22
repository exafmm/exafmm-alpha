#ifndef types_h
#define types_h
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <omp.h>
#include <stack>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>
#include "vec.h"

typedef long                 bigint;                            //!< Big integer type
typedef double               real;                              //!< Real number type on CPU
typedef double               gpureal;                           //!< Real number type on GPU
typedef std::complex<double> complex;                           //!< Complex number type

#ifndef KERNEL
int MPIRANK = 0;                                                //!< MPI comm rank
int MPISIZE = 1;                                                //!< MPI comm size
#else
extern int MPIRANK;                                             //!< MPI comm rank
extern int MPISIZE;                                             //!< MPI comm size
#endif

const int  P       = 10;                                        //!< Order of expansions
const int  NCRIT   = 1000;                                      //!< Number of bodies per cell
const int  MAXBODY = 200000;                                    //!< Maximum number of bodies per GPU kernel
const int  MAXCELL = 10000000;                                  //!< Maximum number of bodies/coefs in cell per GPU kernel
const real CLET    = 2;                                         //!< LET opening critetia
const real EPS2    = 1e-8;                                      //!< Softening parameter
const int  GPUS    = 4;                                         //!< Number of GPUs per node
const int  THREADS = 64;                                        //!< Number of threads per thread-block

const int  NTERM   = P*(P+1)/2;                                 //!< Number of terms for spherical harmonics
const int  NCOEF   = 3 * NTERM;                                 //!< 3-D vector of coefficients

typedef vec<3,real>                            vect;            //!< 3-D vector type
typedef vec<NCOEF,complex>                     coef;            //!< Multipole coefficient type for spherical
typedef std::vector<bigint>                    Bigints;         //!< Vector of big integer types
typedef std::map<std::string,double>           Event;           //!< Map of event name to logged value
typedef std::map<std::string,double>::iterator E_iter;          //!< Iterator for event name map

enum KernelName {Laplace};                                      //!< Kernel name enumeration

//! Structure for source bodies (stuff to send)
struct JBody {
  int         IBODY;                                            //!< Initial body numbering for sorting back
  int         IPROC;                                            //!< Initial process numbering for partitioning back
  bigint      ICELL;                                            //!< Cell index
  vect        X;                                                //!< Position
  vec<4,real> SRC;                                              //!< Source values
  vect        dxdt;
  vect        dxdt2;
  vect        dgdt;
};
//! Structure for bodies
struct Body : JBody {
  vec<4,real> TRG;                                              //!< Target values
  bool operator<(const Body &rhs) const {                       //!< Overload operator for comparing body index
    return this->IBODY < rhs.IBODY;                             //!< Comparison function for body index
  }
};
typedef std::vector<Body>              Bodies;                  //!< Vector of bodies
typedef std::vector<Body>::iterator    B_iter;                  //!< Iterator for body vector
typedef std::vector<JBody>             JBodies;                 //!< Vector of source bodies
typedef std::vector<JBody>::iterator   JB_iter;                 //!< Iterator for source body vector

//! Structure for source cells (stuff to send)
struct JCell {
  bigint ICELL;                                                 //!< Cell index
  coef   M;                                                     //!< Multipole coefficients
};
//! Structure for cells
struct Cell : JCell {
  int    NCHILD;                                                //!< Number of child cells
  int    NLEAF;                                                 //!< Number of leafs
  int    PARENT;                                                //!< Iterator offset of parent cell
  int    CHILD;                                                 //!< Iterator offset of child cells
  B_iter LEAF;                                                  //!< Iterator of first leaf
  vect   X;                                                     //!< Cell center
  real   R;                                                     //!< Cell radius
  coef   L;                                                     //!< Local coefficients
};
typedef std::vector<Cell>              Cells;                   //!< Vector of cells
typedef std::vector<Cell>::iterator    C_iter;                  //!< Iterator for cell vector
typedef std::vector<JCell>             JCells;                  //!< Vector of source cells
typedef std::vector<JCell>::iterator   JC_iter;                 //!< Iterator for source cell vector

typedef std::pair<C_iter,C_iter>       Pair;                    //!< Pair of interacting cells
typedef std::stack<Pair>               Pairs;                   //!< Stack of interacting cell pairs
typedef std::list<C_iter>              List;                    //!< Interaction list
typedef std::list<C_iter>::iterator    L_iter;                  //!< Iterator for interaction list vector
typedef std::vector<List>              Lists;                   //!< Vector of interaction lists
typedef std::map<C_iter,int>           Map;                     //!< Map of interaction lists
typedef std::map<C_iter,int>::iterator M_iter;                  //!< Iterator for interation list map
typedef std::vector<Map>               Maps;                    //!< Vector of map of interaction lists

#endif
