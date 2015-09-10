#include "args.h"
#include "bound_box.h"
#include "dataset.h"
#include "logger.h"
#include "verify.h"

#include "build_tree.h"
#include "kernel.h"
#include "traversal.h"

int main(int argc, char ** argv) {
  Args args(argc,argv);
  Bodies bodies, bodies2, jbodies;
  BoundBox boundBox(args.nspawn);
  Bounds bounds;
  Dataset data;
  Verify verify;
  const int numBodies=args.numBodies;
  kernel::wavek = complex_t(10.,1.) / (2 * M_PI);
  logger::verbose = args.verbose;
  bodies = data.initBodies(args.numBodies, args.distribution, 0);
  bodies.resize(numBodies);
  logger::startTimer("Total FMM");
  bounds = boundBox.getBounds(bodies);
  Bodies buffer(numBodies);
  Cells cells = buildTree(bodies, buffer, bounds);
  evaluate(cells);
  logger::stopTimer("Total FMM");
  jbodies = bodies;
  const int numTargets = 100;
  data.sampleBodies(bodies, numTargets);
  bodies2 = bodies;
  data.initTarget(bodies);
  cells.resize(2);
  C_iter Ci = cells.begin();
  C_iter Cj = cells.begin() + 1;
  Ci->BODY = bodies.begin();
  Ci->NBODY = bodies.size();
  Cj->BODY = jbodies.begin();
  Cj->NBODY = jbodies.size();
  logger::startTimer("Total Direct");
  real_t eps2 = 0;
  vec3 Xperiodic = 0;
  bool mutual = false;
  kernel::P2P(Ci, Cj, eps2, Xperiodic, mutual);
  logger::stopTimer("Total Direct");
  std::complex<double> potDif = verify.getDifScalar(bodies, bodies2);
  std::complex<double> potNrm = verify.getNrmScalar(bodies);
  std::complex<double> accDif = verify.getDifVector(bodies, bodies2);
  std::complex<double> accNrm = verify.getNrmVector(bodies);
  logger::printTitle("FMM vs. direct");
  verify.print("Rel. L2 Error (pot)",std::sqrt(potDif/potNrm));
  verify.print("Rel. L2 Error (acc)",std::sqrt(accDif/accNrm));
}
