#include "dataset.h"
#include "serialfmm.h"

int main() {
  const int numBodies = 10000;
  const int numTarget = 100;
  Bodies bodies(numBodies);
  Cells cells;
  SerialFMM FMM;
  FMM.setKernel("BiotSavart");
  FMM.initialize();
  FMM.IMAGES = 1;
  FMM.IMAGEDIM[0] = 1;
  FMM.IMAGEDIM[1] = 0;
  FMM.IMAGEDIM[2] = 0;
  FMM.THETA = 1 / sqrtf(4);
  bool printNow = true;

  FMM.startTimer("Set bodies   ");
  Dataset Data;
  Data.setKernel("BiotSavart");
  Data.random(bodies);
  FMM.stopTimer("Set bodies   ",printNow);

  FMM.startTimer("Set domain   ");
  FMM.setDomain(bodies);
  FMM.stopTimer("Set domain   ",printNow);

  FMM.startTimer("Build tree   ");
  FMM.bottomup(bodies,cells);
  FMM.stopTimer("Build tree   ",printNow);

  Cells jcells = cells;
  FMM.startTimer("Traversal    ");
  FMM.downward(cells,jcells,1);
  FMM.stopTimer("Traversal    ",printNow);

  FMM.startTimer("Direct sum   ");
  Bodies jbodies = bodies;
  Data.sampleBodies(bodies,numTarget);
  Bodies bodies2 = bodies;
  Data.initTarget(bodies2);
  FMM.evalP2P(bodies2,jbodies);
  FMM.stopTimer("Direct sum   ",printNow);

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  Data.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  if(printNow) Data.printError(diff1,norm1,diff2,norm2);
  FMM.finalize();
}
