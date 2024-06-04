#ifndef serialfmm_h
#define serialfmm_h
#include "bottomup.h"

//! Serial FMM interface
class SerialFMM : public BottomUp {
public:
//! Constructor
  SerialFMM() : BottomUp() {
    preCalculation();
  }
//! Destructor
  ~SerialFMM() {
    postCalculation();
  }

//! Bottomup tree constructor interface. Input: bodies, Output: cells
  void bottomup(Bodies &bodies, Cells &cells) {
    BottomUp::setIndex(bodies);                                 // Set index of cells

    buffer.resize(2*bodies.size());                             // Resize sort buffer
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order

    Cells twigs;                                                // Twigs are cells at the bottom of tree
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs

    Cells sticks;                                               // Sticks are twigs from other processes not twigs here
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
  }
};

#endif
