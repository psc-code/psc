
#ifndef SIMULATION_H
#define SIMULATION_H

#include "vpic_iface.h"

#include "vpic_init.h" // FIXME, bad name for _diag

struct Simulation {
  ~Simulation();
  
  globals_diag *pDiag_;
};

// ======================================================================

Simulation::~Simulation()
{
  delete pDiag_;
}

#endif



