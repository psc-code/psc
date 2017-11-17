
#ifndef SIMULATION_H
#define SIMULATION_H

#include "vpic_iface.h"

#include "vpic_init.h" // FIXME, bad name for _diag
#include "util/rng/rng.h"

// ----------------------------------------------------------------------
// class Simulation

struct Simulation {
  ~Simulation();
  
  globals_diag *pDiag_;
};

// ----------------------------------------------------------------------
// class Rng

#define IN_rng
#include "util/rng/rng_private.h"

struct Rng : rng {
};

// ======================================================================
// Simulation implementation

inline Simulation::~Simulation()
{
  delete pDiag_;
}

#endif



