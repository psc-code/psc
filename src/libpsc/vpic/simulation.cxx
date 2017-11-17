
#include "simulation.h"

extern vpic_simulation *simulation;

// ----------------------------------------------------------------------
// C wrappers

Simulation *Simulation_create()
{
  return new Simulation;
}

void Simulation_delete(Simulation *sim)
{
  delete sim;
}

void Simulation_diagnostics_init(struct Simulation *sim, int interval)
{
  sim->pDiag_ = new globals_diag(interval);
}

void Simulation_diagnostics_setup(struct Simulation *sim)
{
  sim->pDiag_->setup();
}

void Simulation_diagnostics_run(struct Simulation *sim, struct psc_harris *sub)
{
  sim->pDiag_->run();
}




