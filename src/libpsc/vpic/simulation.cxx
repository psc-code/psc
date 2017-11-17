
#include "simulation.h"

extern vpic_simulation *simulation;

struct user_global_t {
  struct globals_diag diag;
};

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
  user_global_t *user_global = (struct user_global_t *) simulation->user_global;
  sim->pDiag_ = new(&user_global->diag) globals_diag(interval);
}

void Simulation_diagnostics_setup(struct Simulation *sim)
{
  species_t *electron = simulation->find_species("electron");
  species_t *ion = simulation->find_species("ion");
  vpic_simulation_setup_diagnostics(simulation, sim->pDiag_, electron, ion);
}

void Simulation_diagnostics_run(struct Simulation *sim, struct psc_harris *sub)
{
  TIC vpic_simulation_diagnostics(simulation, sim->pDiag_); TOC( user_diagnostics, 1 );
}




