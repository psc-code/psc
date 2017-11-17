
#include "simulation.h"

#include <cassert>

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

void Simulation_setup_grid(Simulation *sim, double dx[3], double dt,
			   double cvac, double eps0)
{
  sim->setup_grid(dx, dt, cvac, eps0);
}

void Simulation_define_periodic_grid(struct Simulation *sim, double xl[3],
				     double xh[3], int gdims[3], int np[3])
{
  sim->define_periodic_grid(xl, xh, gdims, np);
}

void Simulation_set_domain_field_bc(struct Simulation *sim, int boundary, int bc)
{
  int fbc;
  switch (bc) {
  case BND_FLD_CONDUCTING_WALL: fbc = pec_fields   ; break;
  case BND_FLD_ABSORBING:       fbc = absorb_fields; break;
  default: assert(0);
  }
  simulation->set_domain_field_bc(boundary, fbc);
}

void Simulation_set_domain_particle_bc(struct Simulation *sim, int boundary, int bc)
{
  int pbc;
  switch (bc) {
  case BND_PART_REFLECTING: pbc = reflect_particles; break;
  case BND_PART_ABSORBING:  pbc = absorb_particles ; break;
  default: assert(0);
  }
  simulation->set_domain_particle_bc(boundary, pbc);
}

struct material *Simulation_define_material(struct Simulation *sim, const char *name,
					    double eps, double mu,
					    double sigma, double zeta)
{
  return sim->define_material(name, eps, mu, sigma, zeta);
}

void Simulation_define_field_array(struct Simulation *sim, struct field_array *fa, double damp)
{
  sim->define_field_array(fa, damp);
}

// ----------------------------------------------------------------------
// diagnostics

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

// ----------------------------------------------------------------------
// Rng

void Simulation_rngPool_seed(struct Simulation *sim, int base)
{
  sim->rng_pool.seed(base, 0);
}

Rng *Simulation_rngPool_get(struct Simulation *sim, int n)
{
  return sim->rng_pool[n];
}

double Rng_uniform(struct Rng *rng, double lo, double hi)
{
  return rng->uniform(lo, hi);
}

double Rng_normal(struct Rng *rng, double mu, double sigma)
{
  return rng->normal(mu, sigma);
}
