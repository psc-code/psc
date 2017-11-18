
#include "simulation.h"

#include <cassert>

extern vpic_simulation *simulation;

// ----------------------------------------------------------------------
// C wrappers

Simulation *Simulation_create()
{
  return new Simulation(simulation);
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
  sim->set_domain_field_bc(boundary, bc);
}

void Simulation_set_domain_particle_bc(struct Simulation *sim, int boundary, int bc)
{
  sim->set_domain_particle_bc(boundary, bc);
}

struct material *Simulation_define_material(struct Simulation *sim, const char *name,
					    double eps, double mu,
					    double sigma, double zeta)
{
  return sim->define_material(name, eps, mu, sigma, zeta);
}

void Simulation_define_field_array(struct Simulation *sim, double damp)
{
  sim->define_field_array(damp);
}

struct species * Simulation_define_species(struct Simulation *sim, const char *name, double q, double m,
					   double max_local_np, double max_local_nm,
					   double sort_interval, double sort_out_of_place)
{
  return sim->define_species(name, q, m, max_local_np, max_local_nm,
			     sort_interval, sort_out_of_place);
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

void Simulation_set_params(Simulation *sim, int num_step, int status_interval,
			   int sync_shared_interval, int clean_div_e_interval,
			   int clean_div_b_interval)
{
  vpic_simulation_t *vpic = sim->simulation_;
  vpic->num_step             = num_step;
  vpic->status_interval      = status_interval;
  vpic->sync_shared_interval = sync_shared_interval;
  vpic->clean_div_e_interval = clean_div_e_interval;
  vpic->clean_div_b_interval = clean_div_b_interval;
}

