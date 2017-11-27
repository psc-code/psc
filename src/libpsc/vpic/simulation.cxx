
#include "vpic_iface.h"

#include <cassert>

// ----------------------------------------------------------------------
// C wrappers

Simulation* Simulation_create()
{
  if( world_rank==0 ) log_printf( "*** Initializing\n" );
  return new Simulation(new SimulationBase);
}

void Simulation_delete(Simulation* sim)
{
  delete sim;
}

void Simulation_setup_grid(Simulation* sim, double dx[3], double dt,
			   double cvac, double eps0)
{
  sim->setup_grid(dx, dt, cvac, eps0);
}

void Simulation_define_periodic_grid(Simulation* sim, double xl[3],
				     double xh[3], int gdims[3], int np[3])
{
  sim->define_periodic_grid(xl, xh, gdims, np);
}

void Simulation_set_domain_field_bc(Simulation* sim, int boundary, int bc)
{
  sim->set_domain_field_bc(boundary, bc);
}

void Simulation_set_domain_particle_bc(Simulation* sim, int boundary, int bc)
{
  sim->set_domain_particle_bc(boundary, bc);
}

struct material *Simulation_define_material(Simulation* sim, const char *name,
					    double eps, double mu,
					    double sigma, double zeta)
{
  return sim->define_material(name, eps, mu, sigma, zeta);
}

void Simulation_define_field_array(Simulation* sim, double damp)
{
  sim->define_field_array(damp);
}

struct species * Simulation_define_species(Simulation* sim, const char *name, double q, double m,
					   double max_local_np, double max_local_nm,
					   double sort_interval, double sort_out_of_place)
{
  return sim->define_species(name, q, m, max_local_np, max_local_nm,
			     sort_interval, sort_out_of_place);
}

// ----------------------------------------------------------------------
// diagnostics

void Simulation_diagnostics_init(Simulation* sim, int interval)
{
  sim->newDiag(interval);
}

void Simulation_diagnostics_setup(Simulation* sim)
{
  sim->setupDiag();
}

void Simulation_diagnostics_run(Simulation* sim)
{
  sim->runDiag();
}

// ----------------------------------------------------------------------
// Rng

void Simulation_rngPool_seed(Simulation* sim, int base)
{
  sim->rng_pool.seed(base, 0);
}

Rng *Simulation_rngPool_get(Simulation* sim, int n)
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

// ----------------------------------------------------------------------
// Simulation_set_params

void Simulation_set_params(Simulation* sim, int num_step, int status_interval,
			   int sync_shared_interval, int clean_div_e_interval,
			   int clean_div_b_interval)
{
  sim->set_params(num_step, status_interval,
		  sync_shared_interval, clean_div_e_interval, clean_div_b_interval);
}

// ----------------------------------------------------------------------
// Simulation_inject_particle

void Simulation_inject_particle(Simulation* sim, Particles *vmprts, int p,
				const struct psc_particle_inject *prt)
{
  assert(p == 0);
  static_cast<ParticlesOps*>(sim)->inject_particle(*vmprts, *sim->accumulator_, *sim->field_array_, prt);
}

// ----------------------------------------------------------------------
// Simulation_accumulate_rho_p

void Simulation_accumulate_rho_p(Simulation *sim, Particles *vmprts, FieldArray *vmflds)
{
  return sim->accumulate_rho_p(*vmprts, *vmflds);
}

// ----------------------------------------------------------------------
// Simulation_initialize

void Simulation_initialize(Simulation *sim, Particles *vmprts, FieldArray *vmflds)
{
  MPI_Comm comm = MPI_COMM_WORLD; // FIXME
  
  mpi_printf(comm, "Uncentering particles\n");
  if (!vmprts->empty()) {
    sim->interpolator_->load(*vmflds);

    for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
      TIC ::uncenter_p(&*sp, sim->interpolator_); TOC(uncenter_p, 1);
    }
  }
}

// ----------------------------------------------------------------------
// Simulation_collision_run

void Simulation_collision_run(Simulation* sim)
{
  sim->collision_run();
}

// ----------------------------------------------------------------------
// Simulation_field_injection

void Simulation_field_injection(Simulation* sim)
{
  sim->field_injection();
}

// ----------------------------------------------------------------------
// Simulation_moments_run

void Simulation_moments_run(Simulation* sim, HydroArray *hydro_array, Particles *vmprts, int kind)
{
  sim->moments_run(hydro_array, vmprts, kind);
}

void Simulation_advance_b(Simulation* sim, FieldArray *vmflds, double frac)
{
  TIC vmflds->advance_b(frac); TOC(advance_b, 1);
}

void Simulation_advance_e(Simulation* sim, FieldArray *vmflds, double frac)
{
  TIC vmflds->advance_e(frac); TOC(advance_e, 1);
}

