
#include "vpic_iface.h"
#include "vpic_diag.h"

#include <cassert>

// ----------------------------------------------------------------------
// C wrappers

Simulation* Simulation_create()
{
  extern vpic_simulation *simulation;
  assert(!simulation);

  if( world_rank==0 ) log_printf( "*** Initializing\n" );
  simulation = new vpic_simulation;
  return new Simulation(simulation);
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
  sim->pDiag_ = new globals_diag(interval);
}

void Simulation_diagnostics_setup(Simulation* sim)
{
  sim->pDiag_->setup(sim);
}

void Simulation_diagnostics_run(Simulation* sim, struct psc_harris *sub)
{
  sim->pDiag_->run(sim);
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

void Simulation_set_params(Simulation* sim, int num_step, int status_interval,
			   int sync_shared_interval, int clean_div_e_interval,
			   int clean_div_b_interval)
{
  vpic_simulation *vpic = sim->simulation_;
  vpic->num_step             = num_step;
  vpic->status_interval      = status_interval;
  vpic->sync_shared_interval = sync_shared_interval;
  vpic->clean_div_e_interval = clean_div_e_interval;
  vpic->clean_div_b_interval = clean_div_b_interval;
}

// ----------------------------------------------------------------------
// Simulation_inject_particle

void Simulation_inject_particle(Simulation* sim, Particles *vmprts, int p,
				const struct psc_particle_inject *prt)
{
  assert(p == 0);
  species_t *sp = find_species_id(prt->kind, vmprts->sl_);

  sim->simulation_->inject_particle(sp, prt->x[0], prt->x[1], prt->x[2],
				    prt->u[0], prt->u[1], prt->u[2], prt->w, 0., 0);
}

// ----------------------------------------------------------------------
// Simulation_collision_run

void Simulation_collision_run(Simulation* sim)
{
  // Note: Particles should not have moved since the last performance sort
  // when calling collision operators.
  // FIXME: Technically, this placement of the collision operators only
  // yields a first order accurate Trotter factorization (not a second
  // order accurate factorization).

  if (sim->simulation_->collision_op_list) {
    // FIXME: originally, vpic_clear_accumulator_array() was called before this.
    // It's now called later, though. I'm not sure why that would be necessary here,
    // but it needs to be checked.
    // The assert() below doesn't unfortunately catch all cases where this might go wrong
    // (ie., it's missing the user_particle_collisions())

    assert(0);
    TIC apply_collision_op_list(sim->simulation_->collision_op_list); TOC(collision_model, 1);
  }
  TIC sim->simulation_->user_particle_collisions(); TOC(user_particle_collisions, 1);
}

// ----------------------------------------------------------------------
// Simulation_emitter

void Simulation_emitter(Simulation* sim)
{
  if (sim->simulation_->emitter_list)
    TIC apply_emitter_list(sim->simulation_->emitter_list); TOC(emission_model, 1);
  TIC sim->simulation_->user_particle_injection(); TOC(user_particle_injection, 1);
}

// ----------------------------------------------------------------------
// Simulation_current_injection

void Simulation_current_injection(Simulation* sim)
{
  TIC sim->simulation_->user_current_injection(); TOC(user_current_injection, 1);
}

// ----------------------------------------------------------------------
// Simulation_field_injection

void Simulation_field_injection(Simulation* sim)
{
  // Let the user add their own contributions to the electric field. It is the
  // users responsibility to insure injected electric fields are consistent
  // across domains.

  TIC sim->simulation_->user_field_injection(); TOC(user_field_injection, 1);
}

// ----------------------------------------------------------------------
// Simulation_moments_run

void Simulation_moments_run(Simulation* sim, HydroArray *hydro_array, Particles *vmprts, int kind)
{
  // This relies on load_interpolator_array() having been called earlier
  clear_hydro_array(hydro_array);
  species_t *sp;
  LIST_FOR_EACH(sp, vmprts->sl_) {
    if (sp->id == kind) {
      accumulate_hydro_p(hydro_array, sp, sim->interpolator_array_);
      break;
    }
  }
  
  synchronize_hydro_array(hydro_array);
}

void Simulation_advance_b(Simulation* sim, FieldArray *fa, double frac)
{
  TIC sim->advance_b(*fa, frac); TOC(advance_b, 1);
}

