
#include "vpic_config.h"
#include "vpic_iface.h"

#include <cassert>
#include <mrc_common.h>

// ----------------------------------------------------------------------
// C wrappers

Simulation* Simulation_create()
{
  mpi_printf(psc_comm_world, "*** Initializing\n" );
  return new Simulation();
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
  return reinterpret_cast<struct material*>(sim->define_material(name, eps, mu, sigma, zeta));
}

void Simulation_define_field_array(Simulation* sim, double damp)
{
  sim->define_field_array(damp);
}

struct species * Simulation_define_species(Simulation* sim, const char *name, double q, double m,
					   double max_local_np, double max_local_nm,
					   double sort_interval, double sort_out_of_place)
{
  return reinterpret_cast<struct species*>(sim->define_species(name, q, m, max_local_np, max_local_nm,
							       sort_interval, sort_out_of_place));
}

Particles * Simulation_get_particles(Simulation *sim)
{
  return &sim->particles_;
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

double Rng_uniform(Rng *rng, double lo, double hi)
{
  return rng->uniform(lo, hi);
}

double Rng_normal(Rng *rng, double mu, double sigma)
{
  return rng->normal(mu, sigma);
}

// ----------------------------------------------------------------------
// Simulation_set_params

void Simulation_set_params(Simulation* sim, int num_step, int status_interval,
			   int sync_shared_interval, int clean_div_e_interval,
			   int clean_div_b_interval)
{
  sim->setParams(num_step, status_interval,
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

int Simulation_mprts_get_nr_particles(Simulation*sim, Particles *vmprts)
{
  return sim->mprts_get_nr_particles(*vmprts);
}

// ----------------------------------------------------------------------
// Simulation_mprts_reserve_all
//
// This is a bit iffy, since we don't really want to reallocate stuff here,
// at least for now, and we wouldn't be able to know how to split this into
// the different species, anyway.

void Simulation_mprts_reserve_all(Simulation* sim, Particles* vmprts, int n_patches,
				 int *n_prts_by_patch)
{
  assert(n_patches == 1);

  for (int p = 0; p < n_patches; p++) {
    int n_prts = 0, n_prts_alloced = 0;
    for (auto sp = vmprts->cbegin(); sp != vmprts->cend(); ++sp) {
      n_prts += sp->np;
      n_prts_alloced += sp->max_np;
    }
#if 0
    if (n_prts_by_patch[p] != n_prts) {
      mprintf("vpic_mparticles_reserve_all: %d (currently %d max %d)\n",
	      n_prts_by_patch[p], n_prts, n_prts_alloced);
    }
#endif
    assert(n_prts_by_patch[p] <= n_prts_alloced);
  }
}

// ----------------------------------------------------------------------
// Simulation_mprts_resize_all
//
// Even more iffy, since can't really resize the per-species arrays, since we don't
// know how the total # of particles we're given should be divided up

void Simulation_mprts_resize_all(Simulation* sim, Particles* vmprts, int n_patches,
				int *n_prts_by_patch)
{
  assert(n_patches == 1);
  
  // we can't resize to the numbers given, unless it's "resize to 0", we'll just do nothing
  // The mparticles conversion function should call resize_all() itself first, resizing to
  // 0, and then using push_back, which will increase the count back to the right value

  if (n_prts_by_patch[0] == 0) {
    for (auto sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
      sp->np = 0;
    }
  } else {
#if 0
    int cur_n_prts_by_patch[n_patches];
    vpic_mparticles_get_size_all(vmprts, n_patches, cur_n_prts_by_patch);

    mprintf("vpic_mparticles_resize_all: ignoring %d -> %d\n",
	    cur_n_prts_by_patch[0], n_prts_by_patch[0]);
#endif
  }
}

// ----------------------------------------------------------------------
// Simulation_mprts_get_size_all

void Simulation_mprts_get_size_all(Simulation* sim, Particles *vmprts, int n_patches,
				   int *n_prts_by_patch)
{
  assert(n_patches == 1);
  n_prts_by_patch[0] = sim->mprts_get_nr_particles(*vmprts);
}

// ----------------------------------------------------------------------
// Simulation_mprts_get_grid_nx_dx

void Simulation_mprts_get_grid_nx_dx(Simulation* sim, Particles* vmprts, int *nx, float *dx)
{
  Grid *g = vmprts->grid();
  nx[0] = g->nx;
  nx[1] = g->ny;
  nx[2] = g->nz;
  dx[0] = g->dx;
  dx[1] = g->dy;
  dx[2] = g->dz;
}

// ----------------------------------------------------------------------
// Simulation_mprts_push_back

void Simulation_mprts_push_back(Simulation* sim, Particles* vmprts, const struct vpic_mparticles_prt *prt)
{
  for (auto sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
    if (sp->id == prt->kind) {
      assert(sp->np < sp->max_np);
      // the below is inject_particle_raw()
      Particles::Particle * RESTRICT p = sp->p + (sp->np++);
      p->dx = prt->dx[0]; p->dy = prt->dx[1]; p->dz = prt->dx[2]; p->i = prt->i;
      p->ux = prt->ux[0]; p->uy = prt->ux[1]; p->uz = prt->ux[2]; p->w = prt->w;
      return;
    }
  }
  mprintf("prt->kind %d not found in species list!\n", prt->kind);
  assert(0);
}


// ======================================================================

double Simulation_mflds_synchronize_tang_e_norm_b(Simulation* sim, FieldArray* vmflds)
{
  double err;
  TIC err = vmflds->synchronize_tang_e_norm_b(); TOC(synchronize_tang_e_norm_b, 1);
  return err;
}

void Simulation_mflds_compute_div_b_err(Simulation* sim, FieldArray* vmflds)
{
  TIC vmflds->compute_div_b_err(); TOC(compute_div_b_err, 1);
}

double Simulation_mflds_compute_rms_div_b_err(Simulation* sim, FieldArray* vmflds)
{
  double err;
  TIC err = vmflds->compute_rms_div_b_err(); TOC(compute_rms_div_b_err, 1);
  return err;
}

void Simulation_mflds_clean_div_b(Simulation* sim, FieldArray* vmflds)
{
  TIC vmflds->clean_div_b(); TOC(clean_div_b, 1);
}

void Simulation_mflds_compute_div_e_err(Simulation* sim, FieldArray* vmflds)
{
  TIC vmflds->compute_div_e_err(); TOC(compute_div_e_err, 1);
}

double Simulation_mflds_compute_rms_div_e_err(Simulation* sim, FieldArray* vmflds)
{
  double err;
  TIC err = vmflds->compute_rms_div_e_err(); TOC(compute_rms_div_e_err, 1);
  return err;
}

void Simulation_mflds_clean_div_e(Simulation* sim, FieldArray* vmflds)
{
  TIC vmflds->clean_div_e(); TOC(clean_div_e, 1);
}


void Simulation_mflds_clear_rhof(Simulation* sim, FieldArray* vmflds)
{
  TIC vmflds->clear_rhof(); TOC(clear_jf, 1);
}

void Simulation_mflds_synchronize_rho(Simulation* sim, FieldArray* vmflds)
{
  TIC vmflds->synchronize_rho(); TOC(compute_curl_b, 1);
}

void Simulation_mflds_compute_rhob(Simulation* sim, FieldArray* vmflds)
{
  TIC vmflds->compute_rhob(); TOC(compute_rhob, 1);
}

void Simulation_mflds_compute_curl_b(Simulation* sim, FieldArray* vmflds)
{
  TIC vmflds->compute_curl_b(); TOC(compute_curl_b, 1);
}

// ======================================================================

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
  sim->uncenter_p(vmprts, vmflds);
}

// ----------------------------------------------------------------------
// Simulation_moments_run

void Simulation_moments_run(Simulation* sim, HydroArray *hydro_array, Particles *vmprts, int kind)
{
  sim->moments_run(hydro_array, vmprts, kind);
}

// ----------------------------------------------------------------------
// substeps of a time integration step

void Simulation_sort_mprts(Simulation *sim, Particles *vmprts, int step)
{
  sim->sort_mprts(*vmprts, step);
}

void Simulation_collision_run(Simulation* sim)
{
  sim->collision_run();
}

void Simulation_push_mprts(Simulation *sim, Particles *vmprts, FieldArray *vmflds)
{
  sim->push_mprts(*vmprts, *vmflds);
}

void Simulation_push_mflds_H(Simulation* sim, FieldArray *vmflds, double frac)
{
  sim->push_mflds_H(*vmflds, frac);
}

void Simulation_push_mflds_E(Simulation* sim, FieldArray *vmflds, double frac)
{
  sim->push_mflds_E(*vmflds, frac);
}

void Simulation_field_injection(Simulation* sim)
{
  sim->field_injection();
}

void Simulation_push_mprts_prep(Simulation *sim, FieldArray *vmflds)
{
  sim->push_mprts_prep(*vmflds);
}

