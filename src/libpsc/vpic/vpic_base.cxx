
#include "vpic_iface.h"

#include "vpic_init.h"
#include "vpic_mfields.h"
#include "vpic_mparticles.h"
#include "vpic_push_particles.h"

#include <vpic.h>

#include <mrc_common.h>

#include <cassert>

extern vpic_simulation *simulation;

// ======================================================================

// ----------------------------------------------------------------------
// vpic_base_init

void vpic_base_init(int *pargc, char ***pargv)
{
  static bool vpic_base_inited = false;

  if (vpic_base_inited) {
    return;
  }
  vpic_base_inited = true;
  
  //  boot_services( &argc, &argv );
  {
    // Start up the checkpointing service.  This should be first.
    
    boot_checkpt( pargc, pargv );
    
    serial.boot( pargc, pargv );
    thread.boot( pargc, pargv );
    
    // Boot up the communications layer
    // See note above about thread-core-affinity
    
    boot_mp( pargc, pargv );
    
    mp_barrier();
    _boot_timestamp = 0;
    _boot_timestamp = uptime();
  }
}

static void
vpic_simulation_get_info(struct vpic_simulation_info *info)
{
  info->num_step = simulation->num_step;

  // grid
  info->dt = simulation->grid->dt;
  info->nx[0] = simulation->grid->nx;
  info->nx[1] = simulation->grid->ny;
  info->nx[2] = simulation->grid->nz;
  info->dx[0] = simulation->grid->dx;
  info->dx[1] = simulation->grid->dy;
  info->dx[2] = simulation->grid->dz;
  info->x0[0] = simulation->grid->x0;
  info->x0[1] = simulation->grid->y0;
  info->x0[2] = simulation->grid->z0;
  info->x1[0] = simulation->grid->x1;
  info->x1[1] = simulation->grid->y1;
  info->x1[2] = simulation->grid->z1;

  // species
  info->n_kinds = num_species(simulation->species_list);
  info->kinds = new vpic_kind_info[info->n_kinds];
  species_t *sp;
  LIST_FOR_EACH( sp, simulation->species_list ) {
    info->kinds[sp->id].q = sp->q;
    info->kinds[sp->id].m = sp->m;
    info->kinds[sp->id].name = sp->name;
  }
  
  // Marder cleaning etc
  info->clean_div_e_interval = simulation->clean_div_e_interval;
  info->clean_div_b_interval = simulation->clean_div_b_interval;
  info->sync_shared_interval = simulation->sync_shared_interval;
  info->num_div_e_round = simulation->num_div_e_round;
  info->num_div_b_round = simulation->num_div_b_round;

  info->status_interval = simulation->status_interval;
}

void vpic_simulation_new()
{
  assert(!simulation);
  if( world_rank==0 ) log_printf( "*** Initializing\n" );
  simulation = new vpic_simulation;
}

void vpic_simulation_init(vpic_simulation_info *info)
{
  // Call the user initialize the simulation
  TIC simulation->user_initialization(0, 0); TOC( user_initialization, 1 );

  vpic_simulation_get_info(info);
}

// ======================================================================
// vpic_diagnostics

void vpic_diagnostics()
{
  // Let the user compute diagnostics
  TIC simulation->user_diagnostics(); TOC( user_diagnostics, 1 );
}

// ======================================================================

void vpic_inc_step(int step)
{
  simulation->grid->step++;
  assert(simulation->grid->step == step);
}


// ======================================================================
// ======================================================================

// ----------------------------------------------------------------------
// vpic_simulation_set_params

void vpic_simulation_set_params(int num_step,
				int status_interval,
				int sync_shared_interval,
				int clean_div_e_interval,
				int clean_div_b_interval)
{
  simulation->num_step             = num_step;
  simulation->status_interval      = status_interval;
  simulation->sync_shared_interval = sync_shared_interval;
  simulation->clean_div_e_interval = clean_div_e_interval;
  simulation->clean_div_b_interval = clean_div_b_interval;
}

// ----------------------------------------------------------------------
// vpic_simulation_setup_grid

void vpic_simulation_setup_grid(double dx[3], double dt, double cvac, double eps0)
{
  grid_t *grid = simulation->grid;

  grid->dx = dx[0];
  grid->dy = dx[1];
  grid->dz = dx[2];
  grid->dt = dt;
  grid->cvac = cvac;
  grid->eps0 = eps0;
}

// ----------------------------------------------------------------------
// vpic_simulation_define_periodic_grid

void vpic_simulation_define_periodic_grid(double xl[3], double xh[3],
					  int gdims[3], int np[3])
{
  simulation->define_periodic_grid(xl[0], xl[1], xl[2], xh[0], xh[1], xh[2],
				   gdims[0], gdims[1], gdims[2], np[0], np[1], np[2]);
}

// ----------------------------------------------------------------------
// vpic_simulation_set_domain_field_bc

// FIXME, replicated
///Possible boundary conditions for fields
enum {
  BND_FLD_OPEN,
  BND_FLD_PERIODIC,
  BND_FLD_UPML,
  BND_FLD_TIME,
  BND_FLD_CONDUCTING_WALL,
  BND_FLD_ABSORBING,
};

///Possible boundary conditions for particles
enum {
  BND_PART_REFLECTING,
  BND_PART_PERIODIC,
  BND_PART_ABSORBING,
  BND_PART_OPEN,
};

void
vpic_simulation_set_domain_field_bc(int boundary, int bc)
{
  int fbc;
  switch (bc) {
  case BND_FLD_CONDUCTING_WALL: fbc = pec_fields   ; break;
  case BND_FLD_ABSORBING:       fbc = absorb_fields; break;
  default: assert(0);
  }
  simulation->set_domain_field_bc(boundary, fbc);
}

void
vpic_simulation_set_domain_particle_bc(int boundary, int bc)
{
  int pbc;
  switch (bc) {
  case BND_PART_REFLECTING: pbc = reflect_particles; break;
  case BND_PART_ABSORBING:  pbc = absorb_particles ; break;
  default: assert(0);
  }
  simulation->set_domain_particle_bc(boundary, pbc);
}

struct material *vpic_simulation_define_material(const char *name, double eps, double mu,
				     double sigma, double zeta)
{
  return simulation->define_material(name, eps, mu, sigma, zeta);
}

void vpic_simulation_define_field_array(field_array_t *fa, double damp)
{
  simulation->define_field_array(fa, damp);
}

void vpic_simulation_set_region_resistive_harris(vpic_harris_params *prm,
						 globals_physics *phys,
						 double dx[3],
						 double thickness,
						 struct material *resistive)
{
  // Define resistive layer surrounding boundary --> set thickness=0
  // to eliminate this feature
#define resistive_layer ((prm->open_bc_x && x < dx[0]*thickness) ||	\
			 (prm->open_bc_x && x > phys->Lx-dx[0]*thickness)	\
                         || z <-phys->Lz/2+dx[2]*thickness  || z > phys->Lz/2-dx[2]*thickness )

  if (thickness > 0) {
    log_printf("Setting resistive layer of thickness %g", thickness);
    // FIXME!!!
    assert(0);
#if 0
#define field simulation->field
    grid_t *grid = simulation->grid;
    set_region_material(resistive_layer, resistive, resistive);
#undef field
#endif
  }
}

struct species *
vpic_simulation_define_species(const char *name, double q, double m,
			       double max_local_np, double max_local_nm,
			       double sort_interval, double sort_out_of_place)
{
  return simulation->define_species(name, q, m, max_local_np, max_local_nm,
				    sort_interval, sort_out_of_place);
}

void vpic_simulation_inject_particle(struct species * sp,
				     double x,  double y,  double z,
				     double ux, double uy, double uz,
				     double w,  double age, bool update_rhob)
{
  simulation->inject_particle(sp, x, y, z, ux, uy, uz, w, age, update_rhob);
}

struct species *vpic_simulation_find_species(const char *name)
{
  return simulation->find_species(name);
}


