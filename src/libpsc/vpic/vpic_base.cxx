
#include "vpic_iface.h"
#include "vpic_mfields.h"
#include "vpic_mparticles.h"
#include "vpic_push_particles.h"

#include <vpic.h>

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


void vpic_simulation_init(struct vpic_simulation_info *info)
{
  if( world_rank==0 ) log_printf( "*** Initializing\n" );
  simulation = new vpic_simulation;

  // Call the user initialize the simulation

  int argc = 0; char **argv = NULL;
  TIC simulation->user_initialization(argc, argv); TOC( user_initialization, 1 );

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

void vpic_simulation_init2(vpic_push_particles *vpushp, vpic_mfields *vmflds,
			   vpic_mparticles *vmprts)
{
  double err;

  // Load fields not initialized by the user

  if( simulation->rank()==0 ) MESSAGE(( "Initializing radiation damping fields" ));
  vmflds->compute_curl_b();

  if( simulation->rank()==0 ) MESSAGE(( "Initializing bound charge density" ));
  vmflds->clear_rhof();
  vmflds->accumulate_rho_p(vmprts);
  vmflds->synchronize_rho();
  vmflds->compute_rhob();

  // Internal sanity checks

  if( simulation->rank()==0 ) MESSAGE(( "Checking electric field divergence" ));

  vmflds->compute_div_e_err();
  err = vmflds->compute_rms_div_e_err();
  if( simulation->rank()==0 ) MESSAGE(( "RMS error = %e (charge/volume)", err ));
  vmflds->clean_div_e();

  if( simulation->rank()==0 ) MESSAGE(( "Rechecking interdomain synchronization" ));
  err = vmflds->synchronize_tang_e_norm_b();
  if( simulation->rank()==0 ) MESSAGE(( "Error = %e (arb units)", err ));
    
  if( vmprts->species_list ) {
    if( simulation->rank()==0 ) MESSAGE(( "Uncentering particles" ));
    vpushp->load_interpolator_array(vmflds);
  }
  vpushp->uncenter_p(vmprts);
}

void vpic_inc_step(int step)
{
  simulation->grid->step++;
  assert(simulation->grid->step == step);
}


