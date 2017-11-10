
#include "vpic_iface.h"
#include "vpic_mfields.h"
#include "vpic_mparticles.h"
#include "vpic_push_particles.h"
#include "vpic_marder.h"

#include <vpic.h>

#include <cassert>

extern vpic_simulation *simulation;

// ======================================================================
// vpic_marder

// ----------------------------------------------------------------------
// vpic_marder_create

struct vpic_marder *
vpic_marder_create()
{
  return new vpic_marder;
}

// ----------------------------------------------------------------------
// vpic_marder_ctor_from_simulation

void
vpic_marder_ctor_from_simulation(struct vpic_marder *vmarder)
{
  vmarder->clean_div_e_interval = simulation->clean_div_e_interval;
  vmarder->clean_div_b_interval = simulation->clean_div_b_interval;
  vmarder->sync_shared_interval = simulation->sync_shared_interval;
  vmarder->num_div_e_round = simulation->num_div_e_round;
  vmarder->num_div_b_round = simulation->num_div_b_round;
}

// ----------------------------------------------------------------------
// vpic_marder_run

void
vpic_marder_run(struct vpic_marder *vmarder, struct vpic_mfields *vmflds,
		struct vpic_mparticles *vmprts, int step)
{
  // Divergence clean e
  int clean_div_e_interval = vmarder->clean_div_e_interval;
  if (clean_div_e_interval > 0 && step % clean_div_e_interval == 0) {
    vpic_marder_clean_div_e(vmarder, vmflds, vmprts);
  }
  
  // Divergence clean b
  int clean_div_b_interval = vmarder->clean_div_b_interval;
  if (clean_div_b_interval > 0 && step % clean_div_b_interval == 0) {
    vpic_marder_clean_div_b(vmarder, vmflds);
  }

  // Synchronize the shared faces
  int sync_shared_interval = vmarder->sync_shared_interval;
  if (sync_shared_interval > 0 && step % sync_shared_interval == 0) {
    vpic_marder_sync_faces(vmarder, vmflds);
  }
}

// ======================================================================

// ----------------------------------------------------------------------
// vpic_base_init

void
vpic_base_init(struct vpic_simulation_info *info)
{
  if (simulation) {
    return;
  }
  
  //  boot_services( &argc, &argv );
  {
    int argc = 0;
    char *_argv[] = {}, **argv = _argv;
    int *pargc = &argc;
    char ***pargv = &argv;
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

  if( world_rank==0 ) log_printf( "*** Initializing\n" );
  simulation = new vpic_simulation;
  simulation->initialize( 0, NULL );

  info->num_step = simulation->num_step;
  info->clean_div_e_interval = simulation->clean_div_e_interval;
  info->clean_div_b_interval = simulation->clean_div_b_interval;
  info->sync_shared_interval = simulation->sync_shared_interval;
  info->status_interval = simulation->status_interval;
}

void vpic_inc_step(int step)
{
  simulation->grid->step++;
  assert(simulation->grid->step == step);
}


