
#include "vpic_iface.h"
#include "vpic_mfields.h"
#include "vpic_mparticles.h"
#include "vpic_push_particles.h"
#include "vpic_marder.h"

#include <vpic.h>

#include <cassert>

extern vpic_simulation *simulation;

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


