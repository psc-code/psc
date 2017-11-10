
#include "vpic_iface.h"

#include <vpic.h>

#include <mpi.h>

#define mprintf(fmt...) do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); { printf("[%d] ", __rank); printf(fmt); } } while(0)

#define MHERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)

// ----------------------------------------------------------------------
// vpic_base_init

extern vpic_simulation *simulation;

void
vpic_base_init(struct vpic_info *info)
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
}

void vpic_step(void)
{
  vpic_performance_sort();
  vpic_clear_accumulator_array();
  vpic_collisions();
  vpic_advance_p();
  vpic_emitter();
  vpic_reduce_accumulator_array();
  vpic_boundary_p();
  vpic_calc_jf();
  vpic_current_injection();

  // Half advance the magnetic field from B_0 to B_{1/2}
  vpic_advance_b(0.5);
  // Advance the electric field from E_0 to E_1
  vpic_advance_e(1.0);
  vpic_field_injection();
  // Half advance the magnetic field from B_{1/2} to B_1
  vpic_advance_b(0.5);

  // Divergence clean e
  if( (simulation->clean_div_e_interval>0) && ((simulation->grid->step % simulation->clean_div_e_interval)==0) ) {
    vpic_clean_div_e();
  }

  // Divergence clean b
  if( (simulation->clean_div_b_interval>0) && ((simulation->grid->step % simulation->clean_div_b_interval)==0) ) {
    vpic_clean_div_b();
  }

  // Synchronize the shared faces
  if( (simulation->sync_shared_interval>0) && ((simulation->grid->step % simulation->sync_shared_interval)==0) ) {
    vpic_sync_faces();
  }

  // Fields are updated ... load the interpolator for next time step and
  // particle diagnostics in user_diagnostics if there are any particle
  // species to worry about
  vpic_load_interpolator_array();

  simulation->grid->step++;

  vpic_print_status();
  vpic_diagnostics();
}

void
vpic_base_integrate()
{
  if( world_rank==0 ) log_printf( "*** Advancing\n" );
  double elapsed = wallclock();
  while (!vpic_done()) {
    vpic_step();
  }
  elapsed = wallclock() - elapsed;
  if (world_rank==0) {
    int  s = (int)elapsed, m  = s/60, h  = m/60, d  = h/24, w = d/ 7;
    /**/ s -= m*60,        m -= h*60, h -= d*24, d -= w*7;
    log_printf( "*** Done (%gs / %iw:%id:%ih:%im:%is elapsed)\n",
                elapsed, w, d, h, m, s );
  }
}
