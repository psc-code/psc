
#include "vpic.h"

#include "vpic_iface.h"

#include <mrc_common.h>

#undef sim_log
#define sim_log(x) do {                                \
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);    \
    if( rank==0 ) {				       \
      std::cerr << "SIM_LOG: " << x << std::endl;      \
      std::cerr.flush();                               \
    }                                                  \
  } while(0)


#define should_dump(x)                                                  \
  (diag->x##_interval>0 && remainder(step, diag->x##_interval) == 0)

void vpic_simulation_diagnostics(vpic_simulation *simulation, VpicDiag *diag)
{
  int64_t step = simulation->step();

  /*--------------------------------------------------------------------------
   * Restart dump
   *------------------------------------------------------------------------*/

  if(step && !(step%diag->restart_interval)) {
    if(!diag->rtoggle) {
      diag->rtoggle = 1;
      //checkpt("restart1/restart", 0);
    }
    else {
      diag->rtoggle = 0;
      //checkpt("restart2/restart", 0);
    } // if
  } // if

  // Dump particle data

  char subdir[36];
  if ( should_dump(eparticle) && step !=0
       && step > 56*(diag->fields_interval)  ) {
    // if ( should_dump(eparticle) && step !=0 ) {
    sprintf(subdir,"particle/T.%lld",step);
    simulation->dump_mkdir(subdir);
    sprintf(subdir,"particle/T.%lld/eparticle",step);
    simulation->dump_particles("electron", subdir);
  }

  if ( should_dump(Hparticle) && step !=0
       && step > 56*(diag->fields_interval)  ) {
    sprintf(subdir,"particle/T.%lld/Hparticle",step);
    simulation->dump_particles("ion", subdir);
  }

}

