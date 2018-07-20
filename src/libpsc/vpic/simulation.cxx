
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

