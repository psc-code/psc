
#include "vpic_iface.h"
#include "vpic_diag.h"

#include "simulation.h"

#include "mrc_common.h"

// ================================================================================
// globals_diag implementation

void vpic_simulation_diagnostics(vpic_simulation *simulation, VpicDiag *diag);
void vpic_simulation_setup_diagnostics(vpic_simulation *simulation, VpicDiag *diag);


VpicDiag::VpicDiag(vpic_simulation *simulation, int interval_)
  : simulation_(simulation)
{
  rtoggle = 0;
  
  interval = interval_;
  fields_interval = interval_;
  ehydro_interval = interval_;
  Hhydro_interval = interval_;
  eparticle_interval = 8 * interval_;
  Hparticle_interval = 8 * interval_;
  restart_interval = 8000;

  energies_interval = 50;

  MPI_Comm comm = MPI_COMM_WORLD;
  mpi_printf(comm, "interval = %d\n", interval);
  mpi_printf(comm, "energies_interval: %d\n", energies_interval);
}

void VpicDiag::setup()
{
  vpic_simulation_setup_diagnostics(simulation_, this);
}

void VpicDiag::run()
{
  TIC vpic_simulation_diagnostics(simulation_, this); TOC(user_diagnostics, 1);
}

