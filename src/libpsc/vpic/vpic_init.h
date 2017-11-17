
#ifndef VPIC_INIT_H
#define VPIC_INIT_H

#include "vpic_iface.h"
#include "vpic.h"

// ----------------------------------------------------------------------
// globals_diag

struct globals_diag {
  int interval;
  int energies_interval;
  int fields_interval;
  int ehydro_interval;
  int Hhydro_interval;
  int eparticle_interval;
  int Hparticle_interval;
  int restart_interval;

  // state
  int rtoggle;               // enables save of last 2 restart dumps for safety
  // Output variables
  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hHdParams;
  std::vector<DumpParameters *> outputParams;
};

// ----------------------------------------------------------------------

void vpic_simulation_diagnostics(vpic_simulation *simulation, globals_diag *diag);
void vpic_simulation_setup_diagnostics(vpic_simulation *simulation, globals_diag *diag,
				       species_t *electron, species_t *ion);


#endif
