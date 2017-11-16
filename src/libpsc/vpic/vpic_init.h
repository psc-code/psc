
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
  double quota_sec;          // Run quota in seconds

  // state
  int rtoggle;               // enables save of last 2 restart dumps for safety
  // Output variables
  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hHdParams;
  std::vector<DumpParameters *> outputParams;

  // Vadim: modified restart machinary
  int write_end_restart; // global flag for all to write restart files
};

// ----------------------------------------------------------------------

void user_init(vpic_simulation *simulation, vpic_params *prm, struct psc_harris *harris,
	       globals_diag *diag);

void vpic_simulation_diagnostics(vpic_simulation *simulation, vpic_params *prm,
				 globals_diag *diag);

#endif
