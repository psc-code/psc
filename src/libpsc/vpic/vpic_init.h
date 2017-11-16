
#ifndef VPIC_INIT_H
#define VPIC_INIT_H

#include "vpic_iface.h"
#include "vpic.h"

// ----------------------------------------------------------------------
// params

struct params : vpic_params, vpic_harris_params {
};

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

struct user_global_t {
  struct params prm;
  struct globals_diag diag;
  struct globals_physics phys;
};

// ----------------------------------------------------------------------

void user_init(vpic_simulation *simulation, user_global_t *user_global, params *prm,
	       globals_physics *phys, globals_diag *diag);

void vpic_simulation_diagnostics(vpic_simulation *simulation, user_global_t *user_global,
				 params *prm, globals_diag *diag);

#endif
