
#ifndef VPIC_DIAG_H
#define VPIC_DIAG_H

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

  globals_diag(int interval_);
  void setup(Simulation *sim);
  void run(Simulation *sim);
};

#endif
