
#ifndef VPIC_DIAG_H
#define VPIC_DIAG_H

#include "vpic_iface.h"
#include "vpic.h"

// ----------------------------------------------------------------------
// VpicDiag

struct VpicDiag {
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

  VpicDiag(vpic_simulation* simulaton, int interval_);
  void setup();
  void run();

private:
  vpic_simulation* simulation_;
};

#endif
