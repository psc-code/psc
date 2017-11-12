
#ifndef VPIC_MARDER_H
#define VPIC_MARDER_H

#include "vpic_iface.h"

#include <vpic.h>

// ======================================================================
// vpic_marder

struct vpic_marder {
  int clean_div_e_interval;
  int clean_div_b_interval;
  int sync_shared_interval;
  int num_div_e_round;
  int num_div_b_round;

  void clean_div_e(vpic_mfields *vmflds, vpic_mparticles *vmprts);
  void clean_div_b(vpic_mfields *vmflds);
};

#endif

