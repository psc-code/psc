
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
};

void vpic_marder_clean_div_e(struct vpic_marder *vmarder, struct vpic_mfields *vmflds,
			     struct vpic_mparticles *vmprts);
void vpic_marder_clean_div_b(struct vpic_marder *vmarder, struct vpic_mfields *vmflds);
void vpic_marder_sync_faces(struct vpic_marder *vmarder, struct vpic_mfields *vmflds);

#endif

