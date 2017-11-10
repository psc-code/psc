
#ifndef PSC_MARDER_VPIC_H
#define PSC_MARDER_VPIC_H

#include "psc_marder_private.h"

// ======================================================================
// psc_marder "vpic"

struct psc_marder_vpic {
  struct vpic_marder *vmarder;
};

#define psc_marder_vpic(marder) mrc_to_subobj(marder, struct psc_marder_vpic)


#endif
