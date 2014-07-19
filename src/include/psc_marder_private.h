
#ifndef PSC_MARDER_PRIVATE_H
#define PSC_MARDER_PRIVATE_H

#include <psc_marder.h>

struct psc_marder {
  struct mrc_obj obj;
  // parameters
  int every_step; //< do Marder correction every so many steps
  double diffusion; //< diffusion coefficient for Marder correction
};

#endif
