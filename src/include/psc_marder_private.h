
#ifndef PSC_MARDER_PRIVATE_H
#define PSC_MARDER_PRIVATE_H

#include <psc_marder.h>

struct psc_marder {
  struct mrc_obj obj;
  // parameters
  int every_step; //< do Marder correction every so many steps
  double diffusion; //< diffusion coefficient for Marder correction
  int loop; //< execute this many relaxation steps in a loop
  bool dump; //< dump div_E, rho

  // state
  struct psc_mfields *div_e;
  struct psc_mfields *rho;
  struct psc_bnd *bnd; //< for filling ghosts on div_e
  struct mrc_io *io; //< for debug dumping
};

#endif
