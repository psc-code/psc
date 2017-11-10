
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

  int clean_div_e_interval;
  int clean_div_b_interval;
  int sync_shared_interval;

  // state
  struct psc_mfields *div_e;
  struct psc_mfields *rho;
  struct psc_bnd *bnd; //< for filling ghosts on rho, div_e
  struct psc_output_fields_item *item_div_e;
  struct psc_output_fields_item *item_rho;
  struct mrc_io *io; //< for debug dumping
};

struct psc_marder_ops {
  MRC_SUBCLASS_OPS(struct psc_marder);
  void (*correct)(struct psc_marder *marder, struct psc_mfields *mflds, struct psc_mfields *f);
};

#define psc_marder_ops(marder) ((struct psc_marder_ops *)(marder->obj.ops))

extern struct psc_marder_ops psc_marder_c_ops;
extern struct psc_marder_ops psc_marder_single_ops;
extern struct psc_marder_ops psc_marder_cuda_ops;

#endif
