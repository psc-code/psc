#include <bitfield3d.h>
#include "psc_domainwindow.h"

struct psc_domainwindow {
  struct mrc_obj obj;
  int np[3];
  struct bitfield3d activepatches;
};

///Interface for custom cases
///
///Take a look at the existing test-cases to see how to overload them.
///@note 
///When you have added an implementation of psc_case_ops, add its declaration
///as extern to the end of psc_case_private.h and psc_case.c
///@sa \ref custom_case
struct psc_domainwindow_ops {
  MRC_SUBCLASS_OPS(struct psc_domainwindow);
  void (*timestep)(struct psc_domainwindow* this, int timestep, double t);
};

extern struct psc_domainwindow_ops psc_movingwindow_z_ops;