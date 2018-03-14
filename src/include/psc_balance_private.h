
#ifndef PSC_BALANCE_PRIVATE_H
#define PSC_BALANCE_PRIVATE_H

#include <psc_balance.h>

struct psc_balance {
  struct mrc_obj obj;
  int every;
  double factor_fields;
  bool print_loads;
  bool write_loads;
};

struct psc_balance_ops {
  MRC_SUBCLASS_OPS(struct psc_balance);
};

#endif
