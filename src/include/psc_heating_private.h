
#ifndef PSC_HEATING_PRIVATE_H
#define PSC_HEATING_PRIVATE_H

#include "psc_heating.h"

struct psc_heating {
  struct mrc_obj obj;

  // params
  int tb; // in terms of step number (FIXME!)
  int te;
  int kind; // the kind of particle to be heated;
  int every_step; // heat every so many steps
};

#define VAR(x) (void *)offsetof(struct psc_heating, x)
static struct param psc_heating_descr[] _mrc_unused = {
  { "tb"                , VAR(tb)                , PARAM_INT(0.)          },
  { "te"                , VAR(te)                , PARAM_INT(0.)          },
  { "kind"              , VAR(kind)              , PARAM_INT(-1)          },
  { "every_step"        , VAR(every_step)        , PARAM_INT(20)          },
  {},
};
#undef VAR

struct psc_heating_ops {
  MRC_SUBCLASS_OPS(struct psc_heating);
  double (*get_H)(struct psc_heating *heating, double x[3]);
};

#define psc_heating_ops(heating) ((struct psc_heating_ops *)((heating)->obj.ops))

#endif
