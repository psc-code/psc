
#ifndef PSC_HEATING_PRIVATE_H
#define PSC_HEATING_PRIVATE_H

#include "psc_heating.h"

struct psc_heating {
  struct mrc_obj obj;

  // params
  double zl; // in internal units (d_e)
  double zh;
  double xc;
  double yc;
  double rH;
  int tb; // in terms of step number (FIXME!)
  int te;
  double T;
  double Mi;
  int kind; // the kind of particle to be heated;
  int every_step; // heat every so many steps

  // state
  double fac;
};

#define VAR(x) (void *)offsetof(struct psc_heating, x)
static struct param psc_heating_descr[] _mrc_unused = {
  { "zl"                , VAR(zl)                , PARAM_DOUBLE(0.)       },
  { "zh"                , VAR(zh)                , PARAM_DOUBLE(0.)       },
  { "xc"                , VAR(xc)                , PARAM_DOUBLE(0.)       },
  { "yc"                , VAR(yc)                , PARAM_DOUBLE(0.)       },
  { "rH"                , VAR(rH)                , PARAM_DOUBLE(0.)       },
  { "tb"                , VAR(tb)                , PARAM_INT(0.)          },
  { "te"                , VAR(te)                , PARAM_INT(0.)          },
  { "T"                 , VAR(T)                 , PARAM_DOUBLE(.04)      },
  { "Mi"                , VAR(Mi)                , PARAM_DOUBLE(1.)       },
  { "kind"              , VAR(kind)              , PARAM_INT(-1)          },
  { "every_step"        , VAR(every_step)        , PARAM_INT(20)          },
  {},
};
#undef VAR

#endif
