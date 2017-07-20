
#ifndef PSC_TARGET_PRIVATE_H
#define PSC_TARGET_PRIVATE_H

#include <psc_target.h>

struct psc_target {
  struct mrc_obj obj;
  
  // params
  double yl;
  double yh;
  double zl;
  double zh;
  double n;
  double Te;
  double Ti;
  int kind_ion;
  int kind_electron;
};

#define VAR(x) (void *)offsetof(struct psc_target, x)
static struct param psc_target_descr[] _mrc_unused = {
  { "yl"           , VAR(yl)           , PARAM_DOUBLE(0.)       },
  { "yh"           , VAR(yh)           , PARAM_DOUBLE(0.)       },
  { "zl"           , VAR(zl)           , PARAM_DOUBLE(0.)       },
  { "zh"           , VAR(zh)           , PARAM_DOUBLE(0.)       },
  { "n"            , VAR(n)            , PARAM_DOUBLE(1.)       },
  { "Te"           , VAR(Te)           , PARAM_DOUBLE(.001)     },
  { "Ti"           , VAR(Ti)           , PARAM_DOUBLE(.001)     },
  { "kind_ion"     , VAR(kind_ion)     , PARAM_INT(-1)          },
  { "kind_electron", VAR(kind_electron), PARAM_INT(-1)          },
  {},
};
#undef VAR

#endif


