
#ifndef PSC_HEATING_PRIVATE_H
#define PSC_HEATING_PRIVATE_H

#include "psc_heating.h"

#include <psc_heating_spot.h>

struct psc_heating {
  struct mrc_obj obj;

  // params
  int tb; // in terms of step number (FIXME!)
  int te;
  int kind; // the kind of particle to be heated
  int every_step; // heat every so many steps

  // state
  struct psc_heating_spot *spot;
};

#define VAR(x) (void *)offsetof(struct psc_heating, x)
static struct param psc_heating_descr[] _mrc_unused = {
  { "tb"                , VAR(tb)                , PARAM_INT(0)                  },
  { "te"                , VAR(te)                , PARAM_INT(0)                  },
  { "kind"              , VAR(kind)              , PARAM_INT(-1)                 },
  { "every_step"        , VAR(every_step)        , PARAM_INT(20)                 },

  { "spot"              , VAR(spot)              , MRC_VAR_OBJ(psc_heating_spot) },
  {},
};
#undef VAR

struct psc_heating_ops {
  MRC_SUBCLASS_OPS(struct psc_heating);
  void (*run)(struct psc_heating *heating, struct psc_mparticles *mprts_base,
	      struct psc_mfields *mflds_base);
};

#define psc_heating_ops(heating) ((struct psc_heating_ops *)((heating)->obj.ops))

extern struct psc_heating_ops psc_heating_ops_single;
extern struct psc_heating_ops psc_heating_ops_cuda;

#endif
