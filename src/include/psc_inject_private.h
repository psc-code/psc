
#ifndef PSC_INJECT_PRIVATE_H
#define PSC_INJECT_PRIVATE_H

#include <psc_inject.h>
#include <psc_target.h>
#include <psc_bnd.h>
#include <psc_output_fields_item.h>

struct psc_inject {
  struct mrc_obj obj;

  // params
  bool do_inject; // whether to inject particles at all
  int every_step; // inject every so many steps
  int tau; // in steps
  int kind_n; // the number of particles to inject are based on this kind's density
  struct psc_target *target;
};

#define VAR(x) (void *)offsetof(struct psc_inject, x)
static struct param psc_inject_descr[] _mrc_unused = {
  { "do_inject"  , VAR(do_inject)  , PARAM_BOOL(true)                    },
  { "every_step" , VAR(every_step) , PARAM_INT(20)                       },
  { "tau"        , VAR(tau)        , PARAM_INT(40)                       },
  { "kind_n"     , VAR(kind_n)     , PARAM_INT(-1)                       },
  { "target"     , VAR(target)     , PARAM_OBJ(psc_target)               },

  {},
};
#undef VAR

struct psc_inject_ops {
  MRC_SUBCLASS_OPS(struct psc_inject);
  void (*run)(struct psc_inject *inject, struct psc_mparticles *mprts_base,
   	      struct psc_mfields *mflds_base);
};

#define psc_inject_ops(inject) ((struct psc_inject_ops *)((inject)->obj.ops))

#endif

