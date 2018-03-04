
#include "psc_collision_private.h"

// ======================================================================
// psc_collision_init

extern struct psc_collision_ops psc_collision_fortran_ops;
extern struct psc_collision_ops psc_collision_single_ops;
extern struct psc_collision_ops psc_collision_double_ops;
extern struct psc_collision_ops psc_collision_vpic_ops;

static void
psc_collision_init()
{
  mrc_class_register_subclass(&mrc_class_psc_collision, &psc_collision_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_collision, &psc_collision_double_ops);
#ifdef USE_FORTRAN
  mrc_class_register_subclass(&mrc_class_psc_collision, &psc_collision_fortran_ops);
#endif
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_collision, &psc_collision_vpic_ops);
#endif
}

// ======================================================================
// psc_collision class

#define VAR(x) (void *)offsetof(struct psc_collision, x)
static struct param psc_collision_descr[] = {
  { "every"         , VAR(every)       , PARAM_INT(0)      },
  { "nu"            , VAR(nu)          , PARAM_DOUBLE(-1.) },
  {},
};
#undef VAR

struct mrc_class_psc_collision_ : mrc_class_psc_collision {
  mrc_class_psc_collision_() {
    name             = "psc_collision";
    size             = sizeof(struct psc_collision);
    param_descr      = psc_collision_descr;
    init             = psc_collision_init;
  }
} mrc_class_psc_collision;

