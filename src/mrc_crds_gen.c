
#include "mrc_crds_gen_private.h"

#include <mrc_crds.h>
#include <mrc_params.h>
#include <assert.h>

// ======================================================================
// mrc_crds_gen

// ----------------------------------------------------------------------
// mrc_crds_gen_run
//
// given parameters, including the number of grid points n,
// create a non-uniform grid in xx, with spacing in dx
// xx and dx are arrays that need to support indexing from -2..n+1

void
mrc_crds_gen_run(struct mrc_crds_gen *gen, double *xx, double *dx)
{
  struct mrc_crds_gen_ops *ops = mrc_crds_gen_ops(gen);
  assert(ops && ops->run);
  assert(gen->d >= 0 && gen->d < 3);
  ops->run(gen, xx, dx);
}

// ----------------------------------------------------------------------
// mrc_crds_gen_init

static void
mrc_crds_gen_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_crds_gen, &mrc_crds_gen_uniform_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds_gen, &mrc_crds_gen_ggcm_x_tanh_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds_gen, &mrc_crds_gen_ggcm_x_cubic_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds_gen, &mrc_crds_gen_ggcm_yz_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds_gen, &mrc_crds_gen_gaussian_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds_gen, &mrc_crds_gen_two_gaussian_ops);
}

// ----------------------------------------------------------------------
// mrc_crds_gen description

#define VAR(x) (void *)offsetof(struct mrc_crds_gen, x)
static struct param mrc_crds_gen_descr[] = {
  { "crds"            , VAR(crds)            , PARAM_OBJ(mrc_crds)   },
  { "d"               , VAR(d)               , PARAM_INT(-1)         },
  { "n"               , VAR(n)               , PARAM_INT(0)          },
  { "sw"              , VAR(sw)              , PARAM_INT(0)          },
  { "xl"              , VAR(xl)              , PARAM_DOUBLE(0.0)      },
  { "xh"              , VAR(xh)              , PARAM_DOUBLE(0.0)      },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mrc_crds_gen class

struct mrc_class_mrc_crds_gen mrc_class_mrc_crds_gen = {
  .name             = "mrc_crds_gen",
  .size             = sizeof(struct mrc_crds_gen),
  .param_descr      = mrc_crds_gen_descr,
  .init             = mrc_crds_gen_init,
};

