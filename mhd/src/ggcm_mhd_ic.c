
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_io.h>
#include <mrc_fld_as_double.h>

#include <assert.h>

// ======================================================================
// ggcm_mhd_ic class

// ----------------------------------------------------------------------
// ggcm_mhd_ic_init_ymask_default

static void
ggcm_mhd_ic_init_ymask_default(struct ggcm_mhd_ic *ic, struct mrc_fld *ymask_base)
{
  struct mrc_fld *ymask = mrc_fld_get_as(ymask_base, FLD_TYPE);

  for (int p = 0; p < mrc_fld_nr_patches(ymask); p++) {
    mrc_fld_foreach(ymask, ix, iy, iz, 2, 2) {
      M3(ymask, 0, ix,iy,iz, p) = 1.;
    } mrc_fld_foreach_end;
  }
    
  mrc_fld_put_as(ymask, ymask_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_run

void
ggcm_mhd_ic_run(struct ggcm_mhd_ic *ic)
{
  assert(ic->mhd);
  struct ggcm_mhd_ic_ops *ops = ggcm_mhd_ic_ops(ic);
  assert(ops && ops->run);
  ops->run(ic);

  if (ops->init_ymask) {
    assert(ic->mhd->ymask);
    ops->init_ymask(ic, ic->mhd->ymask);
  } else {
    if (ic->mhd->ymask) {
      ggcm_mhd_ic_init_ymask_default(ic, ic->mhd->ymask);
    }
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_init

static void
ggcm_mhd_ic_init()
{
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic, x)
static struct param ggcm_mhd_ic_descr[] = {
  { "mhd"             , VAR(mhd)             , PARAM_OBJ(ggcm_mhd)      },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic class description

struct mrc_class_ggcm_mhd_ic mrc_class_ggcm_mhd_ic = {
  .name             = "ggcm_mhd_ic",
  .size             = sizeof(struct ggcm_mhd_ic),
  .param_descr      = ggcm_mhd_ic_descr,
  .init             = ggcm_mhd_ic_init,
};

