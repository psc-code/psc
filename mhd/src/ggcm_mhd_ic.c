
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_private.h"

#include <mrc_io.h>
#include <mrc_fld_as_double_aos.h>

#include <assert.h>

// ======================================================================
// ggcm_mhd_ic class

// ----------------------------------------------------------------------
// ggcm_mhd_ic_init_masks_default

static void
ggcm_mhd_ic_init_masks_default(struct ggcm_mhd_ic *ic)
{
  struct mrc_fld *fld = mrc_fld_get_as(ic->mhd->fld, FLD_TYPE);

  mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {
    F3(fld, _YMASK, ix,iy,iz) = 1.;
  } mrc_fld_foreach_end;
    
  mrc_fld_put_as(fld, ic->mhd->fld);
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

  if (ops->init_masks) {
    ops->init_masks(ic);
  } else {
    ggcm_mhd_ic_init_masks_default(ic);
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

