
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_io.h>

#include <assert.h>

// ======================================================================
// ggcm_mhd_ic class

// ----------------------------------------------------------------------
// ggcm_mhd_ic_init_masks_default

static void
ggcm_mhd_ic_init_masks_default(struct ggcm_mhd_ic *ic)
{
  struct mrc_fld *fld = ic->mhd->fld;

  mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {
    MRC_F3(fld, _YMASK, ix,iy,iz) = 1.;
  } mrc_fld_foreach_end;
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
// ggcm_mhd_ic_ini_b

void
ggcm_mhd_ic_ini_b(struct ggcm_mhd_ic *ic, float b_sw[3])
{
  assert(ic->mhd);
  struct ggcm_mhd_ic_ops *ops = ggcm_mhd_ic_ops(ic);
  assert(ops && ops->ini_b);
  ops->ini_b(ic, b_sw);
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

