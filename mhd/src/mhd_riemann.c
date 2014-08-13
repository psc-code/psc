
#include "mhd_riemann_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <math.h>

// ======================================================================
// mhd_riemann class

// ----------------------------------------------------------------------
// mhd_riemann_setup

static void
_mhd_riemann_setup(struct mhd_riemann *riemann)
{
  assert(riemann->mhd);
}

// ----------------------------------------------------------------------
// mhd_riemann_run

void
mhd_riemann_run(struct mhd_riemann *riemann,
		struct mrc_fld *F1d,
		struct mrc_fld *Ul, struct mrc_fld *Ur,
		struct mrc_fld *Wl, struct mrc_fld *Wr,
		int ldim, int l, int r, int dir)
{
  struct mhd_riemann_ops *ops = mhd_riemann_ops(riemann);
  assert(ops && ops->run);
  ops->run(riemann, F1d, Ul, Ur, Wl, Wr, ldim, l, r, dir);
}

// ----------------------------------------------------------------------
// mhd_riemann_init

static void
mhd_riemann_init()
{
  mrc_class_register_subclass(&mrc_class_mhd_riemann, &mhd_riemann_rusanov_double_ops);
  mrc_class_register_subclass(&mrc_class_mhd_riemann, &mhd_riemann_rusanov_float_ops);
  mrc_class_register_subclass(&mrc_class_mhd_riemann, &mhd_riemann_hll_ops);
  mrc_class_register_subclass(&mrc_class_mhd_riemann, &mhd_riemann_hlld_ops);
  mrc_class_register_subclass(&mrc_class_mhd_riemann, &mhd_riemann_hydro_rusanov_ops);
  mrc_class_register_subclass(&mrc_class_mhd_riemann, &mhd_riemann_hydro_hll_ops); 
  mrc_class_register_subclass(&mrc_class_mhd_riemann, &mhd_riemann_hydro_hllc_ops);
}

// ----------------------------------------------------------------------
// mhd_riemann description

#define VAR(x) (void *)offsetof(struct mhd_riemann, x)
static struct param mhd_riemann_descr[] = {
  { "mhd"             , VAR(mhd)             , PARAM_OBJ(ggcm_mhd)            },

  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mhd_riemann class description

struct mrc_class_mhd_riemann mrc_class_mhd_riemann = {
  .name             = "mhd_riemann",
  .size             = sizeof(struct mhd_riemann),
  .param_descr      = mhd_riemann_descr,
  .init             = mhd_riemann_init,
  .setup            = _mhd_riemann_setup,
};

