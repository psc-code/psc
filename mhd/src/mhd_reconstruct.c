
#include "mhd_reconstruct_private.h"

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "mhd_util.h"

#include <mrc_fld_as_double_aos.h>
#define F1(f, m, i) MRC_D2(f, m, i)
#include <assert.h>

// ======================================================================
// mhd_reconstruct class

// ----------------------------------------------------------------------
// mhd_reconstruct_run

void
mhd_reconstruct_run(struct mhd_reconstruct *mr,
		    struct mrc_fld *Ul, struct mrc_fld *Ur,
		    struct mrc_fld *Wl, struct mrc_fld *Wr,
		    struct mrc_fld *W1d, struct mrc_fld *Bxi,
		    int ldim, int l, int r, int dir)
{
  struct mhd_reconstruct_ops *ops = mhd_reconstruct_ops(mr);
  assert(ops && ops->run);

  ops->run(mr, Ul, Ur, Wl, Wr, W1d, Bxi, ldim, l, r, dir);
}

// ----------------------------------------------------------------------
// mhd_reconstruct_init

static void
mhd_reconstruct_init()
{
  mrc_class_register_subclass(&mrc_class_mhd_reconstruct, &mhd_reconstruct_pcm_double_ops);
  mrc_class_register_subclass(&mrc_class_mhd_reconstruct, &mhd_reconstruct_pcm_float_ops);
  mrc_class_register_subclass(&mrc_class_mhd_reconstruct, &mhd_reconstruct_plm_double_ops);
  mrc_class_register_subclass(&mrc_class_mhd_reconstruct, &mhd_reconstruct_plm_float_ops);
}

// ----------------------------------------------------------------------
// mhd_reconstruct description

#define VAR(x) (void *)offsetof(struct mhd_reconstruct, x)
static struct param mhd_reconstruct_descr[] = {
  { "mhd"             , VAR(mhd)             , PARAM_OBJ(ggcm_mhd)            },

  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mhd_reconstruct class description

struct mrc_class_mhd_reconstruct mrc_class_mhd_reconstruct = {
  .name             = "mhd_reconstruct",
  .size             = sizeof(struct mhd_reconstruct),
  .param_descr      = mhd_reconstruct_descr,
  .init             = mhd_reconstruct_init,
};

