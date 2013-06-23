
#include "ggcm_mhd_commu_private.h"

#include "ggcm_mhd.h"

#include <mrc_params.h>
#include <mrc_io.h>

// ======================================================================
// ggcm_mhd_commu class

// ----------------------------------------------------------------------
// ggcm_mhd_commu_run (forward -> subclass)

void
ggcm_mhd_commu_run(struct ggcm_mhd_commu *commu, struct mrc_fld *fld, int mb, int me)
{
  struct ggcm_mhd_commu_ops *ops = ggcm_mhd_commu_ops(commu);
  ops->run(commu, fld, mb, me);
}

// ----------------------------------------------------------------------
// ggcm_mhd_commu_init

static void
ggcm_mhd_commu_init()
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_commu,
			      &ggcm_mhd_commu_c_ops);
}

// ----------------------------------------------------------------------
// ggcm_mhd_commu description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_commu, x)
static struct param ggcm_mhd_commu_descr[] = {
  { "mhd"             , VAR(mhd)             , PARAM_OBJ(ggcm_mhd)     },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_commu class description

struct mrc_class_ggcm_mhd_commu mrc_class_ggcm_mhd_commu = {
  .name             = "ggcm_mhd_commu",
  .size             = sizeof(struct ggcm_mhd_commu),
  .param_descr      = ggcm_mhd_commu_descr,
  .init             = ggcm_mhd_commu_init,
};

