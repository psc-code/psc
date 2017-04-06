
#include "ggcm_mhd_diag_private.h"

#include "ggcm_mhd.h"

#include <mrc_io.h>
#include <assert.h>

// ======================================================================
// ggcm_mhd_diag class

#define ggcm_mhd_diag_ops(diag) ((struct ggcm_mhd_diag_ops *)(diag->obj.ops))

// ----------------------------------------------------------------------
// ggcm_mhd_diag_run

void
ggcm_mhd_diag_run(struct ggcm_mhd_diag *diag)
{
  struct ggcm_mhd_diag_ops *ops = ggcm_mhd_diag_ops(diag);
  assert(ops && ops->run);
  ops->run(diag);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_run_now

void
ggcm_mhd_diag_run_now(struct ggcm_mhd_diag *diag, struct mrc_fld *fld,
			int diag_type, int itdia)
{
  struct ggcm_mhd_diag_ops *ops = ggcm_mhd_diag_ops(diag);
  assert(ops && ops->run);
  ops->run_now(diag, fld, diag_type, itdia);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_shutdown

void
ggcm_mhd_diag_shutdown(struct ggcm_mhd_diag *diag)
{
  struct ggcm_mhd_diag_ops *ops = ggcm_mhd_diag_ops(diag);
  assert(ops && ops->shutdown);
  ops->shutdown(diag);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_init

static void
ggcm_mhd_diag_init()
{
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_mod_register

void
ggcm_mhd_diag_mod_register(struct ggcm_mhd_diag *mhd_diag, struct mrc_mod *mod)
{
  struct ggcm_mhd_diag_ops *ops = ggcm_mhd_diag_ops(mhd_diag);

  if (ops && ops->mod_register) {
    ops->mod_register(mhd_diag, mod);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_diag, x)
static struct param ggcm_mhd_diag_descr[] = {
  { "mhd"             , VAR(mhd)             , PARAM_OBJ(ggcm_mhd)      },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_diag class

struct mrc_class_ggcm_mhd_diag mrc_class_ggcm_mhd_diag = {
  .name             = "ggcm_mhd_diag",
  .size             = sizeof(struct ggcm_mhd_diag),
  .param_descr      = ggcm_mhd_diag_descr,
  .init             = ggcm_mhd_diag_init,
};

