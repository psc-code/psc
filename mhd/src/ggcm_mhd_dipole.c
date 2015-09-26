
#include "ggcm_mhd_dipole_private.h"

#include "ggcm_mhd_private.h"

#include <assert.h>

// ======================================================================
// ggcm_mhd_dipole class

#define ggcm_mhd_dipole_ops(mhd_dipole) ((struct ggcm_mhd_dipole_ops *) mhd_dipole->obj.ops)

// ----------------------------------------------------------------------
// ggcm_mhd_dipole_setup

static void
_ggcm_mhd_dipole_setup(struct ggcm_mhd_dipole *mhd_dipole)
{
  struct ggcm_mhd *mhd = mhd_dipole->mhd;

  if (mhd) {
    mhd_dipole->bdip = ggcm_mhd_get_3d_fld(mhd_dipole->mhd, 3);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_dipole_destroy

static void
_ggcm_mhd_dipole_destroy(struct ggcm_mhd_dipole *mhd_dipole)
{
  struct ggcm_mhd *mhd = mhd_dipole->mhd;

  if (mhd) {
    ggcm_mhd_put_3d_fld(mhd_dipole->mhd, mhd_dipole->bdip);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_dipole_add_dipole
//
// bx: the field to put bx... by->bx+1 and bz->bx+2
// x0: dipole location
// moment: magnetic moment
// xmir: if x < xmir then A = 0... this is turned off if xmir == 0.0
// keep: B = keep * B + curl A

void
ggcm_mhd_dipole_add_dipole(struct ggcm_mhd_dipole *mhd_dipole, struct mrc_fld *b,
			   float x0[3], float moment[3], float xmir, float keep)
{
  struct ggcm_mhd_dipole_ops *ops = ggcm_mhd_dipole_ops(mhd_dipole);
  assert(ops && ops->add_dipole);
  ops->add_dipole(mhd_dipole, b, x0, moment, xmir, keep);
}

// ----------------------------------------------------------------------
// ggcm_mhd_dipole_set_b_field

void
ggcm_mhd_dipole_set_b_field(struct ggcm_mhd_dipole *mhd_dipole,
			    float moment[3], double diptime)
{
  struct ggcm_mhd_dipole_ops *ops = ggcm_mhd_dipole_ops(mhd_dipole);
  assert(ops);
  if (ops->set_b_field) {
    ops->set_b_field(mhd_dipole, mhd_dipole->bdip, moment, diptime);
  } else {
    ggcm_mhd_dipole_add_dipole(mhd_dipole, mhd_dipole->bdip, (float [3]) { 0., 0., 0. },
			       moment, 0., 0.);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_dipole_update_b_field

void
ggcm_mhd_dipole_update_b_field(struct ggcm_mhd_dipole *mhd_dipole,
			       struct mrc_fld *fld, double dacttime)
{
  struct ggcm_mhd_dipole_ops *ops = ggcm_mhd_dipole_ops(mhd_dipole);
  assert(ops && ops->update_b_field);
  ops->update_b_field(mhd_dipole, mhd_dipole->bdip, fld, dacttime);
}

// ----------------------------------------------------------------------
// ggcm_mhd_dipole_init

static void
ggcm_mhd_dipole_init()
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_dipole, &ggcm_mhd_dipole_none_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_dipole, &ggcm_mhd_dipole_float_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_dipole, &ggcm_mhd_dipole_double_ops);
}

// ----------------------------------------------------------------------
// ggcm_mhd_dipole description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_dipole, x)
static struct param ggcm_mhd_dipole_descr[] = {
  { "mhd"                , VAR(mhd)                , PARAM_OBJ(ggcm_mhd)    },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_dipole class

struct mrc_class_ggcm_mhd_dipole mrc_class_ggcm_mhd_dipole = {
  .name             = "ggcm_mhd_dipole",
  .size             = sizeof(struct ggcm_mhd_dipole),
  .param_descr      = ggcm_mhd_dipole_descr,
  .init             = ggcm_mhd_dipole_init,
  .setup            = _ggcm_mhd_dipole_setup,
  .destroy          = _ggcm_mhd_dipole_destroy,
};

