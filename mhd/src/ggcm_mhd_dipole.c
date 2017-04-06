
#include "ggcm_mhd_dipole_private.h"

#include "ggcm_mhd_private.h"

#include <mrc_physics.h>

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
// ggcm_mhd_dipole_vector_potential

double
ggcm_mhd_dipole_vector_potential(struct ggcm_mhd_dipole *mhd_dipole, int m,
				 double x[3], float x0[3], float moment[3], float xmir)
{
  struct ggcm_mhd_dipole_ops *ops = ggcm_mhd_dipole_ops(mhd_dipole);
  assert(ops && ops->vector_potential);
  return ops->vector_potential(mhd_dipole, m, x, x0, moment, xmir);
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
  { "r1lim"              , VAR(r1lim)              , PARAM_DOUBLE(1.5)      },
  // dipolestrength in external units, nT by default
  { "dipolestrength"     , VAR(dipolestrength)     , PARAM_DOUBLE(C_DIPOLESTRENGTH / 1e-9)   },
  { "dipolestrength_r"   , VAR(dipolestrength_r)   , PARAM_DOUBLE(1.)       },
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

/////////////////////////////////////////////////////////////////////////
// diag items that go with ggcm_mhd_dipole

#include "ggcm_mhd_diag_private.h"
#include "ggcm_mhd_diag_item_private.h"

#include <mrc_domain.h>
#include <mrc_fld_as_float.h>

// ======================================================================
// ggcm_mhd_diag_item subclass "bdip"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_bdip_run

static void
ggcm_mhd_diag_item_bdip_run(struct ggcm_mhd_diag_item *item,
			    struct mrc_io *io, struct mrc_fld *f,
			    int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;
  struct ggcm_mhd_dipole *mhd_dipole = ggcm_mhd_get_var_obj(mhd, "mhd_dipole");
  struct mrc_fld *bdip = mhd_dipole->bdip;

  ggcm_mhd_diag_c_write_one_field(io, bdip, 0, "bdipx", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, bdip, 1, "bdipy", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, bdip, 2, "bdipz", 1., diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "bdip"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_bdip = {
  .name             = "bdip",
  .run              = ggcm_mhd_diag_item_bdip_run,
};

// ======================================================================
// ggcm_mhd_diag_item subclass "bdipcc"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_bdipcc_run

static void
ggcm_mhd_diag_item_bdipcc_run(struct ggcm_mhd_diag_item *item,
			      struct mrc_io *io, struct mrc_fld *fld,
			      int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;
  struct ggcm_mhd_dipole *mhd_dipole = ggcm_mhd_get_var_obj(mhd, "mhd_dipole");

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "bdipccx:bdipccy:bbdipccz");
  mrc_fld_set_type(fld_r, FLD_TYPE);
  mrc_fld_setup(fld_r);

  struct mrc_fld *r = mrc_fld_get_as(fld_r, FLD_TYPE);
  struct mrc_fld *bdip = mrc_fld_get_as(mhd_dipole->bdip, FLD_TYPE);

  int mhd_type;
  mrc_fld_get_param_int(mhd->fld, "mhd_type", &mhd_type);
  if (MT_BGRID(mhd_type) == MT_BGRID_FC_GGCM) {
    for (int p = 0; p < mrc_fld_nr_patches(r); p++) {
      mrc_fld_foreach(r, ix,iy,iz, 0, 0) {
	M3(r, 0, ix,iy,iz, p) = .5f * (M3(bdip, 0, ix,iy,iz, p) + M3(bdip, 0, ix-1,iy,iz, p));
	M3(r, 1, ix,iy,iz, p) = .5f * (M3(bdip, 1, ix,iy,iz, p) + M3(bdip, 1, ix,iy-1,iz, p));
	M3(r, 2, ix,iy,iz, p) = .5f * (M3(bdip, 2, ix,iy,iz, p) + M3(bdip, 2, ix,iy,iz-1, p));
      } mrc_fld_foreach_end;
    }
  } else if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
    for (int p = 0; p < mrc_fld_nr_patches(r); p++) {
      mrc_fld_foreach(r, ix,iy,iz, 0, 0) {
	M3(r, 0, ix,iy,iz, p) = .5f * (M3(bdip, 0, ix,iy,iz, p) + M3(bdip, 0, ix+1,iy,iz, p));
	M3(r, 1, ix,iy,iz, p) = .5f * (M3(bdip, 1, ix,iy,iz, p) + M3(bdip, 1, ix,iy+1,iz, p));
	M3(r, 2, ix,iy,iz, p) = .5f * (M3(bdip, 2, ix,iy,iz, p) + M3(bdip, 2, ix,iy,iz+1, p));
      } mrc_fld_foreach_end;
    }
  } else {
    assert(0);
  }

  mrc_fld_put_as(r, fld_r);
  mrc_fld_put_as(bdip, mhd_dipole->bdip);

  ggcm_mhd_diag_c_write_one_field(io, fld_r, 0, "bdipccx", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld_r, 1, "bdipccy", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld_r, 2, "bdipccz", 1., diag_type, plane);

  mrc_fld_destroy(fld_r);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "bdipcc"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_bdipcc = {
  .name             = "bdipcc",
  .run              = ggcm_mhd_diag_item_bdipcc_run,
};

