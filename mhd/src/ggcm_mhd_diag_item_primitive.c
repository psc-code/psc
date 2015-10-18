
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag_private.h"

#include <mrc_fld_as_double.h>
#include <mrc_domain.h>
#include <stdio.h>
#include <assert.h>

#include "mhd_3d.c"

// ======================================================================
// ggcm_mhd_diag_item subclass "v"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_v_run

static void
ggcm_mhd_diag_item_v_run(struct ggcm_mhd_diag_item *item,
			 struct mrc_io *io, struct mrc_fld *fld,
			 int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;

  int bnd = fld->_nr_ghosts;

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, bnd, "vx:vy:vz");
  mrc_fld_set_type(fld_r, FLD_TYPE);
  mrc_fld_setup(fld_r);

  struct mrc_fld *r = mrc_fld_get_as(fld_r, FLD_TYPE);
  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(f, ix,iy,iz, bnd, bnd) {
      mrc_fld_data_t rri = 1.f / RR_(f, ix,iy,iz, p);
      M3(r, 0, ix,iy,iz, p) = rri * RVX_(f, ix,iy,iz, p);
      M3(r, 1, ix,iy,iz, p) = rri * RVY_(f, ix,iy,iz, p);
      M3(r, 2, ix,iy,iz, p) = rri * RVZ_(f, ix,iy,iz, p);
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(r, fld_r);
  mrc_fld_put_as(f, fld);

  mrc_fld_data_t scale_vv = mhd->vvnorm;
  ggcm_mhd_diag_c_write_one_field(io, fld_r, 0, "vx", scale_vv, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld_r, 1, "vy", scale_vv, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld_r, 2, "vz", scale_vv, diag_type, plane);

  mrc_fld_destroy(fld_r);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "v"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_v = {
  .name             = "v",
  .run              = ggcm_mhd_diag_item_v_run,
};

// ======================================================================
// ggcm_mhd_diag_item subclass "rr"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_rr_run

static void
ggcm_mhd_diag_item_rr_run(struct ggcm_mhd_diag_item *item,
			  struct mrc_io *io, struct mrc_fld *fld,
			  int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;

  mrc_fld_data_t scale_rr = mhd->rrnorm;
  ggcm_mhd_diag_c_write_one_field(io, fld, RR, "rr", scale_rr, diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "rr"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_rr = {
  .name             = "rr",
  .run              = ggcm_mhd_diag_item_rr_run,
};

// ======================================================================
// ggcm_mhd_diag_item subclass "pp"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_pp_run

static void
ggcm_mhd_diag_item_pp_run(struct ggcm_mhd_diag_item *item,
			  struct mrc_io *io, struct mrc_fld *fld,
			  int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;

  int gdims[3];
  mrc_domain_get_global_dims(fld->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);

  int bnd = fld->_nr_ghosts - 1;

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, bnd + 1, "pp");
  mrc_fld_set_type(fld_r, FLD_TYPE);
  mrc_fld_setup(fld_r);

  mrc_fld_data_t gamm = mhd->par.gamm;

  struct mrc_fld *r = mrc_fld_get_as(fld_r, FLD_TYPE);
  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);

  if (mhd_type == MT_SEMI_CONSERVATIVE_GGCM ||
      mhd_type == MT_SEMI_CONSERVATIVE) {
    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, bnd, bnd) {
	mrc_fld_data_t rvv = (sqr(RVX_(f, ix,iy,iz, p)) +
			      sqr(RVY_(f, ix,iy,iz, p)) +
			      sqr(RVZ_(f, ix,iy,iz, p))) / RR_(f, ix,iy,iz, p);
	M3(r,0, ix,iy,iz, p) = (gamm - 1.f) * (UU_(f, ix,iy,iz, p) - .5f * rvv);
      } mrc_fld_foreach_end;
    }
  } else if (mhd_type == MT_FULLY_CONSERVATIVE) {
    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, bnd, bnd) {
	mrc_fld_data_t rvv = (sqr(RVX_(f, ix,iy,iz, p)) +
			      sqr(RVY_(f, ix,iy,iz, p)) +
			      sqr(RVZ_(f, ix,iy,iz, p))) / RR_(f, ix,iy,iz, p);
	mrc_fld_data_t b2  = (sqr(.5f * (BX_(f, ix,iy,iz, p) + BX_(f, ix+dx,iy   ,iz   , p))) +
			      sqr(.5f * (BY_(f, ix,iy,iz, p) + BY_(f, ix   ,iy+dy,iz   , p))) +
			      sqr(.5f * (BZ_(f, ix,iy,iz, p) + BZ_(f, ix   ,iy   ,iz+dz, p))));
	M3(r,0, ix,iy,iz, p) = (gamm - 1.f) * (EE_(f, ix,iy,iz, p) - .5f * rvv - .5f * b2);
      } mrc_fld_foreach_end;
    }
  } else {
    assert(0);
  }

  mrc_fld_put_as(r, fld_r);
  mrc_fld_put_as(f, fld);

  ggcm_mhd_diag_c_write_one_field(io, fld_r, 0, "pp", mhd->ppnorm, diag_type, plane);

  mrc_fld_destroy(fld_r);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "pp"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_pp = {
  .name             = "pp",
  .run              = ggcm_mhd_diag_item_pp_run,
};

// ======================================================================
// ggcm_mhd_diag_item subclass "b" (cell centered b)

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_b_run

static void
ggcm_mhd_diag_item_b_run(struct ggcm_mhd_diag_item *item,
			 struct mrc_io *io, struct mrc_fld *fld,
			 int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;

  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);

  int bnd = fld->_nr_ghosts - 1;

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "bx:by:bz");
  mrc_fld_set_type(fld_r, FLD_TYPE);
  mrc_fld_set_param_int(fld_r, "nr_ghosts", bnd + 1);
  mrc_fld_setup(fld_r);

  struct mrc_fld *r = mrc_fld_get_as(fld_r, FLD_TYPE);
  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);
  if (mhd_type == MT_SEMI_CONSERVATIVE_GGCM) {
    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
	M3(r, 0, ix,iy,iz, p) = .5f * (BX_(f, ix,iy,iz, p) + BX_(f, ix-1,iy,iz, p));
	M3(r, 1, ix,iy,iz, p) = .5f * (BY_(f, ix,iy,iz, p) + BY_(f, ix,iy-1,iz, p));
	M3(r, 2, ix,iy,iz, p) = .5f * (BZ_(f, ix,iy,iz, p) + BZ_(f, ix,iy,iz-1, p));
      } mrc_fld_foreach_end;
    }
  } else if (mhd_type == MT_SEMI_CONSERVATIVE ||
	     mhd_type == MT_FULLY_CONSERVATIVE) {
    compute_Bt_cc(mhd, fld_r, f, bnd, bnd);
  } else {
    assert(0);
  }

  mrc_fld_put_as(r, fld_r);
  mrc_fld_put_as(f, fld);

  mrc_fld_data_t scale_bb = mhd->bbnorm;
  ggcm_mhd_diag_c_write_one_field(io, fld_r, 0, "bx", scale_bb, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld_r, 1, "by", scale_bb, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld_r, 2, "bz", scale_bb, diag_type, plane);

  mrc_fld_destroy(fld_r);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "b"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_b = {
  .name             = "b",
  .run              = ggcm_mhd_diag_item_b_run,
};

