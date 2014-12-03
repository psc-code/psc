
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

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "vx:vy:vz");
  mrc_fld_set_type(fld_r, FLD_TYPE);
  mrc_fld_setup(fld_r);

  struct mrc_fld *r = mrc_fld_get_as(fld_r, FLD_TYPE);
  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
      float rri = 1.f / RR_(f, ix,iy,iz, p);
      M3(r, 0, ix,iy,iz, p) = rri * RVX_(f, ix,iy,iz, p);
      M3(r, 1, ix,iy,iz, p) = rri * RVY_(f, ix,iy,iz, p);
      M3(r, 2, ix,iy,iz, p) = rri * RVZ_(f, ix,iy,iz, p);
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(r, fld_r);
  mrc_fld_put_as(f, fld);

  float scale_vv = mhd->par.vvnorm;
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

  float scale_rr = mhd->par.rrnorm;
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

  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "pp");
  mrc_fld_set_type(fld_r, FLD_TYPE);
  mrc_fld_setup(fld_r);

  float gamm = mhd->par.gamm;

  struct mrc_fld *r = mrc_fld_get_as(fld_r, FLD_TYPE);
  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);

  if (mhd_type == MT_SEMI_CONSERVATIVE_GGCM ||
      mhd_type == MT_SEMI_CONSERVATIVE) {
    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	float rvv = (sqr(RVX_(f, ix,iy,iz, p)) +
		     sqr(RVY_(f, ix,iy,iz, p)) +
		     sqr(RVZ_(f, ix,iy,iz, p))) / RR_(f, ix,iy,iz, p);
	M3(r,0, ix,iy,iz, p) = (gamm - 1.f) * (UU_(f, ix,iy,iz, p) - .5f * rvv);
      } mrc_fld_foreach_end;
    }
  } else if (mhd_type == MT_FULLY_CONSERVATIVE) {
    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
	float rvv = (sqr(RVX_(f, ix,iy,iz, p)) +
		     sqr(RVY_(f, ix,iy,iz, p)) +
		     sqr(RVZ_(f, ix,iy,iz, p))) / RR_(f, ix,iy,iz, p);
	float b2  = (sqr(.5f * (BX_(f, ix,iy,iz, p) + BX_(f, ix+1,iy  ,iz  , p))) +
		     sqr(.5f * (BY_(f, ix,iy,iz, p) + BY_(f, ix  ,iy+1,iz  , p))) +
		     sqr(.5f * (BZ_(f, ix,iy,iz, p) + BZ_(f, ix  ,iy  ,iz+1, p))));
	M3(r,0, ix,iy,iz, p) = (gamm - 1.f) * (EE_(f, ix,iy,iz, p) - .5f * rvv - .5f * b2);
      } mrc_fld_foreach_end;
    }
  } else {
    assert(0);
  }

  mrc_fld_put_as(r, fld_r);
  mrc_fld_put_as(f, fld);

  ggcm_mhd_diag_c_write_one_field(io, fld_r, 0, "pp", mhd->par.ppnorm, diag_type, plane);

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

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "bx:by:bz");
  mrc_fld_set_type(fld_r, FLD_TYPE);
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
    compute_B_cc(fld_r, f, 0, 0);
  } else {
    assert(0);
  }

  mrc_fld_put_as(r, fld_r);
  mrc_fld_put_as(f, fld);

  float scale_bb = mhd->par.bbnorm;
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

