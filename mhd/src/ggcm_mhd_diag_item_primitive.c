
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag_private.h"

#include <mrc_domain.h>
#include <stdio.h>
#include <assert.h>

// ======================================================================
// ggcm_mhd_diag_item subclass "v"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_v_run

static void
ggcm_mhd_diag_item_v_run(struct ggcm_mhd_diag_item *item,
			 struct mrc_io *io, struct mrc_fld *f,
			 int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;

  primvar1_c(mhd);
  float scale_vv = mhd->par.vvnorm;
  ggcm_mhd_diag_c_write_one_field(io, f, _VX, "vx", scale_vv, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, f, _VY, "vy", scale_vv, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, f, _VZ, "vz", scale_vv, diag_type, plane);
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
			  struct mrc_io *io, struct mrc_fld *f,
			  int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;

  primvar1_c(mhd);
  float scale_rr = mhd->par.rrnorm;
  ggcm_mhd_diag_c_write_one_field(io, f, _RR, "rr", scale_rr, diag_type, plane);
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

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "pp");
  mrc_fld_setup(fld_r);

  float gamm = mhd->par.gamm;

  struct mrc_fld *r = mrc_fld_get_as(fld_r, "float");
  struct mrc_fld *f = mrc_fld_get_as(fld, "float");

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    float rvv = (sqr(MRC_F3(f, _RV1X, ix,iy,iz)) +
		 sqr(MRC_F3(f, _RV1Y, ix,iy,iz)) +
		 sqr(MRC_F3(f, _RV1Z, ix,iy,iz))) / MRC_F3(f, _RR1, ix,iy,iz);
    MRC_F3(r,0, ix,iy,iz) = (gamm - 1.f) * (MRC_F3(f,_UU1, ix,iy,iz) - .5f * rvv);
  } mrc_fld_foreach_end;

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
// ggcm_mhd_diag_item subclass "pp_full"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_pp_full_run

#define F3 MRC_F3 // FIXME

static void
ggcm_mhd_diag_item_pp_full_run(struct ggcm_mhd_diag_item *item,
			       struct mrc_io *io, struct mrc_fld *fld,
			       int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "pp_full");
  mrc_fld_setup(fld_r);

  float gamm = mhd->par.gamm;

  struct mrc_fld *r = mrc_fld_get_as(fld_r, "float");
  struct mrc_fld *f = mrc_fld_get_as(fld, "mhd_fc_float");

  mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
    float rvv = (sqr(MRC_F3(f, _RV1X, ix,iy,iz)) +
		 sqr(MRC_F3(f, _RV1Y, ix,iy,iz)) +
		 sqr(MRC_F3(f, _RV1Z, ix,iy,iz))) / MRC_F3(f, _RR1, ix,iy,iz);
    float b2  = (sqr(.5*(B1X(f, ix,iy,iz) + B1X(f, ix+1,iy  ,iz  ))) +
		 sqr(.5*(B1Y(f, ix,iy,iz) + B1Y(f, ix  ,iy+1,iz  ))) +
		 sqr(.5*(B1Z(f, ix,iy,iz) + B1Z(f, ix  ,iy  ,iz+1))));
    MRC_F3(r,0, ix,iy,iz) = (gamm - 1.f) * (MRC_F3(f,_UU1, ix,iy,iz) - .5f * rvv - .5f * b2);
  } mrc_fld_foreach_end;

  mrc_fld_put_as(r, fld_r);
  mrc_fld_put_as(f, fld);

  ggcm_mhd_diag_c_write_one_field(io, fld_r, 0, "pp_full", mhd->par.ppnorm, diag_type, plane);

  mrc_fld_destroy(fld_r);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "pp_full"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_pp_full = {
  .name             = "pp_full",
  .run              = ggcm_mhd_diag_item_pp_full_run,
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

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "pp_full");
  mrc_fld_set_param_int(fld_r, "nr_comps", 3);
  mrc_fld_setup(fld_r);

  struct mrc_fld *r = mrc_fld_get_as(fld_r, "float");
  struct mrc_fld *f = mrc_fld_get_as(fld, "float");

  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    MRC_F3(fld_r, 0, ix,iy,iz) = .5f*(MRC_F3(f,_B1X, ix,iy,iz) +
				      MRC_F3(f,_B1X, ix-1,iy,iz));
    MRC_F3(fld_r, 1, ix,iy,iz) = .5f*(MRC_F3(f,_B1Y, ix,iy,iz) +
				      MRC_F3(f,_B1Y, ix,iy-1,iz));
    MRC_F3(fld_r, 2, ix,iy,iz) = .5f*(MRC_F3(f,_B1Z, ix,iy,iz) +
				      MRC_F3(f,_B1Z, ix,iy,iz-1));
  } mrc_fld_foreach_end;

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


