
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag_private.h"

#include <mrc_fld_as_double.h>
#include <mrc_domain.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

// ======================================================================
// ggcm_mhd_diag_item subclass "rr1"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_rr1_run

static void
ggcm_mhd_diag_item_rr1_run(struct ggcm_mhd_diag_item *item,
			   struct mrc_io *io, struct mrc_fld *fld,
			   int diag_type, float plane)
{
  float scale_rr = 1.;
  ggcm_mhd_diag_c_write_one_field(io, fld, RR, "rr1", scale_rr, diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "rr1"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_rr1 = {
  .name             = "rr1",
  .run              = ggcm_mhd_diag_item_rr1_run,
};

// ======================================================================
// ggcm_mhd_diag_item subclass "uu1"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_uu1_run

static void
ggcm_mhd_diag_item_uu1_run(struct ggcm_mhd_diag_item *item,
			   struct mrc_io *io, struct mrc_fld *fld,
			   int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;

  int gdims[3];
  mrc_domain_get_global_dims(fld->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);

  int bnd = fld->_nr_ghosts;

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "uu1");
  mrc_fld_set_type(fld_r, FLD_TYPE);
  mrc_fld_set_param_int(fld_r, "nr_ghosts", bnd);
  mrc_fld_setup(fld_r);

  struct mrc_fld *r = mrc_fld_get_as(fld_r, FLD_TYPE);
  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);

  if (MT_FORMULATION(mhd_type) == MT_FORMULATION_SCONS) {
    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, bnd, bnd) {
	M3(r, 0, ix,iy,iz, p) = UU_(f, ix,iy,iz, p);
      } mrc_fld_foreach_end;
    }
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_FCONS) {
    if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
      for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
	mrc_fld_foreach(f, ix,iy,iz, bnd - 1, bnd - 1) {
	  float b2  = (sqr(.5f * (BX_(f, ix,iy,iz, p) + BX_(f, ix+dx,iy   ,iz   , p))) +
		       sqr(.5f * (BY_(f, ix,iy,iz, p) + BY_(f, ix   ,iy+dy,iz   , p))) +
		       sqr(.5f * (BZ_(f, ix,iy,iz, p) + BZ_(f, ix   ,iy   ,iz+dz, p))));
	  M3(r, 0, ix,iy,iz, p) = EE_(f, ix,iy,iz, p) -.5f * b2 / mhd->par.mu0_code;
	} mrc_fld_foreach_end;
      }
    } else if (MT_BGRID(mhd_type) == MT_BGRID_CC) {
      for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
	mrc_fld_foreach(f, ix,iy,iz, bnd - 1, bnd - 1) {
	  float b2  = (sqr(BX_(f, ix,iy,iz, p)) +
		       sqr(BY_(f, ix,iy,iz, p)) +
		       sqr(BZ_(f, ix,iy,iz, p)));
	  M3(r, 0, ix,iy,iz, p) = EE_(f, ix,iy,iz, p) -.5f * b2 / mhd->par.mu0_code;
	} mrc_fld_foreach_end;
      }
    } else {
      assert(0);
    }
  } else {
    assert(0);
  }

  float max = 0.;
  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
      max = mrc_fld_max(max, mrc_fld_abs(M3(r, 0, ix,iy,iz, p)));
      if (!isfinite(UU_(f, ix,iy,iz, p))) max = 9999.;
    } mrc_fld_foreach_end;
  }
  float max_uu1;
  MPI_Allreduce(&max, &max_uu1, 1, MPI_FLOAT, MPI_MAX, ggcm_mhd_comm(mhd));
  mpi_printf(ggcm_mhd_comm(mhd), "max uu1 = %g\n", max_uu1);

  mrc_fld_put_as(r, fld_r);
  mrc_fld_put_as(f, fld);

  float scale_uu = 1.;
  ggcm_mhd_diag_c_write_one_field(io, fld_r, 0, "uu1", scale_uu, diag_type, plane);
  mrc_fld_destroy(fld_r);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "uu1"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_uu1 = {
  .name             = "uu1",
  .run              = ggcm_mhd_diag_item_uu1_run,
};

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_ee1_run

static void
ggcm_mhd_diag_item_ee1_run(struct ggcm_mhd_diag_item *item,
			   struct mrc_io *io, struct mrc_fld *fld,
			   int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;

  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);

  int bnd = fld->_nr_ghosts;

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "ee1");
  mrc_fld_set_type(fld_r, FLD_TYPE);
  mrc_fld_set_param_int(fld_r, "nr_ghosts", bnd);
  mrc_fld_setup(fld_r);

  struct mrc_fld *r = mrc_fld_get_as(fld_r, FLD_TYPE);
  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);

  if (mhd_type == MT_FCONS_CC) {
    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, bnd - 1, bnd - 1) {
	M3(r, 0, ix,iy,iz, p) = EE_(f, ix,iy,iz, p);
      } mrc_fld_foreach_end;
    }
  } else {
    assert(0);
  }

  mrc_fld_put_as(r, fld_r);
  mrc_fld_put_as(f, fld);

  float scale_ee = 1.;
  ggcm_mhd_diag_c_write_one_field(io, fld_r, 0, "ee1", scale_ee, diag_type, plane);
  mrc_fld_destroy(fld_r);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "ee1"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_ee1 = {
  .name             = "ee1",
  .run              = ggcm_mhd_diag_item_ee1_run,
};

// ======================================================================
// ggcm_mhd_diag_item subclass "rv1"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_rv1_run

static void
ggcm_mhd_diag_item_rv1_run(struct ggcm_mhd_diag_item *item,
			   struct mrc_io *io, struct mrc_fld *fld,
			   int diag_type, float plane)
{
  float scale_rv = 1.;
  ggcm_mhd_diag_c_write_one_field(io, fld, RVX, "rv1x", scale_rv, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, RVY, "rv1y", scale_rv, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, RVZ, "rv1z", scale_rv, diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "rv1"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_rv1 = {
  .name             = "rv1",
  .run              = ggcm_mhd_diag_item_rv1_run,
};

// ======================================================================
// ggcm_mhd_diag_item subclass "b1"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_b1_run

static void
ggcm_mhd_diag_item_b1_run(struct ggcm_mhd_diag_item *item,
			  struct mrc_io *io, struct mrc_fld *fld,
			  int diag_type, float plane)
{
  float scale_bb = 1.;
  ggcm_mhd_diag_c_write_one_field(io, fld, BX, "b1x", scale_bb, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, BY, "b1y", scale_bb, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, BZ, "b1z", scale_bb, diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "b1"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_b1 = {
  .name             = "b1",
  .run              = ggcm_mhd_diag_item_b1_run,
};

