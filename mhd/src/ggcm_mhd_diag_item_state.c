
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag_private.h"

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
  struct ggcm_mhd *mhd = item->diag->mhd;
  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "divB");
  mrc_fld_set_param_int(fld_r, "nr_comps", _NR_FLDS);
  mrc_fld_setup(fld_r);

  struct mrc_fld *f = mrc_fld_get_as(fld, mrc_fld_type(fld));
  struct mrc_fld *r = mrc_fld_get_as(fld_r, "float");

  mrc_fld_foreach(r, ix,iy,iz, 0, 0) {
    MRC_F3(r, _RR1, ix,iy,iz) = MRC_F3(f, _RR1, ix,iy,iz);
  } mrc_fld_foreach_end;

  mrc_fld_put_as(f, fld);
  mrc_fld_put_as(r, fld_r);

  float scale_rr = 1.;
  ggcm_mhd_diag_c_write_one_field(io, fld_r, _RR1, "rr1", scale_rr, diag_type, plane);

  mrc_fld_destroy(fld_r);
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
  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "uu1");
  mrc_fld_set_param_int(fld_r, "nr_comps", _NR_FLDS);
  mrc_fld_setup(fld_r);
  //  mrc_fld_copy(fld_r, fld);

  struct mrc_fld *f = mrc_fld_get_as(fld, mrc_fld_type(fld));
  struct mrc_fld *r = mrc_fld_get_as(fld_r, "float");

  float max = 0.;
  mrc_fld_foreach(r, ix,iy,iz, 0, 0) {
    MRC_F3(r,_UU1, ix,iy,iz) = MRC_F3(f,_UU1, ix,iy,iz);
    max = fmaxf(max, fabsf(MRC_F3(r,_UU1, ix,iy,iz)));
    if (!isfinite(MRC_F3(r,_UU1, ix,iy,iz))) max = 9999.;
  } mrc_fld_foreach_end;

  mrc_fld_put_as(f, fld);
  mrc_fld_put_as(r, fld_r);

  float scale_uu = 1.;
  ggcm_mhd_diag_c_write_one_field(io, fld_r, _UU1, "uu1", scale_uu, diag_type, plane);

  mrc_fld_destroy(fld_r);
  mprintf("max uu1 = %g\n", max);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "uu1"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_uu1 = {
  .name             = "uu1",
  .run              = ggcm_mhd_diag_item_uu1_run,
};

// ======================================================================
// ggcm_mhd_diag_item subclass "rv1"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_rv1_run

static void
ggcm_mhd_diag_item_rv1_run(struct ggcm_mhd_diag_item *item,
			   struct mrc_io *io, struct mrc_fld *f,
			   int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;
  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "divB");
  mrc_fld_set_param_int(fld_r, "nr_comps", mrc_fld_nr_comps(f));
  mrc_fld_setup(fld_r);
  mrc_fld_copy(fld_r, f);

  float scale_rv = 1.;
  ggcm_mhd_diag_c_write_one_field(io, fld_r, _RV1X, "rv1x", scale_rv, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld_r, _RV1Y, "rv1y", scale_rv, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld_r, _RV1Z, "rv1z", scale_rv, diag_type, plane);

  mrc_fld_destroy(fld_r);
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
			  struct mrc_io *io, struct mrc_fld *f,
			  int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;
  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "divB");
  mrc_fld_set_param_int(fld_r, "nr_comps", mrc_fld_nr_comps(f));
  mrc_fld_setup(fld_r);
  mrc_fld_copy(fld_r, f);

  float scale_bb = 1.;
  ggcm_mhd_diag_c_write_one_field(io, fld_r, _B1X, "b1x", scale_bb, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld_r, _B1Y, "b1y", scale_bb, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld_r, _B1Z, "b1z", scale_bb, diag_type, plane);

  mrc_fld_destroy(fld_r);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "b1"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_b1 = {
  .name             = "b1",
  .run              = ggcm_mhd_diag_item_b1_run,
};

