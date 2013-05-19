
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag_private.h"

#include <stdio.h>
#include <assert.h>

// ======================================================================
// ggcm_mhd_diag_item subclass "rr1"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_rr1_run

static void
ggcm_mhd_diag_item_rr1_run(struct ggcm_mhd_diag_item *item,
			   struct mrc_io *io, struct mrc_fld *f,
			   int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;

  float scale_rr = mhd->par.rrnorm;
  ggcm_mhd_diag_c_write_one_field(io, f, _RR1, "rr1", scale_rr, diag_type, plane);
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
			   struct mrc_io *io, struct mrc_fld *f,
			   int diag_type, float plane)
{
  float scale_uu = 1.;
  ggcm_mhd_diag_c_write_one_field(io, f, _UU1, "uu1", scale_uu, diag_type, plane);
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
  float scale_rv = 1.;
  ggcm_mhd_diag_c_write_one_field(io, f, _RV1X, "rv1x", scale_rv, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, f, _RV1Y, "rv1y", scale_rv, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, f, _RV1Z, "rv1z", scale_rv, diag_type, plane);
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
  float scale_rv = 1.;
  ggcm_mhd_diag_c_write_one_field(io, f, _B1X, "b1x", scale_rv, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, f, _B1Y, "b1y", scale_rv, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, f, _B1Z, "b1z", scale_rv, diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "b1"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_b1 = {
  .name             = "b1",
  .run              = ggcm_mhd_diag_item_b1_run,
};

