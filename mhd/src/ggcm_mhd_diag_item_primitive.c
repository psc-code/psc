
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag_private.h"

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
			  struct mrc_io *io, struct mrc_fld *f,
			  int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;

  primvar1_c(mhd);
  float scale_pp = mhd->par.ppnorm;
  ggcm_mhd_diag_c_write_one_field(io, f, _PP, "pp", scale_pp, diag_type, plane);
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
			 struct mrc_io *io, struct mrc_fld *f,
			 int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;

  assert(f == mhd->fld);
  primbb_c(mhd);
  float scale_bb = mhd->par.bbnorm;
  ggcm_mhd_diag_c_write_one_field(io, f, _BX, "bx", scale_bb, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, f, _BY, "by", scale_bb, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, f, _BZ, "bz", scale_bb, diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "b"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_b = {
  .name             = "b",
  .run              = ggcm_mhd_diag_item_b_run,
};


