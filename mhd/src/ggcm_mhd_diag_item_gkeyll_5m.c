
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag_private.h"

#include <mrc_fld_as_double.h>
#include <mrc_domain.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

// FIXME, it's kinda lazy to provide just the groups (e, i, em), and there
// really should be a generic way to output state fields just by
// component name in the first place

// ======================================================================
// ggcm_mhd_diag_item subclass "gkeyll_e"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_gkeyll_e_run

static void
ggcm_mhd_diag_item_gkeyll_e_run(struct ggcm_mhd_diag_item *item,
				struct mrc_io *io, struct mrc_fld *fld,
				int diag_type, float plane)
{
  ggcm_mhd_diag_c_write_one_field(io, fld, 0, "rr_e", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 1, "rvx_e", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 2, "rvy_e", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 3, "rvz_e", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 4, "uu_e", 1., diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "gkeyll_e"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_gkeyll_e = {
  .name             = "gkeyll_e",
  .run              = ggcm_mhd_diag_item_gkeyll_e_run,
};

// ======================================================================
// ggcm_mhd_diag_item subclass "gkeyll_i"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_gkeyll_i_run

static void
ggcm_mhd_diag_item_gkeyll_i_run(struct ggcm_mhd_diag_item *item,
				struct mrc_io *io, struct mrc_fld *fld,
				int diag_type, float plane)
{
  ggcm_mhd_diag_c_write_one_field(io, fld, 5, "rr_i", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 6, "rvx_i", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 7, "rvy_i", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 8, "rvz_i", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 9, "uu_i", 1., diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "gkeyll_i"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_gkeyll_i = {
  .name             = "gkeyll_i",
  .run              = ggcm_mhd_diag_item_gkeyll_i_run,
};

// ======================================================================
// ggcm_mhd_diag_item subclass "gkeyll_em"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_gkeyll_em_run

static void
ggcm_mhd_diag_item_gkeyll_em_run(struct ggcm_mhd_diag_item *item,
				 struct mrc_io *io, struct mrc_fld *fld,
				 int diag_type, float plane)
{
  ggcm_mhd_diag_c_write_one_field(io, fld, 10, "ex", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 11, "ey", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 12, "ez", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 13, "bx", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 14, "by", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 15, "bz", 1., diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "gkeyll_em"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_gkeyll_em = {
  .name             = "gkeyll_em",
  .run              = ggcm_mhd_diag_item_gkeyll_em_run,
};

