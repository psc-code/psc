
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
  struct ggcm_mhd *mhd = item->diag->mhd;

  ggcm_mhd_diag_c_write_one_field(io, fld, 10, "ex", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 11, "ey", 1., diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld, 12, "ez", 1., diag_type, plane);

  if (mhd->b0) {
    int bnd = fld->_nr_ghosts - 1;
    
    struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, SW_2, "bx:by:bz");
    mrc_fld_set_type(fld_r, FLD_TYPE);
    mrc_fld_set_param_int(fld_r, "nr_ghosts", bnd + 1);
    mrc_fld_setup(fld_r);
    
    struct mrc_fld *r = mrc_fld_get_as(fld_r, FLD_TYPE);
    struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);
    struct mrc_fld *b0 = mrc_fld_get_as(mhd->b0, FLD_TYPE);
   
    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
				M3(r, 0, ix,iy,iz, p) = M3(f, 13, ix,iy,iz, p) + M3(b0, 0, ix,iy,iz, p);
				M3(r, 1, ix,iy,iz, p) = M3(f, 14, ix,iy,iz, p) + M3(b0, 1, ix,iy,iz, p);
				M3(r, 2, ix,iy,iz, p) = M3(f, 15, ix,iy,iz, p) + M3(b0, 2, ix,iy,iz, p);
      } mrc_fld_foreach_end;
    }

    mrc_fld_put_as(r, fld_r);
    mrc_fld_put_as(f, fld);
    mrc_fld_put_as(b0, mhd->b0);
    
    ggcm_mhd_diag_c_write_one_field(io, fld_r, 0, "bx", 1., diag_type, plane);
    ggcm_mhd_diag_c_write_one_field(io, fld_r, 1, "by", 1., diag_type, plane);
    ggcm_mhd_diag_c_write_one_field(io, fld_r, 2, "bz", 1., diag_type, plane);

    mrc_fld_destroy(fld_r);
  } else {
    ggcm_mhd_diag_c_write_one_field(io, fld, 13, "bx", 1., diag_type, plane);
    ggcm_mhd_diag_c_write_one_field(io, fld, 14, "by", 1., diag_type, plane);
    ggcm_mhd_diag_c_write_one_field(io, fld, 15, "bz", 1., diag_type, plane);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "gkeyll_em"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_gkeyll_em = {
  .name             = "gkeyll_em",
  .run              = ggcm_mhd_diag_item_gkeyll_em_run,
};

