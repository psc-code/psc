
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag_private.h"
#include "ggcm_mhd_step.h"

#include <mrc_domain.h>
#include <mrc_fld_as_float.h>

#include <stdio.h>
#include <assert.h>

// ======================================================================
// ggcm_mhd_diag_item subclass "e_ec"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_e_ec_run

static void
ggcm_mhd_diag_item_e_ec_run(struct ggcm_mhd_diag_item *item,
                            struct mrc_io *io, struct mrc_fld *f,
                            int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;
  float scale_ee = mhd->eenorm;
  struct mrc_fld *E = mrc_domain_fld_create(mhd->domain, SW_2, "ex_ec:ey_ec:ez_ec");
  mrc_fld_set_type(E, FLD_TYPE);
  mrc_fld_setup(E);

  ggcm_mhd_step_get_e_ec(mhd->step, E, f);
  ggcm_mhd_diag_c_write_one_field(io, E, 0, "ex_ec", scale_ee, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, E, 1, "ey_ec", scale_ee, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, E, 2, "ez_ec", scale_ee, diag_type, plane);
  mrc_fld_destroy(E);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "e_ec"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_e_ec = {
  .name             = "e_ec",
  .run              = ggcm_mhd_diag_item_e_ec_run,
};


// ======================================================================
// ggcm_mhd_diag_item subclass "e_cc"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_e_cc_run

static void
ggcm_mhd_diag_item_e_cc_run(struct ggcm_mhd_diag_item *item,
                            struct mrc_io *io, struct mrc_fld *f,
                            int diag_type, float plane)
{
  int gdims[3];
  mrc_domain_get_global_dims(f->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  int mhd_type;
  mrc_fld_get_param_int(f, "mhd_type", &mhd_type);

  struct ggcm_mhd *mhd = item->diag->mhd;
  float scale_ee = mhd->eenorm;
  struct mrc_fld *Eec = mrc_domain_fld_create(mhd->domain, SW_2, "ex_ec:ey_ec:ez_ec");
  struct mrc_fld *Ecc = mrc_domain_fld_create(mhd->domain, SW_2, "ex_cc:ey_cc:ez_cc");
  mrc_fld_set_type(Eec, FLD_TYPE);
  mrc_fld_set_type(Ecc, FLD_TYPE);
  mrc_fld_setup(Eec);
  mrc_fld_setup(Ecc);

  ggcm_mhd_step_get_e_ec(mhd->step, Eec, f);

  // average ec -> cc
  if (MT_BGRID(mhd_type) == MT_BGRID_FC_GGCM) {
    for (int p = 0; p < mrc_fld_nr_patches(Ecc); p++) {
      mrc_fld_foreach(Eec, ix,iy,iz, SW_2 - 1, SW_2) {
	M3(Ecc, 0, ix, iy, iz, p) =
	  0.25f * (M3(Eec, 0, ix  , iy   , iz   , p) + M3(Eec, 0, ix   , iy-dy, iz   , p) +
		   M3(Eec, 0, ix  , iy   , iz-dz, p) + M3(Eec, 0, ix   , iy-dy, iz-dz, p));
	M3(Ecc, 1, ix, iy, iz, p) =
	  0.25f * (M3(Eec, 1, ix  , iy   , iz   , p) + M3(Eec, 1, ix-dx, iy   , iz   , p) +
		   M3(Eec, 1, ix  , iy   , iz-dz, p) + M3(Eec, 1, ix-dx, iy   , iz-dz, p));
	M3(Ecc, 2, ix, iy, iz, p) =
	  0.25f * (M3(Eec, 2, ix  , iy   , iz   , p) + M3(Eec, 2, ix-dx, iy   , iz   , p) +
		   M3(Eec, 2, ix  , iy-dy, iz   , p) + M3(Eec, 2, ix-dx, iy-dy, iz   , p));
      } mrc_fld_foreach_end;
    }
  } else if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
    for (int p = 0; p < mrc_fld_nr_patches(Ecc); p++) {
      mrc_fld_foreach(Eec, ix,iy,iz, SW_2, SW_2 - 1) {
	M3(Ecc, 0, ix, iy, iz, p) =
	  0.25f * (M3(Eec, 0, ix  , iy   , iz   , p) + M3(Eec, 0, ix   , iy+dy, iz   , p) +
		   M3(Eec, 0, ix  , iy   , iz+dz, p) + M3(Eec, 0, ix   , iy+dy, iz+dz, p));
	M3(Ecc, 1, ix, iy, iz, p) =
	  0.25f * (M3(Eec, 1, ix  , iy   , iz   , p) + M3(Eec, 1, ix+dx, iy   , iz   , p) +
		   M3(Eec, 1, ix  , iy   , iz+dz, p) + M3(Eec, 1, ix+dx, iy   , iz+dz, p));
	M3(Ecc, 2, ix, iy, iz, p) =
	  0.25f * (M3(Eec, 2, ix  , iy   , iz   , p) + M3(Eec, 2, ix+dx, iy   , iz   , p) +
		   M3(Eec, 2, ix  , iy+dy, iz   , p) + M3(Eec, 2, ix+dx, iy+dy, iz   , p));
      } mrc_fld_foreach_end;
    }
  } else {
    assert(0);
  }

  ggcm_mhd_diag_c_write_one_field(io, Ecc, 0, "ex_cc", scale_ee, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, Ecc, 1, "ey_cc", scale_ee, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, Ecc, 2, "ez_cc", scale_ee, diag_type, plane);

  mrc_fld_destroy(Eec);
  mrc_fld_destroy(Ecc);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "e_cc"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_e_cc = {
  .name             = "e_cc",
  .run              = ggcm_mhd_diag_item_e_cc_run,
};

