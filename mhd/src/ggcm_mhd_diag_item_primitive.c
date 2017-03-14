
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag_private.h"
#include "ggcm_mhd_gkeyll.h"

#include <mrc_fld_as_double.h>
#include <mrc_domain.h>
#include <stdio.h>
#include <assert.h>

#include "mhd_3d.c"

#include "ggcm_mhd_convert.h"

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

  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);

  int bnd = fld->_nr_ghosts;

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, bnd, "vx:vy:vz");
  mrc_fld_set_type(fld_r, FLD_TYPE);
  mrc_fld_setup(fld_r);

  struct mrc_fld *r = mrc_fld_get_as(fld_r, FLD_TYPE);
  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);

  if (mhd_type == MT_GKEYLL) {
    int nr_fluids = mhd->par.gk_nr_fluids;
    int nr_moments = mhd->par.gk_nr_moments;

    assert(nr_moments == 5);
    int idx[nr_fluids];
    ggcm_mhd_gkeyll_fluid_species_index_all(mhd, idx);

    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, bnd, bnd) {
        M3(r, 0, ix,iy,iz, p) = 0.;
        mrc_fld_data_t rr = 0.;
        for (int s = 0; s < nr_fluids; s++) {
          M3(r, 0, ix,iy,iz, p) += M3(f, idx[s]+G5M_RVXS, ix,iy,iz, p);
          M3(r, 1, ix,iy,iz, p) += M3(f, idx[s]+G5M_RVYS, ix,iy,iz, p);
          M3(r, 2, ix,iy,iz, p) += M3(f, idx[s]+G5M_RVZS, ix,iy,iz, p);
          rr += M3(f, idx[s]+G5M_RRS, ix,iy,iz, p);
        }
        for (int d = 0; d < 3; d++) {
          M3(r, d, ix,iy,iz, p) /= rr;
        }
      } mrc_fld_foreach_end;
    }
  } else {
    for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
      mrc_fld_foreach(f, ix,iy,iz, bnd, bnd) {
        mrc_fld_data_t rri = 1.f / RR_(f, ix,iy,iz, p);
        M3(r, 0, ix,iy,iz, p) = rri * RVX_(f, ix,iy,iz, p);
        M3(r, 1, ix,iy,iz, p) = rri * RVY_(f, ix,iy,iz, p);
        M3(r, 2, ix,iy,iz, p) = rri * RVZ_(f, ix,iy,iz, p);
      } mrc_fld_foreach_end;
    }
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

  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);

  if (mhd_type == MT_GKEYLL) {
    int bnd = fld->_nr_ghosts - 1;

    struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, bnd + 1, "rr");
    mrc_fld_set_type(fld_r, FLD_TYPE);
    mrc_fld_setup(fld_r);

    struct mrc_fld *r = mrc_fld_get_as(fld_r, FLD_TYPE);
    struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);

    int nr_fluids = mhd->par.gk_nr_fluids;
    int nr_moments = mhd->par.gk_nr_moments;

    assert(nr_moments == 5);
    int idx[nr_fluids];
    ggcm_mhd_gkeyll_fluid_species_index_all(mhd, idx);

    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, bnd, bnd) {
        M3(r, 0, ix,iy,iz, p) = 0.;
        for (int s = 0; s < nr_fluids; s++)
          M3(r, 0, ix,iy,iz, p) += M3(f, idx[s]+G5M_RRS, ix,iy,iz, p);
      } mrc_fld_foreach_end;
    }

    mrc_fld_put_as(r, fld_r);
    mrc_fld_put_as(f, fld);

    ggcm_mhd_diag_c_write_one_field(io, r, 0, "rr", scale_rr, diag_type, plane);
  
    mrc_fld_destroy(fld_r);
  } else {
    ggcm_mhd_diag_c_write_one_field(io, fld, RR, "rr", scale_rr, diag_type, plane);
  }
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

  static bool is_setup = false;
  if (!is_setup) {
    ggcm_mhd_convert_setup(mhd);
  }
  
  int gdims[3];
  mrc_domain_get_global_dims(fld->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);

  int bnd = fld->_nr_ghosts - 1;

  struct mrc_fld *fld_r = mrc_domain_fld_create(mhd->domain, bnd + 1, "pp");
  mrc_fld_set_type(fld_r, FLD_TYPE);
  mrc_fld_setup(fld_r);

  struct mrc_fld *r = mrc_fld_get_as(fld_r, FLD_TYPE);
  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);

  if (MT_FORMULATION(mhd_type) == MT_FORMULATION_SCONS) {
    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, bnd, bnd) {
	mrc_fld_data_t prim[5], state[5];
	for (int m = 0; m < 5; m++) {
	  state[m] = M3(f, m, ix,iy,iz, p);
	}
	convert_prim_from_state_scons(prim, state);
	M3(r,0, ix,iy,iz, p) = prim[PP];
      } mrc_fld_foreach_end;
    }
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_FCONS) {
    if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
      for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
	mrc_fld_foreach(f, ix,iy,iz, bnd, bnd) {
	  mrc_fld_data_t prim[5], state[8];
	  for (int m = 0; m < 5; m++) {
	    state[m] = M3(f, m, ix,iy,iz, p);
	  }
	  state[BX] = .5f * (BX_(f, ix,iy,iz, p) + BX_(f, ix+dx,iy   ,iz   , p));
	  state[BY] = .5f * (BY_(f, ix,iy,iz, p) + BY_(f, ix   ,iy+dy,iz   , p));
	  state[BZ] = .5f * (BZ_(f, ix,iy,iz, p) + BZ_(f, ix   ,iy   ,iz+dz, p));
	  convert_prim_from_state_fcons(prim, state);
	  M3(r,0, ix,iy,iz, p) = prim[PP];
	} mrc_fld_foreach_end;
      }
    } else if (MT_BGRID(mhd_type) == MT_BGRID_CC) {
      for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
	mrc_fld_foreach(f, ix,iy,iz, bnd, bnd) {
	  mrc_fld_data_t prim[5], state[8];
	  for (int m = 0; m < 5; m++) {
	    state[m] = M3(f, m, ix,iy,iz, p);
	  }
	  state[BX] = BX_(f, ix,iy,iz, p);
	  state[BY] = BY_(f, ix,iy,iz, p);
	  state[BZ] = BZ_(f, ix,iy,iz, p);
	  convert_prim_from_state_fcons(prim, state);
	  M3(r,0, ix,iy,iz, p) = prim[PP];
	} mrc_fld_foreach_end;
      }
    }
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_GKEYLL) {
    mrc_fld_data_t gamm = mhd->par.gamm;

    int nr_fluids = mhd->par.gk_nr_fluids;
    int nr_moments = mhd->par.gk_nr_moments;

    assert(nr_moments == 5);
    int idx[nr_fluids];
    ggcm_mhd_gkeyll_fluid_species_index_all(mhd, idx);

    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, bnd, bnd) {
        M3(r, 0, ix,iy,iz, p) = 0.;
        for (int s = 0; s < nr_fluids; s++) {
          M3(r, 0, ix,iy,iz, p) += 
            (gamm-1.) * ( M3(f, idx[s]+G5M_UUS, ix,iy,iz, p)
            - .5 * (sqr(M3(f, idx[s]+G5M_RVXS, ix,iy,iz, p))
                  + sqr(M3(f, idx[s]+G5M_RVYS, ix,iy,iz, p))
                  + sqr(M3(f, idx[s]+G5M_RVZS, ix,iy,iz, p)))
                  / M3(f, idx[s]+G5M_RRS, ix,iy,iz, p));
        }
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
  if (MT_BGRID(mhd_type) == MT_BGRID_FC_GGCM) {
    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
	M3(r, 0, ix,iy,iz, p) = .5f * (BX_(f, ix,iy,iz, p) + BX_(f, ix-1,iy,iz, p));
	M3(r, 1, ix,iy,iz, p) = .5f * (BY_(f, ix,iy,iz, p) + BY_(f, ix,iy-1,iz, p));
	M3(r, 2, ix,iy,iz, p) = .5f * (BZ_(f, ix,iy,iz, p) + BZ_(f, ix,iy,iz-1, p));
      } mrc_fld_foreach_end;
    }
  } else if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
    compute_Bt_cc(mhd, fld_r, f, 1, 1);
  } else if (MT_BGRID(mhd_type) == MT_BGRID_CC) {
    struct mrc_fld *b0 = mhd->b0;
    for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
      mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
	M3(r, 0, ix,iy,iz, p) = BT(f, 0, ix,iy,iz, p);
	M3(r, 1, ix,iy,iz, p) = BT(f, 1, ix,iy,iz, p);
	M3(r, 2, ix,iy,iz, p) = BT(f, 2, ix,iy,iz, p);
      } mrc_fld_foreach_end;
    }
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

