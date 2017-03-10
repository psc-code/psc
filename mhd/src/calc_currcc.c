
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_gkeyll.h"

#include <mrc_domain.h>
#include <mrc_fld_as_double.h>

#include <math.h>
#include <string.h>

// ----------------------------------------------------------------------
// ggcm_mhd_calc_curcc_fc_ggcm

static void
ggcm_mhd_calc_currcc_fc_ggcm(struct ggcm_mhd *mhd, struct mrc_fld *fld, int m,
			     struct mrc_fld *currcc)
{
  int gdims[3];
  mrc_domain_get_global_dims(fld->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  struct mrc_fld *tmp = mrc_fld_duplicate(currcc);

  struct mrc_fld *f = ggcm_mhd_fld_get_as(fld, FLD_TYPE, MT_SEMI_CONSERVATIVE_GGCM,
					  BX, BX + 3);
  struct mrc_fld *t = mrc_fld_get_as(tmp, FLD_TYPE);
  struct mrc_fld *c = mrc_fld_get_as(currcc, FLD_TYPE);
  
  for (int p = 0; p < mrc_fld_nr_patches(t); p++) {
    float *bd4x = ggcm_mhd_crds_get_crd_p(mhd->crds, 0, BD4, p);
    float *bd4y = ggcm_mhd_crds_get_crd_p(mhd->crds, 1, BD4, p);
    float *bd4z = ggcm_mhd_crds_get_crd_p(mhd->crds, 2, BD4, p);

    mrc_fld_foreach(t,ix,iy,iz, 1, 1) {
      // compute current on edge first
      M3(t,0,ix,iy,iz, p) =
	(M3(f, m+2, ix,iy+dy,iz, p) - M3(f, m+2, ix,iy,iz, p)) * bd4y[iy] -
	(M3(f, m+1, ix,iy,iz+dz, p) - M3(f, m+1, ix,iy,iz, p)) * bd4z[iz];
      M3(t,1,ix,iy,iz, p) =
	(M3(f, m  , ix,iy,iz+dz, p) - M3(f, m  , ix,iy,iz, p)) * bd4z[iz] -
	(M3(f, m+2, ix+dx,iy,iz, p) - M3(f, m+2, ix,iy,iz, p)) * bd4x[ix];
      M3(t,2, ix,iy,iz, p) =
	(M3(f,m+1 , ix+dx,iy,iz, p) - M3(f, m+1, ix,iy,iz, p)) * bd4x[ix] -
	(M3(f,m   , ix,iy+dy,iz, p) - M3(f, m  , ix,iy,iz, p)) * bd4y[iy];
    } mrc_fld_foreach_end;
  }

  for (int p = 0; p < mrc_fld_nr_patches(c); p++) {
    mrc_fld_foreach(c, ix,iy,iz, 1, 1) {
      // average to the center
      M3(c, 0, ix,iy,iz, p) = 0.25f *
	(M3(t, 0, ix,iy,iz   , p) + M3(t, 0, ix,iy-dy,iz   , p) +
	 M3(t, 0, ix,iy,iz-dz, p) + M3(t, 0, ix,iy-dy,iz-dz, p));
      M3(c, 1, ix,iy,iz, p) = 0.25f *
	(M3(t, 1, ix,iy,iz   , p) + M3(t, 1, ix-dx,iy,iz   , p) +
	 M3(t, 1, ix,iy,iz-dz, p) + M3(t, 1, ix-dx,iy,iz-dz, p));
      M3(c, 2, ix,iy,iz, p) = 0.25f *
	(M3(t, 2, ix,iy   ,iz, p) + M3(t, 2, ix-dx,iy   ,iz, p) +
	 M3(t, 2, ix,iy-dy,iz, p) + M3(t, 2, ix-dx,iy-dy,iz, p));
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(c, currcc);
  mrc_fld_put_as(t, tmp);
  ggcm_mhd_fld_put_as(f, fld, 0, 0);
  mrc_fld_destroy(tmp); 
}

// ----------------------------------------------------------------------
// ggcm_mhd_calc_curcc_cc

static void
ggcm_mhd_calc_currcc_cc(struct ggcm_mhd *mhd, struct mrc_fld *fld, int m,
			struct mrc_fld *currcc)
{
  int gdims[3];
  mrc_domain_get_global_dims(fld->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);
  struct mrc_fld *c = mrc_fld_get_as(currcc, FLD_TYPE);
  
  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    float *fd1x = ggcm_mhd_crds_get_crd_p(mhd->crds, 0, FD1, p);
    float *fd1y = ggcm_mhd_crds_get_crd_p(mhd->crds, 1, FD1, p);
    float *fd1z = ggcm_mhd_crds_get_crd_p(mhd->crds, 2, FD1, p);

    mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
      M3(c, 0, ix,iy,iz, p) =
	(M3(f, m+2, ix,iy+dy,iz, p) - M3(f, m+2, ix,iy-dy,iz, p)) * .5f*fd1y[iy] -
	(M3(f, m+1, ix,iy,iz+dz, p) - M3(f, m+1, ix,iy,iz-dz, p)) * .5f*fd1z[iz];
      M3(c, 1, ix,iy,iz, p) =
	(M3(f, m+0, ix,iy,iz+dz, p) - M3(f, m+0, ix,iy,iz-dz, p)) * .5f*fd1z[iz] -
	(M3(f, m+2, ix+dx,iy,iz, p) - M3(f, m+2, ix-dx,iy,iz, p)) * .5f*fd1x[ix];
      M3(c, 2, ix,iy,iz, p) =
	(M3(f, m+1, ix+dx,iy,iz, p) - M3(f, m+1, ix-dx,iy,iz, p)) * .5f*fd1x[ix] -
	(M3(f, m+0, ix,iy+dy,iz, p) - M3(f, m+0, ix,iy-dy,iz, p)) * .5f*fd1y[iy];
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(c, currcc);
  mrc_fld_put_as(f, fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_calc_curcc_gkeyll

static void
ggcm_mhd_calc_currcc_gkeyll(struct ggcm_mhd *mhd, struct mrc_fld *fld, int m,
			     struct mrc_fld *currcc)
{
  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);
  struct mrc_fld *c = mrc_fld_get_as(currcc, FLD_TYPE);

  int nr_fluids = mhd->par.gk_nr_fluids;
  int nr_moments = mhd->par.gk_nr_moments;
  float *mass = mhd->par.gk_mass.vals;
  float *charge = mhd->par.gk_charge.vals;

  float q_m[nr_fluids];
  for (int s = 0; s < nr_fluids; s++)
    q_m[s] = charge[s] / mass[s];

  assert(nr_moments == 5);
  int idx[nr_fluids];
  ggcm_mhd_gkeyll_fluid_species_index_all(mhd, idx);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
      for (int d = 0; d < 3; d++)
        M3(c, d, ix,iy,iz, p) = 0.;
      for (int s = 0; s < nr_fluids; s++) {
        M3(c, 0, ix,iy,iz, p) += M3(f, idx[s]+G5M_RVXS, ix,iy,iz, p) * q_m[s];
        M3(c, 1, ix,iy,iz, p) += M3(f, idx[s]+G5M_RVYS, ix,iy,iz, p) * q_m[s];
        M3(c, 2, ix,iy,iz, p) += M3(f, idx[s]+G5M_RVZS, ix,iy,iz, p) * q_m[s];
      }
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(c, currcc);
  mrc_fld_put_as(f, fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_calc_currcc

void
ggcm_mhd_calc_currcc(struct ggcm_mhd *mhd, struct mrc_fld *fld, int m,
		     struct mrc_fld *currcc)
{
  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);

  if (MT_FORMULATION(mhd_type) == MT_FORMULATION_GKEYLL) {
    return ggcm_mhd_calc_currcc_gkeyll(mhd, fld, m, currcc);
  } else { // MHD
    switch (MT_BGRID(mhd_type)) {
    case MT_BGRID_CC:
      return ggcm_mhd_calc_currcc_cc(mhd, fld, m, currcc);
    case MT_BGRID_FC_GGCM:
      return ggcm_mhd_calc_currcc_fc_ggcm(mhd, fld, m, currcc);
    }
  }

  assert(0);
}
