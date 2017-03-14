
#include "ggcm_mhd_private.h"

#include <mrc_fld_as_double.h>

void ggcm_mhd_calc_pp_gkeyll(struct ggcm_mhd *mhd, struct mrc_fld *pp_base,
			     struct mrc_fld *fld_base);

// ----------------------------------------------------------------------
// ggcm_mhd_calc_pp

void
ggcm_mhd_calc_pp(struct ggcm_mhd *mhd, struct mrc_fld *pp_base,
		 struct mrc_fld *fld_base)
{
  int mhd_type;
  mrc_fld_get_param_int(fld_base, "mhd_type", &mhd_type);

  struct mrc_fld *pp = mrc_fld_get_as(pp_base, FLD_TYPE);
  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);

  if (MT_FORMULATION(mhd_type) == MT_FORMULATION_SCONS) {
    ggcm_mhd_calc_pp_scons(mhd, pp, fld);
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_FCONS) {
    if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
      ggcm_mhd_calc_pp_fcons_fc(mhd, pp, fld);
    } else if (MT_BGRID(mhd_type) == MT_BGRID_CC) {
      ggcm_mhd_calc_pp_fcons_cc(mhd, pp, fld);
    }
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_GKEYLL) {
    ggcm_mhd_calc_pp_gkeyll(mhd, pp, fld);
  } else {
    assert(0);
  }

  mrc_fld_put_as(pp, pp_base);
  mrc_fld_put_as(fld, fld_base);
}


#include "ggcm_mhd_gkeyll.h"

// ----------------------------------------------------------------------
// ggcm_mhd_calc_pp_gkeyll

void
ggcm_mhd_calc_pp_gkeyll(struct ggcm_mhd *mhd, struct mrc_fld *pp,
			struct mrc_fld *fld)
{
  mrc_fld_data_t gamm = mhd->par.gamm;

  int nr_fluids = mhd->par.gk_nr_fluids;
  int nr_moments = mhd->par.gk_nr_moments;
  
  assert(nr_moments == 5);
  int idx[nr_fluids];
  ggcm_mhd_gkeyll_fluid_species_index_all(mhd, idx);
  
  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      M3(pp, 0, ix,iy,iz, p) = 0.;
      for (int s = 0; s < nr_fluids; s++) {
	M3(pp, 0, ix,iy,iz, p) += 
	  (gamm-1.) * ( M3(fld, idx[s]+G5M_UUS, ix,iy,iz, p)
			- .5 * (sqr(M3(fld, idx[s]+G5M_RVXS, ix,iy,iz, p)) +
				sqr(M3(fld, idx[s]+G5M_RVYS, ix,iy,iz, p)) +
				sqr(M3(fld, idx[s]+G5M_RVZS, ix,iy,iz, p)))
			/ M3(fld, idx[s]+G5M_RRS, ix,iy,iz, p));
      }
    } mrc_fld_foreach_end;
  }
}
