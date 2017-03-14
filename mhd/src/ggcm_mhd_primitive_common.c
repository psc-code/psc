
#include <mrc_domain.h>

#include "ggcm_mhd_convert.h"

#if MT_FORMULATION(MT) == MT_FORMULATION_SCONS

void
ggcm_mhd_calc_pp_scons(struct ggcm_mhd *mhd, struct mrc_fld *pp, struct mrc_fld *fld)
{
  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      mrc_fld_data_t prim[5], state[5];
      for (int m = 0; m < 5; m++) {
	state[m] = M3(fld, m, ix,iy,iz, p);
      }
      convert_prim_from_state_scons(prim, state);
      M3(pp,0, ix,iy,iz, p) = prim[PP];
    } mrc_fld_foreach_end;
  }
}

#elif MT_FORMULATION(MT) == MT_FORMULATION_FCONS

#if MT_BGRID(MT) == MT_BGRID_FC

void
ggcm_mhd_calc_pp_fcons_fc(struct ggcm_mhd *mhd, struct mrc_fld *pp, struct mrc_fld *fld)
{
  int gdims[3];
  mrc_domain_get_global_dims(fld->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      mrc_fld_data_t prim[5], state[8];
      for (int m = 0; m < 5; m++) {
	state[m] = M3(fld, m, ix,iy,iz, p);
      }
      state[BX] = .5f * (BX_(fld, ix,iy,iz, p) + BX_(fld, ix+dx,iy   ,iz   , p));
      state[BY] = .5f * (BY_(fld, ix,iy,iz, p) + BY_(fld, ix   ,iy+dy,iz   , p));
      state[BZ] = .5f * (BZ_(fld, ix,iy,iz, p) + BZ_(fld, ix   ,iy   ,iz+dz, p));
      convert_prim_from_state_fcons(prim, state);
      M3(pp,0, ix,iy,iz, p) = prim[PP];
    } mrc_fld_foreach_end;
  }
}

#elif MT_BGRID(MT) == MT_BGRID_CC

void
ggcm_mhd_calc_pp_fcons_cc(struct ggcm_mhd *mhd, struct mrc_fld *pp, struct mrc_fld *fld)
{
  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      mrc_fld_data_t prim[5], state[8];
      for (int m = 0; m < 5; m++) {
	state[m] = M3(fld, m, ix,iy,iz, p);
      }
      state[BX] = BX_(fld, ix,iy,iz, p);
      state[BY] = BY_(fld, ix,iy,iz, p);
      state[BZ] = BZ_(fld, ix,iy,iz, p);
      convert_prim_from_state_fcons(prim, state);
      M3(pp,0, ix,iy,iz, p) = prim[PP];
    } mrc_fld_foreach_end;
  }
}

#endif

#endif

