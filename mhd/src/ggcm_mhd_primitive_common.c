
#include <mrc_domain.h>
#include <mrc_common.h>

#include "pde/pde_defs.h"
#include "pde/pde_mhd_convert.c"

void
SFX(ggcm_mhd_calc_pp)(struct ggcm_mhd *mhd, struct mrc_fld *pp, struct mrc_fld *fld)
{
  static bool is_setup = false;
  if (!is_setup) {
    pde_mhd_setup(mhd, mrc_fld_nr_comps(mhd->fld));
  }
  
  int gdims[3];
  mrc_domain_get_global_dims(fld->_domain, gdims);
  int dx _mrc_unused = (gdims[0] > 1), dy _mrc_unused = (gdims[1] > 1), dz _mrc_unused = (gdims[2] > 1);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, i,j,k, 0, 0) {
      mrc_fld_data_t state[s_n_state], prim[8];
#if MT_FORMULATION(MT) == MT_FORMULATION_GKEYLL
      convert_get_state_from_3d(state, fld, i,j,k, p);
#else
      for (int m = 0; m < 5; m++) {
	state[m] = M3(fld, m, i,j,k, p);
      }
#if MT_BGRID(MT) == MT_BGRID_FC
      state[BX] = .5f * (BX_(fld, i,j,k, p) + BX_(fld, i+dx,j   ,k   , p));
      state[BY] = .5f * (BY_(fld, i,j,k, p) + BY_(fld, i   ,j+dy,k   , p));
      state[BZ] = .5f * (BZ_(fld, i,j,k, p) + BZ_(fld, i   ,j   ,k+dz, p));
#elif MT_BGRID(MT) == MT_BGRID_CC
      state[BX] = BX_(fld, i,j,k, p);
      state[BY] = BY_(fld, i,j,k, p);
      state[BZ] = BZ_(fld, i,j,k, p);
#endif
#endif
      convert_prim_from_state(prim, state);
      M3(pp, 0, i,j,k, p) = prim[PP];
    } mrc_fld_foreach_end;
  }

  if (0) { // FIXME, this one just kills the warning
    pde_free();
  }
}
