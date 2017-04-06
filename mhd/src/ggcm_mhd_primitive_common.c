
#include <mrc_domain.h>
#include <mrc_common.h>

#include "pde/pde_defs.h"
#include "pde/pde_mhd_convert.c"

static bool is_setup = false;

void
SFX(ggcm_mhd_calc_rr)(struct ggcm_mhd *mhd, struct mrc_fld *pp, struct mrc_fld *fld)
{
  if (!is_setup) {
    pde_mhd_setup(mhd, mrc_fld_nr_comps(mhd->fld));
  }
  
  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, i,j,k, 0, 0) {
      mrc_fld_data_t state[s_n_state], prim[8];
      convert_get_cc_state_from_3d(state, fld, i,j,k, p);
      convert_prim_from_state(prim, state);
      M3(pp, 0, i,j,k, p) = prim[RR];
    } mrc_fld_foreach_end;
  }

  if (0) { // FIXME, this one just kills the warning
    pde_free();
  }
}

void
SFX(ggcm_mhd_calc_v)(struct ggcm_mhd *mhd, struct mrc_fld *v, struct mrc_fld *fld)
{
  if (!is_setup) {
    pde_mhd_setup(mhd, mrc_fld_nr_comps(mhd->fld));
  }
  
  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, i,j,k, 0, 0) {
      mrc_fld_data_t state[s_n_state], prim[8];
      convert_get_cc_state_from_3d(state, fld, i,j,k, p);
      convert_prim_from_state(prim, state);
      M3(v, 0, i,j,k, p) = prim[VX];
      M3(v, 1, i,j,k, p) = prim[VY];
      M3(v, 2, i,j,k, p) = prim[VZ];
    } mrc_fld_foreach_end;
  }

  if (0) { // FIXME, this one just kills the warning
    pde_free();
  }
}

void
SFX(ggcm_mhd_calc_pp)(struct ggcm_mhd *mhd, struct mrc_fld *pp, struct mrc_fld *fld)
{
  if (!is_setup) {
    pde_mhd_setup(mhd, mrc_fld_nr_comps(mhd->fld));
  }
  
  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, i,j,k, 0, 0) {
      mrc_fld_data_t state[s_n_state], prim[8];
      convert_get_cc_state_from_3d(state, fld, i,j,k, p);
      convert_prim_from_state(prim, state);
      M3(pp, 0, i,j,k, p) = prim[PP];
    } mrc_fld_foreach_end;
  }

  if (0) { // FIXME, this one just kills the warning
    pde_free();
  }
}

