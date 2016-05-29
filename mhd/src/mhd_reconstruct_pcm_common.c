
#include "mhd_reconstruct_private.h"

#include "ggcm_mhd_defs.h"
#include "mhd_util.h"

#include "pde/pde_setup.c"
#include "pde/pde_mhd_convert.c"
#include "pde/pde_mhd_reconstruct.c"

// ======================================================================
// mhd_reconstruct subclass "pcm"

// ----------------------------------------------------------------------
// mhd_reconstruct_pcm_run

static void
mhd_reconstruct_pcm_run(struct mhd_reconstruct *mr,
			struct mrc_fld *U_l, struct mrc_fld *U_r,
			struct mrc_fld *W_l, struct mrc_fld *W_r,
			struct mrc_fld *W, struct mrc_fld *Bxi,
			int ldim, int l, int r, int dir)
{
  int nr_comps = mrc_fld_dims(U_l)[0];
  switch (nr_comps) {
  case 5: return mhd_reconstruct_pcm_run_sc(mr->mhd, (fld1d_state_t) { .mrc_fld = U_l }, (fld1d_state_t) { .mrc_fld = U_r}, (fld1d_state_t) { .mrc_fld = W_l }, (fld1d_state_t) { .mrc_fld = W_r }, (fld1d_state_t) { .mrc_fld = W }, Bxi, ldim, l, r, dir);
  case 8: return mhd_reconstruct_pcm_run_fc(mr->mhd, (fld1d_state_t) { .mrc_fld = U_l }, (fld1d_state_t) { .mrc_fld = U_r}, (fld1d_state_t) { .mrc_fld = W_l }, (fld1d_state_t) { .mrc_fld = W_r }, (fld1d_state_t) { .mrc_fld = W }, Bxi, ldim, l, r, dir);
  default: assert(0);
  }
}

