
#include "mhd_reconstruct_private.h"

#include "ggcm_mhd_defs.h"
#include "mhd_util.h"

#include "pde/pde_setup.c"
#include "pde/pde_mhd_convert.c"

#include "mhd_1d.c"

// ======================================================================
// mhd_reconstruct subclass "pcm"

// ----------------------------------------------------------------------
// mhd_reconstruct_pcm_run_sc

static void
mhd_reconstruct_pcm_run_sc(struct mhd_reconstruct *mr,
			   struct mrc_fld *U_l, struct mrc_fld *U_r,
			   struct mrc_fld *W_l, struct mrc_fld *W_r,
			   struct mrc_fld *W_1d, struct mrc_fld *Bxi,
			   int ldim, int l, int r, int dim)
{
  for (int i = -l; i < ldim + r; i++) {
    for (int m = 0; m < 5; m++) {
      F1(W_l, m, i) = F1(W_1d, m, i-1);
    }
    for (int m = 0; m < 5; m++) {
      F1(W_r, m, i) = F1(W_1d, m, i  );
    }
  }

  mhd_sc_from_prim(mr->mhd, (fld1d_state_t) { .mrc_fld = U_l }, (fld1d_state_t) { .mrc_fld = W_l }, ldim, l, r);
  mhd_sc_from_prim(mr->mhd, (fld1d_state_t) { .mrc_fld = U_r }, (fld1d_state_t) { .mrc_fld = W_r }, ldim, l, r);
}

// ----------------------------------------------------------------------
// mhd_reconstruct_pcm_run_fc

static void
mhd_reconstruct_pcm_run_fc(struct mhd_reconstruct *mr,
			   struct mrc_fld *Ul, struct mrc_fld *Ur,
			   struct mrc_fld *Wl, struct mrc_fld *Wr,
			   struct mrc_fld *W1d, struct mrc_fld *Bxi,
			   int ldim, int l, int r, int dir)
{
  for (int i = -l; i < ldim + r; i++) {
    for (int m = 0; m < 8; m++) {
      F1(Wl, m, i) = F1(W1d, m, i-1);
    }
    for (int m = 0; m < 8; m++) {
      F1(Wr, m, i) = F1(W1d, m, i  );
    }
  }

  // CHECKME, seems inconsistent to use cell-centered Bx here, then replace it afterwards
  mhd_fc_from_prim(mr->mhd, Ul, Wl, ldim, l, r);
  mhd_fc_from_prim(mr->mhd, Ur, Wr, ldim, l, r);

  if (Bxi) {
    for (int i = -l; i < ldim + r; i++) {
      F1(Wl, BX, i) = F1(Bxi, 0, i);
      F1(Wr, BX, i) = F1(Bxi, 0, i);
    }
  }
}

// ----------------------------------------------------------------------
// mhd_reconstruct_pcm_run

static void
mhd_reconstruct_pcm_run(struct mhd_reconstruct *mr,
			struct mrc_fld *U_l, struct mrc_fld *U_r,
			struct mrc_fld *W_l, struct mrc_fld *W_r,
			struct mrc_fld *W_1d, struct mrc_fld *Bxi,
			int ldim, int l, int r, int dir)
{
  int nr_comps = mrc_fld_dims(U_l)[0];
  switch (nr_comps) {
  case 5: return mhd_reconstruct_pcm_run_sc(mr, U_l, U_r, W_l, W_r, W_1d, Bxi, ldim, l, r, dir);
  case 8: return mhd_reconstruct_pcm_run_fc(mr, U_l, U_r, W_l, W_r, W_1d, Bxi, ldim, l, r, dir);
  default: assert(0);
  }
}

