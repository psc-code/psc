
#include "mhd_reconstruct_private.h"

#include "ggcm_mhd_defs.h"

#include "mhd_1d.c"

// ======================================================================
// mhd_reconstruct subclass "plm"

// ----------------------------------------------------------------------
// mhd_reconstruct_plm_run_sc

// FIXME bnd
static void
mhd_reconstruct_plm_run_sc(struct mhd_reconstruct *mr,
			   struct mrc_fld *Ul, struct mrc_fld *Ur,
			   struct mrc_fld *Wl, struct mrc_fld *Wr,
			   struct mrc_fld *W1d, struct mrc_fld *Bxi,
			   int ldim, int l, int r, int dir)
{
  mrc_fld_data_t dWc[5], dWl[5];
  mrc_fld_data_t dWr[5], dWg[5];
  mrc_fld_data_t dWm[5];

  for (int i = -1; i < ldim + 1; i++) {
    for (int n = 0; n < 5; n++) {
      dWc[n] = F1(W1d, n, i+1) - F1(W1d, n, i-1);
      dWl[n] = F1(W1d, n, i  ) - F1(W1d, n, i-1);
      dWr[n] = F1(W1d, n, i+1) - F1(W1d, n, i  );
      dWg[n] = (dWl[n] * dWr[n] > 0.) ? 2. * dWl[n] * dWr[n] / (dWl[n] + dWr[n]) : 0.;
    }

    // generalized minmod
    for (int n = 0; n < 5; n++) {
      if (dWl[n] * dWr[n] > 0.0) {
        mrc_fld_data_t min_lr = fmin(    fabs(dWl[n]), fabs(dWr[n]));
        mrc_fld_data_t min_cg = fmin(0.5*fabs(dWc[n]), fabs(dWg[n]));
        dWm[n] = copysign(fmin(2.0 * min_lr, min_cg), dWc[n]);
      } else {
	dWm[n] = 0.0;
      }
    }

    for (int n = 0; n < 5; n++) {
      F1(Wl, n, i+1) = F1(W1d, n, i) + .5 * dWm[n];
      F1(Wr, n, i  ) = F1(W1d, n, i) - .5 * dWm[n];
    }
  }

  mhd_sc_from_prim(mr->mhd, Ul, Wl, ldim, 0, 1);
  mhd_sc_from_prim(mr->mhd, Ur, Wr, ldim, 0, 1);
}

// ----------------------------------------------------------------------
// mhd_reconstruct_plm_run

static void
mhd_reconstruct_plm_run(struct mhd_reconstruct *mr,
			struct mrc_fld *U_l, struct mrc_fld *U_r,
			struct mrc_fld *W_l, struct mrc_fld *W_r,
			struct mrc_fld *W_1d, struct mrc_fld *Bxi,
			int ldim, int l, int r, int dir)
{
  int nr_comps = mrc_fld_dims(U_l)[0];

  switch (nr_comps) {
  case 5: return mhd_reconstruct_plm_run_sc(mr, U_l, U_r, W_l, W_r, W_1d, Bxi, ldim, l, r, dir);
  default: assert(0);
  }
}

