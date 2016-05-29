
// ----------------------------------------------------------------------
// mhd_reconstruct_pcm_run_sc

static void _mrc_unused
mhd_reconstruct_pcm_run_sc(struct ggcm_mhd *mhd,
			   fld1d_state_t U_l, fld1d_state_t U_r,
			   fld1d_state_t W_l, fld1d_state_t W_r,
			   fld1d_state_t W, struct mrc_fld *Bxi,
			   int ldim, int l, int r, int dim)
{
  for (int i = -l; i < ldim + r; i++) {
    for (int m = 0; m < 5; m++) {
      F1S(W_l, m, i) = F1S(W, m, i-1);
    }
    for (int m = 0; m < 5; m++) {
      F1S(W_r, m, i) = F1S(W, m, i  );
    }
  }

  mhd_sc_from_prim(mhd, U_l, W_l, ldim, l, r);
  mhd_sc_from_prim(mhd, U_r, W_r, ldim, l, r);
}

// ----------------------------------------------------------------------
// mhd_reconstruct_pcm_run_fc

static void _mrc_unused
mhd_reconstruct_pcm_run_fc(struct ggcm_mhd *mhd,
			   fld1d_state_t U_l, fld1d_state_t U_r,
			   fld1d_state_t W_l, fld1d_state_t W_r,
			   fld1d_state_t W, struct mrc_fld *Bxi,
			   int ldim, int l, int r, int dir)
{
  for (int i = -l; i < ldim + r; i++) {
    for (int m = 0; m < 8; m++) {
      F1S(W_l, m, i) = F1S(W, m, i-1);
    }
    for (int m = 0; m < 8; m++) {
      F1S(W_r, m, i) = F1S(W, m, i  );
    }
  }

  // CHECKME, seems inconsistent to use cell-centered Bx here, then replace it afterwards
  mhd_fc_from_prim(mhd, U_l, W_l, ldim, l, r);
  mhd_fc_from_prim(mhd, U_r, W_r, ldim, l, r);

  if (Bxi) {
    for (int i = -l; i < ldim + r; i++) {
      F1S(W_l, BX, i) = F1(Bxi, 0, i);
      F1S(W_r, BX, i) = F1(Bxi, 0, i);
    }
  }
}

// FIXME too much duplication between sc and fc

// ----------------------------------------------------------------------
// mhd_reconstruct_plm_run_sc

// FIXME bnd
static void _mrc_unused
mhd_reconstruct_plm_run_sc(struct ggcm_mhd *mhd,
			   fld1d_state_t U_l, fld1d_state_t U_r,
			   fld1d_state_t W_l, fld1d_state_t W_r,
			   fld1d_state_t W, struct mrc_fld *Bxi,
			   int ldim, int l, int r, int dir)
{
  mrc_fld_data_t dWc[5], dWl[5];
  mrc_fld_data_t dWr[5], dWg[5];
  mrc_fld_data_t dWm[5];

  for (int i = -1; i < ldim + 1; i++) {
    for (int n = 0; n < 5; n++) {
      dWc[n] = F1S(W, n, i+1) - F1S(W, n, i-1);
      dWl[n] = F1S(W, n, i  ) - F1S(W, n, i-1);
      dWr[n] = F1S(W, n, i+1) - F1S(W, n, i  );
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
      F1S(W_l, n, i+1) = F1S(W, n, i) + .5 * dWm[n];
      F1S(W_r, n, i  ) = F1S(W, n, i) - .5 * dWm[n];
    }
  }

  mhd_sc_from_prim(mhd, U_l, W_l, ldim, 0, 1);
  mhd_sc_from_prim(mhd, U_r, W_r, ldim, 0, 1);
}

// ----------------------------------------------------------------------
// mhd_reconstruct_plm_run_fc

// FIXME bnd
static void _mrc_unused
mhd_reconstruct_plm_run_fc(struct ggcm_mhd *mhd,
			   fld1d_state_t U_l, fld1d_state_t U_r,
			   fld1d_state_t W_l, fld1d_state_t W_r,
			   fld1d_state_t W, struct mrc_fld *Bxi,
			   int ldim, int l, int r, int dir)
{
  mrc_fld_data_t dWc[8], dWl[8];
  mrc_fld_data_t dWr[8], dWg[8];
  mrc_fld_data_t dWm[8];

  // These limits are such that Wl[l..r) and Wr[l..r) both get set,
  // (typically, r = l + 1, so that [l..l] are all set)
  for (int i = -l - 1; i < ldim + r; i++) {
    for (int n = 0; n < 8; n++) {
      dWc[n] = F1S(W, n, i+1) - F1S(W, n, i-1);
      dWl[n] = F1S(W, n, i  ) - F1S(W, n, i-1);
      dWr[n] = F1S(W, n, i+1) - F1S(W, n, i  );
      dWg[n] = (dWl[n] * dWr[n] > 0.) ? 2. * dWl[n] * dWr[n] / (dWl[n] + dWr[n]) : 0.;
    }

    for (int n = 0; n < 8; n++) {
      if (dWl[n] * dWr[n] > 0.0) {
        mrc_fld_data_t lim_slope1 = fmin(    fabs(dWl[n]), fabs(dWr[n]));
        mrc_fld_data_t lim_slope2 = fmin(0.5*fabs(dWc[n]), fabs(dWg[n]));
        dWm[n] = copysign(fmin(2.0 * lim_slope1, lim_slope2), dWc[n]);
      } else {
	dWm[n] = 0.0;
      }
    }

    for (int n = 0; n < 8; n++) {
      F1S(W_l, n, i+1) = F1S(W, n, i) + .5 * dWm[n];
      F1S(W_r, n, i  ) = F1S(W, n, i) - .5 * dWm[n];
    }
    if (Bxi) {
      F1S(W_l, BX, i) = F1(Bxi, 0, i);
      F1S(W_r, BX, i) = F1(Bxi, 0, i);
    }
  }

  mhd_fc_from_prim(mhd, U_l, W_l, ldim, l, r+1);
  mhd_fc_from_prim(mhd, U_r, W_r, ldim, l, r+1);
}

