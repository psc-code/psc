
// ----------------------------------------------------------------------
// mhd_reconstruct_pcm

static void _mrc_unused
mhd_reconstruct_pcm(fld1d_state_t U_l, fld1d_state_t U_r,
		    fld1d_state_t W_l, fld1d_state_t W_r,
		    fld1d_state_t W, struct mrc_fld *Bxi,
		    int ldim, int l, int r, int dir)
{
  for (int i = -l; i < ldim + r; i++) {
    for (int m = 0; m < s_n_comps; m++) {
      F1S(W_l, m, i) = F1S(W, m, i-1);
    }
    for (int m = 0; m < s_n_comps; m++) {
      F1S(W_r, m, i) = F1S(W, m, i  );
    }
  }

  // CHECKME, seems inconsistent to use cell-centered Bx here, then replace it afterwards
  mhd_cons_from_prim(U_l, W_l, ldim, l, r);
  mhd_cons_from_prim(U_r, W_r, ldim, l, r);

  // if not doing fully conservative, Bxi will be NULL, so the following will be skipped
  if (Bxi) {
    for (int i = -l; i < ldim + r; i++) {
      F1S(W_l, BX, i) = F1(Bxi, 0, i);
      F1S(W_r, BX, i) = F1(Bxi, 0, i);
    }
  }
}

// ----------------------------------------------------------------------
// mhd_reconstruct_plm

// FIXME bnd
static void _mrc_unused
mhd_reconstruct_plm(fld1d_state_t U_l, fld1d_state_t U_r,
		    fld1d_state_t W_l, fld1d_state_t W_r,
		    fld1d_state_t W, struct mrc_fld *Bxi,
		    int ldim, int l, int r, int dir)
{
  mrc_fld_data_t dWc[s_n_comps], dWl[s_n_comps];
  mrc_fld_data_t dWr[s_n_comps], dWg[s_n_comps];
  mrc_fld_data_t dWm[s_n_comps];

  // These limits are such that Wl[l..r) and Wr[l..r) both get set,
  // (typically, r = l + 1, so that [l..l] are all set)
  for (int i = -l - 1; i < ldim + r; i++) {
    for (int n = 0; n < s_n_comps; n++) {
      dWc[n] = F1S(W, n, i+1) - F1S(W, n, i-1);
      dWl[n] = F1S(W, n, i  ) - F1S(W, n, i-1);
      dWr[n] = F1S(W, n, i+1) - F1S(W, n, i  );
      dWg[n] = (dWl[n] * dWr[n] > 0.) ? 2. * dWl[n] * dWr[n] / (dWl[n] + dWr[n]) : 0.;
    }

    for (int n = 0; n < s_n_comps; n++) {
      if (dWl[n] * dWr[n] > 0.0) {
        mrc_fld_data_t min_lr = fmin(    fabs(dWl[n]), fabs(dWr[n]));
        mrc_fld_data_t min_cg = fmin(0.5*fabs(dWc[n]), fabs(dWg[n]));
        dWm[n] = copysign(fmin(2.f * min_lr, min_cg), dWc[n]);
      } else {
	dWm[n] = 0.0;
      }
    }

    for (int n = 0; n < s_n_comps; n++) {
      F1S(W_l, n, i+1) = F1S(W, n, i) + .5 * dWm[n];
      F1S(W_r, n, i  ) = F1S(W, n, i) - .5 * dWm[n];
    }

    if (Bxi) {
      F1S(W_l, BX, i) = F1(Bxi, 0, i);
      F1S(W_r, BX, i) = F1(Bxi, 0, i);
    }
  }

  mhd_cons_from_prim(U_l, W_l, ldim, l, r+1);
  mhd_cons_from_prim(U_r, W_r, ldim, l, r+1);
}

// ----------------------------------------------------------------------
// mhd_reconstruct

static void _mrc_unused
mhd_reconstruct(fld1d_state_t U_l, fld1d_state_t U_r,
		fld1d_state_t W_l, fld1d_state_t W_r,
		fld1d_state_t W, struct mrc_fld *Bxi,
		int ldim, int l, int r, int dir)
{
  if (s_opt_limiter == OPT_LIMITER_FLAT) {
    mhd_reconstruct_pcm(U_l, U_r, W_l, W_r, W, Bxi, ldim, l, r, dir);
  } else if (s_opt_limiter == OPT_LIMITER_GMINMOD) {
    mhd_reconstruct_plm(U_l, U_r, W_l, W_r, W, Bxi, ldim, l, r, dir);
  } else {
    assert(0);
  }
}
