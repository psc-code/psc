
// ----------------------------------------------------------------------
// mhd_reconstruct_pcm

static void _mrc_unused
mhd_reconstruct_pcm(fld1d_state_t U_l, fld1d_state_t U_r,
		    fld1d_state_t W_l, fld1d_state_t W_r,
		    fld1d_state_t W, fld1d_t bx,
		    int ib, int ie)
{
  for (int i = ib; i < ie; i++) {
    for (int m = 0; m < s_n_comps; m++) {
      F1S(W_l, m, i) = F1S(W, m, i-1);
    }
    for (int m = 0; m < s_n_comps; m++) {
      F1S(W_r, m, i) = F1S(W, m, i  );
    }
  }

  // CHECKME, seems inconsistent to use cell-centered Bx here, then replace it afterwards
  mhd_cons_from_prim(U_l, W_l, ib, ie);
  mhd_cons_from_prim(U_r, W_r, ib, ie);

  // if not doing fully conservative, bx will be NULL, so the following will be skipped
  if (fld1d_is_setup(bx)) {
    for (int i = ib; i < ie; i++) {
      F1S(W_l, BX, i) = F1(bx, i);
      F1S(W_r, BX, i) = F1(bx, i);
    }
  }
}

// ----------------------------------------------------------------------
// mhd_reconstruct_plm_gminmod

static void _mrc_unused
mhd_reconstruct_plm_gminmod(fld1d_state_t U_l, fld1d_state_t U_r,
			    fld1d_state_t W_l, fld1d_state_t W_r,
			    fld1d_state_t W, fld1d_t bx, int ib, int ie)
{
  mrc_fld_data_t dWc[s_n_comps], dWl[s_n_comps];
  mrc_fld_data_t dWr[s_n_comps], dWg[s_n_comps];
  mrc_fld_data_t dWm[s_n_comps];

  // These limits are such that Wl[ib..ie) and Wr[ib..ie) both get set,
  for (int i = ib - 1; i < ie; i++) {
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

    if (fld1d_is_setup(bx)) {
      F1S(W_l, BX, i) = F1(bx, i);
      F1S(W_r, BX, i) = F1(bx, i);
    }
  }

  mhd_cons_from_prim(U_l, W_l, ib, ie);
  mhd_cons_from_prim(U_r, W_r, ib, ie);
}

// ----------------------------------------------------------------------
// minmod

static inline mrc_fld_data_t
minmod(mrc_fld_data_t a, mrc_fld_data_t b)
{
  if (a * b > 0.) {
    return mrc_fld_abs(a) < mrc_fld_abs(b) ? a : b;
  } else {
    return 0.;
  }
}

// ----------------------------------------------------------------------
// limit_minmod
//
// minmod limiter

static inline mrc_fld_data_t
limit_minmod(mrc_fld_data_t dvm, mrc_fld_data_t dvp)
{
  return minmod(dvm, dvp);
}

// ----------------------------------------------------------------------
// limit_mc
//
// monotonized central limiter

static inline mrc_fld_data_t
limit_mc(mrc_fld_data_t dvm, mrc_fld_data_t dvp)
{
  mrc_fld_data_t dvc = .5 * (dvm + dvp);
  return minmod(dvc, s_limiter_mc_beta * minmod(dvm, dvp));
}

// ----------------------------------------------------------------------
// limit_slope

static void
limit_slope(mrc_fld_data_t dW[], mrc_fld_data_t dWm[], mrc_fld_data_t dWp[])
{
  if (s_opt_limiter == OPT_LIMITER_MINMOD) {
    for (int m = 0; m < s_n_comps; m++) {
      dW[m] = limit_minmod(dWm[m], dWp[m]);
    }
  } else if (s_opt_limiter == OPT_LIMITER_MC) {
    for (int m = 0; m < s_n_comps; m++) {
      dW[m] = limit_mc(dWm[m], dWp[m]);
    }
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// reconstruct_plm_prim
//
// piecewise linear slope-limited reconstruction on primitive variables
//
// out: reconstructed prim face states Wl, Wr
//      reconstructed cons face states Ul, Ur
// in:  prim variable vector W
//      flag
//      [ib, ie[ interval where the face-centered Wl, Wr, Ul, Ur are calculated, ie., to get
//      proper l/r states including end faces, this should be [0,mx+1[
//      FIXME: SHIFT is a hack to determine how faces are indexed

static void
mhd_reconstruct_plm_prim(fld1d_state_t Ul, fld1d_state_t Ur,
			 fld1d_state_t Wl, fld1d_state_t Wr,
			 fld1d_state_t W, int ib, int ie)
{
  for (int i = ib - 1; i < ie; i++) {
    // one-sided differences after geometric correction
    mrc_fld_data_t dWm[s_n_comps], dWp[s_n_comps];
    for (int m = 0; m < s_n_comps; m++) {
      dWm[m] = F1S(W, m, i  ) - F1S(W, m, i-1);
      dWp[m] = F1S(W, m, i+1) - F1S(W, m, i  );
    }

    // find limited slope
    mrc_fld_data_t dW[s_n_comps];
    limit_slope(dW, dWm, dWp);

    // l/r states based on limited slope
    for (int m = 0; m < s_n_comps; m++) {
      F1S(Wl, m, i+1) = F1S(W, m, i) + .5f * dW[m];
      F1S(Wr, m, i  ) = F1S(W, m, i) - .5f * dW[m];
    }
  }

  // set conservative states, too
  mhd_cons_from_prim(Ul, Wl, ib, ie);
  mhd_cons_from_prim(Ur, Wr, ib, ie);
}

// ----------------------------------------------------------------------
// mhd_reconstruct

static void _mrc_unused
mhd_reconstruct(fld1d_state_t U_l, fld1d_state_t U_r,
		fld1d_state_t W_l, fld1d_state_t W_r,
		fld1d_state_t W, fld1d_t bx, int ib, int ie)
{
  if (s_opt_limiter == OPT_LIMITER_FLAT) {
    mhd_reconstruct_pcm(U_l, U_r, W_l, W_r, W, bx, ib, ie);
  } else if (s_opt_limiter == OPT_LIMITER_GMINMOD) {
    mhd_reconstruct_plm_gminmod(U_l, U_r, W_l, W_r, W, bx, ib, ie);
  } else {
    mhd_reconstruct_plm_prim(U_l, U_r, W_l, W_r, W, ib, ie);
  }
}
