
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
limit_minmod(mrc_fld_data_t dWm, mrc_fld_data_t dWp)
{
  return minmod(dWm, dWp);
}

// ----------------------------------------------------------------------
// limit_mc
//
// monotonized central limiter

static inline mrc_fld_data_t
limit_mc(mrc_fld_data_t dWm, mrc_fld_data_t dWp)
{
  mrc_fld_data_t dWc = .5f * (dWm + dWp);
  return minmod(dWc, s_limiter_mc_beta * minmod(dWm, dWp));
}

// ----------------------------------------------------------------------
// limit_gminmod
//
// generalized minmod (athena)

static inline mrc_fld_data_t
limit_gminmod(mrc_fld_data_t dWm, mrc_fld_data_t dWp)
{
  mrc_fld_data_t dWc = .5f * (dWm + dWp);
  mrc_fld_data_t dWg = 2. * dWm * dWp / (dWm + dWp);
  return minmod(minmod(dWc, dWg), 2.f * minmod(dWm, dWp));
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
  } else if (s_opt_limiter == OPT_LIMITER_GMINMOD) {
    for (int m = 0; m < s_n_comps; m++) {
      dW[m] = limit_gminmod(dWm[m], dWp[m]);
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
			 fld1d_state_t W, fld1d_t bx, int ib, int ie)
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
    if (F1(s_aux.bnd_mask, i) == 2.f) {
      // force constant reconstruction next to boundary
      for (int m = 0; m < s_n_comps; m++) {
	dW[m] = 0.f;
      }
    } else {
      limit_slope(dW, dWm, dWp);
    }

    // l/r states based on limited slope
    for (int m = 0; m < s_n_comps; m++) {
      F1S(Wl, m, i+1) = F1S(W, m, i) + .5f * dW[m];
      F1S(Wr, m, i  ) = F1S(W, m, i) - .5f * dW[m];
    }

    if (fld1d_is_setup(bx)) {
      F1S(Wl, BX, i) = F1(bx, i);
      F1S(Wr, BX, i) = F1(bx, i);
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
  } else {
    mhd_reconstruct_plm_prim(U_l, U_r, W_l, W_r, W, bx, ib, ie);
  }
}
