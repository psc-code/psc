
// ----------------------------------------------------------------------
// fluxes_sc

static void // FIXME, duplicated
fluxes_sc(mrc_fld_data_t F[5], mrc_fld_data_t U[5], mrc_fld_data_t W[5])
{
  F[RR]  = W[RR] * W[VX];
  F[RVX] = W[RR] * W[VX] * W[VX];
  F[RVY] = W[RR] * W[VY] * W[VX];
  F[RVZ] = W[RR] * W[VZ] * W[VX];
  F[UU]  = (U[UU] + W[PP]) * W[VX];
}

// ----------------------------------------------------------------------
// fluxes_rusanov_sc

static void
fluxes_rusanov_sc(mrc_fld_data_t F[5], mrc_fld_data_t Ul[5], mrc_fld_data_t Ur[5],
		  mrc_fld_data_t Wl[5], mrc_fld_data_t Wr[5], mrc_fld_data_t gamma)
{
  mrc_fld_data_t Fl[5], Fr[5];
  
  fluxes_sc(Fl, Ul, Wl);
  fluxes_sc(Fr, Ur, Wr);

  mrc_fld_data_t vv, cs2;
  vv = sqr(Wl[VX]) + sqr(Wl[VY]) + sqr(Wl[VZ]);
  cs2 = gamma * Wl[PP] / Wl[RR];
  mrc_fld_data_t cmsv_l = sqrtf(vv) + sqrtf(cs2);

  vv = sqr(Wr[VX]) + sqr(Wr[VY]) + sqr(Wr[VZ]);
  cs2 = gamma * Wr[PP] / Wr[RR];
  mrc_fld_data_t cmsv_r = sqrtf(vv) + sqrtf(cs2);

  mrc_fld_data_t lambda = .5 * (cmsv_l + cmsv_r);
  
  for (int m = 0; m < 5; m++) {
    F[m] = .5f * ((Fr[m] + Fl[m]) - lambda * (Ur[m] - Ul[m]));
  }
}

// ----------------------------------------------------------------------
// mhd_riemann_rusanov_run_sc

static void _mrc_unused
mhd_riemann_rusanov_run_sc(struct ggcm_mhd *mhd, fld1d_state_t F,
			   fld1d_state_t U_l, fld1d_state_t U_r,
			   fld1d_state_t W_l, fld1d_state_t W_r,
			   int ldim, int l, int r, int dim)
{
  mrc_fld_data_t gamma = mhd->par.gamm;
  for (int i = -l; i < ldim + r; i++) {
    fluxes_rusanov_sc(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
		      &F1S(W_l, 0, i), &F1S(W_r, 0, i), gamma);
  }
}

