
// ----------------------------------------------------------------------
// constants

static mrc_fld_data_t Gamma _mrc_unused;

// ----------------------------------------------------------------------
// fluxes_fc

static inline void
fluxes_fc(mrc_fld_data_t F[8], mrc_fld_data_t U[8], mrc_fld_data_t W[8])
{
  mrc_fld_data_t b2 = sqr(W[BX]) + sqr(W[BY]) + sqr(W[BZ]);

  F[RR]  = W[RR] * W[VX];
  F[RVX] = W[RR] * W[VX] * W[VX] + W[PP] + .5 * b2 - W[BX] * W[BX];
  F[RVY] = W[RR] * W[VY] * W[VX]                   - W[BY] * W[BX];
  F[RVZ] = W[RR] * W[VZ] * W[VX]                   - W[BZ] * W[BX];
  F[EE] = (U[EE] + W[PP] + .5 * b2) * W[VX]
    - W[BX] * (W[BX] * W[VX] + W[BY] * W[VY] + W[BZ] * W[VZ]);
  F[BX] = 0;
  F[BY] = W[BY] * W[VX] - W[BX] * W[VY];
  F[BZ] = W[BZ] * W[VX] - W[BX] * W[VZ]; 
}

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
// wavespeed_fc
//
// calculate speed of fastest (fast magnetosonic) wave

static inline mrc_fld_data_t
wavespeed_fc(mrc_fld_data_t U[8], mrc_fld_data_t W[8])
{
  mrc_fld_data_t cs2 = Gamma * W[PP] / W[RR];
  mrc_fld_data_t b2 = sqr(W[BX]) + sqr(W[BY]) + sqr(W[BZ]);
  mrc_fld_data_t as2 = b2 / W[RR]; 
  mrc_fld_data_t cf2 = .5f * (cs2 + as2 + 
			      mrc_fld_sqrt(sqr(as2 + cs2) - (4.f * sqr(sqrt(cs2) * W[BX]) / W[RR])));
  return mrc_fld_sqrt(cf2);
}

// ----------------------------------------------------------------------
// fluxes_rusanov_fc

static void
fluxes_rusanov_fc(mrc_fld_data_t F[8], mrc_fld_data_t Ul[8], mrc_fld_data_t Ur[8],
		  mrc_fld_data_t Wl[8], mrc_fld_data_t Wr[8])
{
  mrc_fld_data_t Fl[8], Fr[8];
  mrc_fld_data_t cf;

  cf = wavespeed_fc(Ul, Wl);
  mrc_fld_data_t cp_l = Wl[VX] + cf;
  mrc_fld_data_t cm_l = Wl[VX] - cf; 
  fluxes_fc(Fl, Ul, Wl);

  cf = wavespeed_fc(Ur, Wr);
  mrc_fld_data_t cp_r = Wr[VX] + cf;
  mrc_fld_data_t cm_r = Wr[VX] - cf; 
  fluxes_fc(Fr, Ur, Wr);

  mrc_fld_data_t c_l = mrc_fld_max(mrc_fld_abs(cm_l), mrc_fld_abs(cp_l)); 
  mrc_fld_data_t c_r = mrc_fld_max(mrc_fld_abs(cm_r), mrc_fld_abs(cp_r)); 
  mrc_fld_data_t c_max = mrc_fld_max(c_l, c_r);

  for (int m = 0; m < 8; m++) {
    F[m] = .5f * (Fl[m] + Fr[m] - c_max * (Ur[m] - Ul[m]));
  }
}

// ----------------------------------------------------------------------
// fluxes_rusanov_sc

static void
fluxes_rusanov_sc(mrc_fld_data_t F[5], mrc_fld_data_t Ul[5], mrc_fld_data_t Ur[5],
		  mrc_fld_data_t Wl[5], mrc_fld_data_t Wr[5])
{
  mrc_fld_data_t Fl[5], Fr[5];
  
  fluxes_sc(Fl, Ul, Wl);
  fluxes_sc(Fr, Ur, Wr);

  mrc_fld_data_t vv, cs2;
  vv = sqr(Wl[VX]) + sqr(Wl[VY]) + sqr(Wl[VZ]);
  cs2 = Gamma * Wl[PP] / Wl[RR];
  mrc_fld_data_t cmsv_l = sqrtf(vv) + sqrtf(cs2);

  vv = sqr(Wr[VX]) + sqr(Wr[VY]) + sqr(Wr[VZ]);
  cs2 = Gamma * Wr[PP] / Wr[RR];
  mrc_fld_data_t cmsv_r = sqrtf(vv) + sqrtf(cs2);

  mrc_fld_data_t lambda = .5 * (cmsv_l + cmsv_r);
  
  for (int m = 0; m < 5; m++) {
    F[m] = .5f * ((Fr[m] + Fl[m]) - lambda * (Ur[m] - Ul[m]));
  }
}

// ----------------------------------------------------------------------
// mhd_riemann_rusanov_run_fc

static void _mrc_unused
mhd_riemann_rusanov_run_fc(struct ggcm_mhd *mhd, fld1d_state_t F,
			   fld1d_state_t U_l, fld1d_state_t U_r,
			   fld1d_state_t W_l, fld1d_state_t W_r,
			   int ldim, int l, int r, int dim)
{
  Gamma = mhd->par.gamm;

  for (int i = -l; i < ldim + r; i++) {
    fluxes_rusanov_fc(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
		      &F1S(W_l, 0, i), &F1S(W_r, 0, i));
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
  Gamma = mhd->par.gamm;

  for (int i = -l; i < ldim + r; i++) {
    fluxes_rusanov_sc(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
		      &F1S(W_l, 0, i), &F1S(W_r, 0, i));
  }
}

