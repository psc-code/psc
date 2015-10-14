
#include "mhd_riemann_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <math.h>

// ----------------------------------------------------------------------
// fluxes_cc

static void // FIXME, duplicated
fluxes_cc(mrc_fld_data_t F[5], mrc_fld_data_t U[5], mrc_fld_data_t W[5])
{
  F[RR]  = W[RR] * W[VX];
  F[RVX] = W[RR] * W[VX] * W[VX];
  F[RVY] = W[RR] * W[VY] * W[VX];
  F[RVZ] = W[RR] * W[VZ] * W[VX];
  F[UU]  = (U[UU] + W[PP]) * W[VX];
}

// ----------------------------------------------------------------------
// fluxes_rusanov

static void
fluxes_rusanov(mrc_fld_data_t F[5], mrc_fld_data_t Ul[5], mrc_fld_data_t Ur[5],
	       mrc_fld_data_t Wl[5], mrc_fld_data_t Wr[5], mrc_fld_data_t gamma)
{
  mrc_fld_data_t Fl[5], Fr[5];
  
  fluxes_cc(Fl, Ul, Wl);
  fluxes_cc(Fr, Ur, Wr);

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
// mhd_riemann_rusanov_run

static void
mhd_riemann_rusanov_run(struct mhd_riemann *riemann, struct mrc_fld *F,
			struct mrc_fld *U_l, struct mrc_fld *U_r,
			struct mrc_fld *W_l, struct mrc_fld *W_r,
			int ldim, int l, int r, int dim)
{
  mrc_fld_data_t gamma = riemann->mhd->par.gamm;
  for (int i = -l; i < ldim + r; i++) {
    fluxes_rusanov(&F1(F, 0, i), &F1(U_l, 0, i), &F1(U_r, 0, i),
		   &F1(W_l, 0, i), &F1(W_r, 0, i), gamma);
  }
}

