
#include "mhd_riemann_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_fld_as_double_aos.h>
#define F1(f, m, i) MRC_D2(f, m, i)

#include <math.h>

// ----------------------------------------------------------------------
// constants

static mrc_fld_data_t Gamma;

// ----------------------------------------------------------------------
// fluxes

static inline void
fluxes(mrc_fld_data_t F[8], mrc_fld_data_t U[8], mrc_fld_data_t W[8])
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
// wavespeed
//
// calculate speed of fastest (fast magnetosonic) wave

static inline mrc_fld_data_t
wavespeed(mrc_fld_data_t U[8], mrc_fld_data_t W[8])
{
  mrc_fld_data_t cs2 = Gamma * W[PP] / W[RR];
  mrc_fld_data_t b2 = sqr(W[BX]) + sqr(W[BY]) + sqr(W[BZ]);
  mrc_fld_data_t as2 = b2 / W[RR]; 
  mrc_fld_data_t cf2 = .5f * (cs2 + as2 + 
			      mrc_fld_sqrt(sqr(as2 + cs2) - (4.f * sqr(sqrt(cs2) * W[BX]) / W[RR])));
  return mrc_fld_sqrt(cf2);
}

// ----------------------------------------------------------------------
// fluxes_hll

static void
fluxes_hll(mrc_fld_data_t F[8], mrc_fld_data_t Ul[8], mrc_fld_data_t Ur[8],
	       mrc_fld_data_t Wl[8], mrc_fld_data_t Wr[8])
{
  mrc_fld_data_t Fl[8], Fr[8];
  mrc_fld_data_t cf;

  cf = wavespeed(Ul, Wl);
  mrc_fld_data_t cp_l = Wl[VX] + cf;
  mrc_fld_data_t cm_l = Wl[VX] - cf; 
  fluxes(Fl, Ul, Wl);

  cf = wavespeed(Ur, Wr);
  mrc_fld_data_t cp_r = Wr[VX] + cf;
  mrc_fld_data_t cm_r = Wr[VX] - cf; 
  fluxes(Fr, Ur, Wr);

  mrc_fld_data_t c_l =  mrc_fld_min(mrc_fld_min(cm_l, cm_r), 0.); 
  mrc_fld_data_t c_r =  mrc_fld_max(mrc_fld_max(cp_l, cp_r), 0.); 

  for (int m = 0; m < 8; m++) {
    F[m] = ((c_r * Fl[m] - c_l * Fr[m]) + (c_r * c_l * (Ur[m] - Ul[m]))) / (c_r - c_l);
  }
}

// ----------------------------------------------------------------------
// mhd_riemann_hll_run

static void
mhd_riemann_hll_run(struct mhd_riemann *riemann, struct mrc_fld *F,
			struct mrc_fld *U_l, struct mrc_fld *U_r,
			struct mrc_fld *W_l, struct mrc_fld *W_r,
			int ldim, int l, int r, int dim)
{
  Gamma = riemann->mhd->par.gamm;

  for (int i = -l; i < ldim + r; i++) {
    fluxes_hll(&F1(F, 0, i), &F1(U_l, 0, i), &F1(U_r, 0, i),
		   &F1(W_l, 0, i), &F1(W_r, 0, i));
  }
}

// ----------------------------------------------------------------------
// mhd_riemann_hll_ops

struct mhd_riemann_ops mhd_riemann_hll_ops = {
  .name             = "hll",
  .run              = mhd_riemann_hll_run,
};

