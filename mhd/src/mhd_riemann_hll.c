
#include "mhd_riemann_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_fld_as_double_aos.h>
#define F1(f, m, i) MRC_D2(f, m, i)

#include <math.h>

// ----------------------------------------------------------------------
// fluxes_cc

static void
fluxes_cc(mrc_fld_data_t F[8], mrc_fld_data_t U[8], mrc_fld_data_t W[8], 
	  mrc_fld_data_t gamma, mrc_fld_data_t bb)
{

  mrc_fld_data_t mb = W[BX] * W[VX] + W[BY] * W[VY] + W[BZ] * W[VZ];
  F[RR]  = W[RR] * W[VX];
  F[RVX] = W[RR] * W[VX] * W[VX] + W[PP] + .5 * bb - sqr(W[BX]);
  F[RVY] = W[RR] * W[VY] * W[VX] - W[BX] * W[BY];
  F[RVZ] = W[RR] * W[VZ] * W[VX] - W[BX] * W[BZ];
  F[BY] = W[BY] * W[VX] - W[BX] * W[VY];
  F[BZ] = W[BZ] * W[VX] - W[BX] * W[VZ]; 
  F[UU] = (U[EE] + W[PP] + .5 * bb) * W[VX] - W[BX] * mb ; 
}

// ----------------------------------------------------------------------
// fluxes_hll

static void
fluxes_hll(mrc_fld_data_t F[8], mrc_fld_data_t Ul[8], mrc_fld_data_t Ur[8],
	       mrc_fld_data_t Wl[8], mrc_fld_data_t Wr[8], mrc_fld_data_t gamma)
{
  mrc_fld_data_t Fl[8], Fr[8];
  mrc_fld_data_t bb, cs2, as2, cf;

  bb = sqr(Wl[BX]) + sqr(Wl[BY]) + sqr(Wl[BZ]);
  cs2 = gamma * Wl[PP] / Wl[RR];
  as2 = bb / Wl[RR]; 
  cf = sqrtf(.5 * (cs2 + as2 + sqrtf(sqr(as2 + cs2)
       - (4. * sqr(sqrt(cs2) * Wl[BX]) / Wl[RR]))));       

  fluxes_cc(Fl, Ul, Wl, gamma, bb);
  mrc_fld_data_t cpv_l = Wl[VX] + cf;
  mrc_fld_data_t cmv_l = Wl[VX] - cf; 
  
  bb = sqr(Wr[BX]) + sqr(Wr[BY]) + sqr(Wr[BZ]);
  cs2 = gamma * Wr[PP] / Wr[RR];
  as2 = bb / Wr[RR]; 
  cf = sqrtf(.5 * (cs2 + as2 + sqrtf(sqr(as2 + cs2)
       - (4. * sqr(sqrt(cs2) * Wr[BX]) / Wr[RR]))));     

  fluxes_cc(Fr, Ur, Wr, gamma, bb);
  mrc_fld_data_t cpv_r = Wr[VX] + cf;
  mrc_fld_data_t cmv_r = Wr[VX] - cf; 

  mrc_fld_data_t SR =  fmaxf(fmaxf(cpv_l, cpv_r), 0.); 
  mrc_fld_data_t SL =  fminf(fminf(cmv_l, cmv_r), 0.); 

  for (int m = 0; m <8; m++) {
    F[m] = ((SR * Fl[m] - SL * Fr[m]) + (SR * SL * (Ur[m] - Ul[m]))) / (SR - SL);
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
  mrc_fld_data_t gamma = riemann->mhd->par.gamm;
  for (int i = -l; i < ldim + r; i++) {
    fluxes_hll(&F1(F, 0, i), &F1(U_l, 0, i), &F1(U_r, 0, i),
		   &F1(W_l, 0, i), &F1(W_r, 0, i), gamma);
  }
}

// ----------------------------------------------------------------------
// mhd_riemann_hll_ops

struct mhd_riemann_ops mhd_riemann_hll_ops = {
  .name             = "hll",
  .run              = mhd_riemann_hll_run,
};

