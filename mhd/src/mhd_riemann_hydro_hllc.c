
#include "mhd_riemann_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include <mrc_fld_as_double_aos.h>
#define F1(f, m, i) MRC_D2(f, m, i)

#include <math.h>

// ----------------------------------------------------------------------
// fluxes_cc

static void
fluxes_cc(mrc_fld_data_t F[5], mrc_fld_data_t U[5], mrc_fld_data_t W[5])
{
  F[RR]  = W[RR] * W[VX];
  F[RVX] = W[RR] * W[VX] * W[VX] + W[PP];
  F[RVY] = W[RR] * W[VY] * W[VX];
  F[RVZ] = W[RR] * W[VZ] * W[VX];
  F[UU]  = (U[UU] + W[PP]) * W[VX];
}

// ----------------------------------------------------------------------
// fluxes_hllc

static void
fluxes_hllc(mrc_fld_data_t F[5], mrc_fld_data_t Ul[5], mrc_fld_data_t Ur[5],
	       mrc_fld_data_t Wl[5], mrc_fld_data_t Wr[5], mrc_fld_data_t gamma)
{
  mrc_fld_data_t Fl[5], Fr[5];

  fluxes_cc(Fl, Ul, Wl);
  fluxes_cc(Fr, Ur, Wr);

  mrc_fld_data_t vv, cs2;

  cs2 = gamma * Wl[PP] / Wl[RR];
  mrc_fld_data_t cpv_l = Wl[VX] + sqrtf(cs2);
  mrc_fld_data_t cmv_l = Wl[VX] - sqrtf(cs2); 

  cs2 = gamma * Wr[PP] / Wr[RR];
  mrc_fld_data_t cpv_r = Wr[VX] + sqrtf(cs2);
  mrc_fld_data_t cmv_r = Wr[VX] - sqrtf(cs2); 

  mrc_fld_data_t SR =  fmaxf(fmaxf(cpv_l, cpv_r), 0.); 
  mrc_fld_data_t SL =  fminf(fminf(cmv_l, cmv_r), 0.); 

  mrc_fld_data_t SRmUR = SR - Wr[VX];
  mrc_fld_data_t SLmUL = SL - Wl[VX];
  mrc_fld_data_t SM =
    (SRmUR * Wr[RR] * Wr[VX] - SLmUL * Wl[RR] * Wl[VX] - Wr[PP] + Wl[PP]) / 
    (SRmUR * Wr[RR] - SLmUL * Wl[RR]);

  mrc_fld_data_t spT= Wr[PP] + (Wr[RR] * SRmUR * (SM - Wr[VX]));

  mrc_fld_data_t sUR[5];
  mrc_fld_data_t sUL[5]; 
  
  sUR[0] = Wr[RR] * SRmUR / ( SR - SM );
  sUL[0] = Wl[RR] * SLmUL / ( SL - SM ); 

  sUR[1] = sUR[0] * SM;
  sUL[1] = sUL[0] * SM;
  sUR[2] = sUR[0] * Wr[VY]; 
  sUL[2] = sUL[0] * Wl[VY]; 
  sUR[3] = sUR[0] * Wr[VZ]; 
  sUL[3] = sUL[0] * Wl[VZ];

  sUR[4] = ((SR - Wr[VX]) * Ur[UU] - Wr[PP] * Wr[VX] + spT * SM) / (SR - SM); 
  sUL[4] = ((SL - Wl[VX]) * Ul[UU] - Wl[PP] * Wl[VX] + spT * SM) / (SL - SM); 

 for (int m = 0; m < 5; m++) {
   if ( SL > 0 ) {
     F[m] = Fl[m];
   } else if (( SL <= 0 ) && ( SM >= 0 )) {  
     F[m] = (SL * (sUL[m]-Ul[m])) + Fl[m];
   } else if (( SR >= 0 ) && ( SM <= 0 )) {
     F[m] = (SR * (sUR[m]-Ur[m])) + Fr[m];
   } else if ( SR < 0 ) { 	  
     F[m] = Fr[m];
   }
 }
}

// ----------------------------------------------------------------------
// mhd_riemann_hllc_run

static void
mhd_riemann_hllc_run(struct mhd_riemann *riemann, struct mrc_fld *F,
			struct mrc_fld *U_l, struct mrc_fld *U_r,
			struct mrc_fld *W_l, struct mrc_fld *W_r,
			int ldim, int l, int r, int dim)
{
  mrc_fld_data_t gamma = riemann->mhd->par.gamm;
  for (int i = -l; i < ldim + r; i++) {
    fluxes_hllc(&F1(F, 0, i), &F1(U_l, 0, i), &F1(U_r, 0, i),
		   &F1(W_l, 0, i), &F1(W_r, 0, i), gamma);
  }
}

// ----------------------------------------------------------------------
// mhd_riemann_hydro_hllc_ops

struct mhd_riemann_ops mhd_riemann_hydro_hllc_ops = {
  .name             = "hydro_hllc",
  .run              = mhd_riemann_hllc_run,
};

