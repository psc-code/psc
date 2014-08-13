
#include "mhd_riemann_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_fld_as_double_aos.h>
#define F1(f, m, i) MRC_D2(f, m, i)
#include <math.h>
#define sign(x) (( x > 0. ) - ( x < 0. ))


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
// calc_state_s

static void
calc_state_s(mrc_fld_data_t Ws[8], mrc_fld_data_t U[8], mrc_fld_data_t W[8], 
	     mrc_fld_data_t S, mrc_fld_data_t SM, mrc_fld_data_t sPt, 
	     mrc_fld_data_t SmU)
{
  // MK eq 43 
  Ws[RR] = W[RR] * SmU / (S - SM);

  // MK eq 39
  Ws[VX] = SM;

  mrc_fld_data_t cden  = 1. / ((W[RR] * SmU * (S-SM))-sqr(W[BX]));

  // MK eq 44 & eq 46
  Ws[VY] = W[VY] - W[BX] * W[BY] * (SM - W[VX]) * cden;
  Ws[VZ] = W[VZ] - W[BX] * W[BZ] * (SM - W[VX]) * cden;
  
  // MK eq 45 & eq 47 
  Ws[BX] = U[BX];   
  Ws[BY] = W[BY] * (W[RR] * sqr(SmU) - sqr(W[BX])) * cden;
  Ws[BZ] = W[BZ] * (W[RR] * sqr(SmU) - sqr(W[BX])) * cden; 
    
  // MK eq 48
  mrc_fld_data_t vb = W[VX] * W[BX] + W[VY] * W[BY] + W[VZ] * W[BZ]; 
  mrc_fld_data_t vsbs = Ws[VX] * Ws[BX] + Ws[VY] * Ws[BY] + Ws[VZ] * Ws[BZ]; 
  Ws[EE] = (SmU * U[EE] - W[PP] * W[VX] + sPt * SM + W[BX]* ( vb -vsbs )) / (S - SM);
}


// ----------------------------------------------------------------------
// fluxes_hlld
// 
// Miyoshi & Kusano (2005)

static void
fluxes_hlld(mrc_fld_data_t F[8], mrc_fld_data_t Ul[8], mrc_fld_data_t Ur[8],
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
    mrc_fld_data_t SRmUR = SR - Wr[VX];
    mrc_fld_data_t SLmUL = SL - Wl[VX];
    
    // MK eq. 38
    mrc_fld_data_t SM =
      (SRmUR * Wr[RR] * Wr[VX] - SLmUL * Wl[RR] * Wl[VX] - Wr[PP] + Wl[PP]) / 
      (SRmUR * Wr[RR] - SLmUL * Wl[RR]);
    
    // MK eq. 41
    mrc_fld_data_t sPt= (SRmUR * Wr[RR] * Wl[PP] - SLmUL * Wl[RR] * Wr[PP] + 
			 Wl[RR] * Wr[RR] * SRmUR * SLmUL * (Wr[VX] - Wl[VX])) / 
      (SRmUR * Wr[RR] - SLmUL * Wl[RR]);
    
    mrc_fld_data_t Urs[8], Uls[8], Wls[8], Wrs[8]; 
    mrc_fld_data_t Urss[8], Ulss[8], Wlss[8], Wrss[8];
    
    calc_state_s(Wls, Ul, Wl, SL, SM, sPt, SLmUL);
    calc_state_s(Wrs, Ur, Wr, SR, SM, sPt, SRmUR);
    
    // MK eq. 49 
    Wlss[RR] = Wls[RR];
    Wrss[RR] = Wrs[RR];
    
    // MK eq. 50 
    mrc_fld_data_t ssPt = sPt; 
    
    // MK eq. 59-63
    mrc_fld_data_t cden = 1./(sqrt(Wls[RR]) + sqrt(Wrs[RR]));
    
    Wlss[VY] = (sqrt(Wls[RR]) * Wls[VY] + sqrt(Wrs[RR]) * Wrs[VY] 
		+ (Wrs[BY] - Wls[BY]) * sign(Wl[BX])) * cden; 
    Wlss[VZ] = (sqrt(Wls[RR]) * Wls[VZ] + sqrt(Wrs[RR]) * Wrs[VZ] 
		+ (Wrs[BZ] - Wls[BZ]) * sign(Wl[BX])) * cden; 
    Wlss[BY] = (sqrt(Wls[RR]) * Wrs[BY] + sqrt(Wrs[RR]) * Wls[BY] 
		+ sqrt( Wls[RR] * Wrs[RR] ) *(Wrs[VY] - Wls[VY])
		* sign(Wl[BX])) * cden; 
    Wlss[BZ] = (sqrt(Wls[RR]) * Wrs[BZ] + sqrt(Wrs[RR]) * Wls[BZ] 
		+ sqrt( Wls[RR] * Wrs[RR] ) *(Wrs[VZ] - Wls[VZ])
		* sign(Wl[BX])) * cden; 
    
    Wrss[VY] = Wlss[VY];
    Wrss[VZ] = Wlss[VZ];    
    Wrss[BY] = Wlss[BY];
    Wrss[BZ] = Wlss[BZ]; 

    mrc_fld_data_t vbs = 
      Wls[VX] * Wls[BX] + Wls[VY] * Wls[BY] + Wls[VZ] * Wls[BZ]; 
    mrc_fld_data_t vbss = 
      Wlss[VX] * Wlss[BX] + Wlss[VY] * Wlss[BY] + Wlss[VZ] * Wlss[BZ];   
    Wlss[EE] = Wls[EE] - sqrt(Wls[RR]) * ( vbs - vbss ) * sign(Wls[BX]);   
    vbs  = Wrs[VX] * Wrs[BX] + Wrs[VY] * Wrs[BY] + Wrs[VZ] * Wrs[BZ]; 
     Wrss[EE] = Wrs[EE] + sqrt(Wrs[RR]) * ( vbs - vbss ) * sign(Wrs[BX]);    
    
    Urs[RR] = Wrs[RR];
    Uls[RR] = Wls[RR];
    Urs[RVX] = Wrs[RR] * SM;
    Uls[RVX] = Wls[RR] * SM;
    Urs[RVY] = Wrs[RR] * Wrs[VY]; 
    Uls[RVY] = Wls[RR] * Wls[VY]; 
    Urs[RVZ] = Wrs[RR] * Wrs[VZ]; 
    Uls[RVZ] = Wls[RR] * Wls[VZ];
    
    Urs[BX] = Wrs[BX];
    Uls[BX] = Wls[BX]; 
    Urs[BY] = Wrs[BY];
    Uls[BY] = Wls[BY];  
    Urs[BZ] = Wrs[BZ];
    Uls[BZ] = Wls[BZ];   
    Uls[EE] = Wls[EE]; 
    Urs[EE] = Wrs[EE];
    
    Urss[RR] = Wrss[RR];
    Ulss[RR] = Wlss[RR];
    Urss[RVX] = Wrss[RR] * SM;
    Ulss[RVX] = Wlss[RR] * SM;
    Urss[RVY] = Wrss[RR] * Wrss[VY]; 
    Ulss[RVY] = Wlss[RR] * Wlss[VY]; 
    Urss[RVZ] = Wrss[RR] * Wrss[VZ]; 
    Ulss[RVZ] = Wlss[RR] * Wlss[VZ];
    
    Urss[BX] = Wrs[BX];
    Ulss[BX] = Wls[BX]; 
    Urss[BY] = Wrss[BY];
    Ulss[BY] = Wlss[BY];  
    Urss[BZ] = Wrss[BZ];
    Ulss[BZ] = Wlss[BZ];   
    Ulss[EE] = Wlss[EE]; 
    Urss[EE] = Wrss[EE];
    
    // MK eq. 51
    mrc_fld_data_t SLs = SM - fabs(Wl[BX]) / sqrt(Wls[RR]) ; 
    mrc_fld_data_t SRs = SM + fabs(Wr[BX]) / sqrt(Wrs[RR]) ;
        
    for (int m = 0; m < 8; m++) {
      if ( SL > 0 ) {
	F[m] = Fl[m];
      } else if (( SL <= 0 ) && ( SLs >= 0 )) {  
	F[m] = Fl[m] + (SL * (Uls[m] - Ul[m]));
      } else if (( SLs <= 0 ) && ( SM >= 0 )) { 
	F[m] = Fl[m] + SLs * Ulss[m] - (SLs - SL) * Uls[m] - SL * Ul[m];
      } else if (( SRs >= 0 ) && ( SM <= 0 )) {
	F[m] = Fr[m] + SRs * Urss[m] - (SRs - SR) * Urs[m] - SR * Ur[m];
      } else if (( SRs <= 0 ) && ( SR >= 0 )) {
	F[m] = Fr[m] + (SR * (Urs[m] - Ur[m]));
      } else if ( SR < 0 ) { 	  
	F[m] = Fr[m];
      } else {
	assert(0);
      }
    }		      
}

// ----------------------------------------------------------------------
// mhd_riemann_hlld_run

static void
mhd_riemann_hlld_run(struct mhd_riemann *riemann, struct mrc_fld *F,
			struct mrc_fld *U_l, struct mrc_fld *U_r,
			struct mrc_fld *W_l, struct mrc_fld *W_r,
			int ldim, int l, int r, int dim)
{
  mrc_fld_data_t gamma = riemann->mhd->par.gamm;
  for (int i = -l; i < ldim + r; i++) {
    fluxes_hlld(&F1(F, 0, i), &F1(U_l, 0, i), &F1(U_r, 0, i),
		   &F1(W_l, 0, i), &F1(W_r, 0, i), gamma);
  }
}

// ----------------------------------------------------------------------
// mhd_riemann_hlld_ops

struct mhd_riemann_ops mhd_riemann_hlld_ops = {
  .name             = "hlld",
  .run              = mhd_riemann_hlld_run,
};

