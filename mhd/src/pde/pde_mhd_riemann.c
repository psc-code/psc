
#ifndef PDE_MHD_RIEMANN_C
#define PDE_MHD_RIEMANN_C

#include "pde/pde_mhd_divb_glm.c"

// FIXME, at least uppercase
#define sign(x) (( x > 0. ) - ( x < 0. ))

// ----------------------------------------------------------------------
// fluxes_mhd_fcons

static inline void
fluxes_mhd_fcons(mrc_fld_data_t F[], mrc_fld_data_t U[], mrc_fld_data_t W[], int i)
{
  mrc_fld_data_t B0X, B0Y, B0Z;
  if (s_opt_background) {
    mrc_fld_data_t *B0 = &F1V(s_aux.b0, 0, i);
    B0X = B0[0]; B0Y = B0[1]; B0Z = B0[2];
  } else {
    B0X = B0Y = B0Z = 0.f;
  }
  mrc_fld_data_t BTX = B0X + W[BX], BTY = B0Y + W[BY], BTZ = B0Z + W[BZ];
  mrc_fld_data_t b2 = sqr(W[BX]) + sqr(W[BY]) + sqr(W[BZ]);
  mrc_fld_data_t ptot = W[PP] + s_mu0_inv * (.5f * b2 + B0X*W[BX] + B0Y*W[BY] + B0Z*W[BZ]);
  mrc_fld_data_t v_dot_B = (W[BX] * W[VX] + W[BY] * W[VY] + W[BZ] * W[VZ]);

  F[RR]  = W[RR] * W[VX];
  F[RVX] = W[RR] * W[VX] * W[VX] + ptot - s_mu0_inv * (BTX * W[BX] + W[BX] * B0X);
  F[RVY] = W[RR] * W[VY] * W[VX]        - s_mu0_inv * (BTX * W[BY] + W[BX] * B0Y);
  F[RVZ] = W[RR] * W[VZ] * W[VX]        - s_mu0_inv * (BTX * W[BZ] + W[BX] * B0Z);
  F[EE] = (U[EE] + ptot) * W[VX]        - s_mu0_inv * BTX * v_dot_B;
  F[BX] = 0;
  F[BY]  = W[VX] * BTY - W[VY] * BTX;
  F[BZ]  = W[VX] * BTZ - W[VZ] * BTX; 

  if (s_opt_hall == OPT_HALL_CONST) {
    mrc_fld_data_t *j = &F1V(s_aux.j, 0, i);
    F[BY] -= s_d_i * (j[0] * BTY - j[1] * BTX);
    F[BZ] -= s_d_i * (j[0] * BTZ - j[2] * BTX);
    // FIXME, energy contribution
  } else if (s_opt_hall == OPT_HALL_YES) {
    mrc_fld_data_t *j = &F1V(s_aux.j, 0, i);
    F[BY] -= s_d_i / W[RR] * (j[0] * BTY - j[1] * BTX);
    F[BZ] -= s_d_i / W[RR] * (j[0] * BTZ - j[2] * BTX);
    // FIXME, energy contribution
  }

  if (s_opt_divb == OPT_DIVB_GLM) {
    F[BX ] = W[PSI];
    F[PSI] = sqr(s_divb_glm_ch) * W[BX];
  }
}

// ----------------------------------------------------------------------
// fluxes_mhd_scons

static void
fluxes_mhd_scons(mrc_fld_data_t F[5], mrc_fld_data_t U[5], mrc_fld_data_t W[5], int i)
{
  F[RR]  = W[RR] * W[VX];
  F[RVX] = W[RR] * W[VX] * W[VX];
  F[RVY] = W[RR] * W[VY] * W[VX];
  F[RVZ] = W[RR] * W[VZ] * W[VX];
  F[UU]  = (U[UU] + W[PP]) * W[VX];
}

// ----------------------------------------------------------------------
// fluxes_hd

// FIXME, mhd_scons vs hd is a kinda arbitrary distinction, but the former does not
// contain pressure, because that's what's needed for Jimmy-MHD ("c3")

static void
fluxes_hd(mrc_fld_data_t F[5], mrc_fld_data_t U[5], mrc_fld_data_t W[5], int i)
{
  F[RR]  = W[RR] * W[VX];
  F[RVX] = W[RR] * W[VX] * W[VX] + W[PP];
  F[RVY] = W[RR] * W[VY] * W[VX];
  F[RVZ] = W[RR] * W[VZ] * W[VX];
  F[UU]  = (U[UU] + W[PP]) * W[VX];
}

// ----------------------------------------------------------------------
// fluxes

static void
fluxes(mrc_fld_data_t F[], mrc_fld_data_t U[], mrc_fld_data_t W[], int i)
{
  if (s_opt_eqn == OPT_EQN_MHD_FCONS) {
    fluxes_mhd_fcons(F, U, W, i);
  } else if (s_opt_eqn == OPT_EQN_MHD_SCONS) {
    fluxes_mhd_scons(F, U, W, i);
  } else if (s_opt_eqn == OPT_EQN_HD) {
    fluxes_hd(F, U, W, i);
  }
}

// ----------------------------------------------------------------------
// wavespeed_mhd_fcons
//
// calculate speed of fastest (fast magnetosonic) wave

static inline mrc_fld_data_t
wavespeed_mhd_fcons(mrc_fld_data_t U[], mrc_fld_data_t W[], int i)
{
  // FIXME, replicated from fluxes() above
  mrc_fld_data_t B0X, B0Y, B0Z;
  if (s_opt_background) {
    mrc_fld_data_t *B0 = &F1V(s_aux.b0, 0, i);
    B0X = B0[0]; B0Y = B0[1]; B0Z = B0[2];
  } else {
    B0X = B0Y = B0Z = 0.f;
  }
  mrc_fld_data_t BTX = B0X + W[BX], BTY = B0Y + W[BY], BTZ = B0Z + W[BZ];

  // OPT: 1/rr can be factored out, and inner square root can be written in terms
  // of By^2 + Bz^2
  mrc_fld_data_t cs2 = s_gamma * W[PP] / W[RR];
  mrc_fld_data_t bt2 = sqr(BTX) + sqr(BTY) + sqr(BTZ);
  mrc_fld_data_t vA2 = bt2 / W[RR] * s_mu0_inv; 
  mrc_fld_data_t cf2 = .5f * (cs2 + vA2 + 
			      mrc_fld_sqrt(sqr(vA2 + cs2) - (4.f * cs2 * s_mu0_inv * sqr(BTX) / W[RR])));
  mrc_fld_data_t cf = mrc_fld_sqrt(cf2);

  if (s_opt_hall == OPT_HALL_CONST) {
    mrc_fld_data_t cw = s_d_i * mrc_fld_sqrt(bt2) * s_mu0_inv * M_PI * PDE_INV_DS(i);
    cf += cw;
  } else if (s_opt_hall == OPT_HALL_YES) {
    mrc_fld_data_t cw = s_d_i / W[RR] * mrc_fld_sqrt(bt2) * s_mu0_inv * M_PI * PDE_INV_DS(i);
    cf += cw;
  }

  return cf;
}

// ----------------------------------------------------------------------
// wavespeed_mhd_scons
//
// calculate speed of fastest wave (soundspeed)

static inline mrc_fld_data_t
wavespeed_mhd_scons(mrc_fld_data_t U[], mrc_fld_data_t W[], int i)
{
  return sqrtf(s_gamma * W[PP] / W[RR]);
}

// ----------------------------------------------------------------------
// wavespeed

static inline mrc_fld_data_t
wavespeed(mrc_fld_data_t U[], mrc_fld_data_t W[], int i)
{
  if (s_opt_eqn == OPT_EQN_MHD_FCONS) {
    return wavespeed_mhd_fcons(U, W, i);
  } else if (s_opt_eqn == OPT_EQN_MHD_SCONS ||
	     s_opt_eqn == OPT_EQN_HD) {
    return wavespeed_mhd_scons(U, W, i);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// fluxes_rusanov
//
// FIXME? scons/hydro do things weirdly, IIRC to match what original OpenGGCM
// is doing.

static void
fluxes_rusanov(mrc_fld_data_t F[], mrc_fld_data_t Ul[], mrc_fld_data_t Ur[],
	       mrc_fld_data_t Wl[], mrc_fld_data_t Wr[], int i)
{
  mrc_fld_data_t Fl[s_n_comps], Fr[s_n_comps];
  mrc_fld_data_t cf, c_l, c_r, c_max;

  cf = wavespeed(Ul, Wl, i);
  if (s_opt_eqn == OPT_EQN_MHD_FCONS) {
    mrc_fld_data_t cp_l = Wl[VX] + cf;
    mrc_fld_data_t cm_l = Wl[VX] - cf; 
    c_l = mrc_fld_max(mrc_fld_abs(cm_l), mrc_fld_abs(cp_l)); 
  } else if (s_opt_eqn == OPT_EQN_MHD_SCONS ||
	     s_opt_eqn == OPT_EQN_HD) {
    mrc_fld_data_t vv = sqr(Wl[VX]) + sqr(Wl[VY]) + sqr(Wl[VZ]);
    c_l = sqrtf(vv) + cf;
  }

  cf = wavespeed(Ur, Wr, i);
  if (s_opt_eqn == OPT_EQN_MHD_FCONS) {
    mrc_fld_data_t cp_r = Wr[VX] + cf;
    mrc_fld_data_t cm_r = Wr[VX] - cf; 
    c_r = mrc_fld_max(mrc_fld_abs(cm_r), mrc_fld_abs(cp_r)); 
  } else if (s_opt_eqn == OPT_EQN_MHD_SCONS ||
	     s_opt_eqn == OPT_EQN_HD) {
    mrc_fld_data_t vv = sqr(Wr[VX]) + sqr(Wr[VY]) + sqr(Wr[VZ]);
    c_r = sqrtf(vv) + cf;
  }

  if (s_opt_eqn == OPT_EQN_MHD_FCONS) {
    c_max = mrc_fld_max(c_l, c_r);
  } else if (s_opt_eqn == OPT_EQN_MHD_SCONS ||
	     s_opt_eqn == OPT_EQN_HD) {
    c_max = .5 * (c_l + c_r);
  }

  fluxes(Fl, Ul, Wl, i);
  fluxes(Fr, Ur, Wr, i);

  for (int m = 0; m < s_n_comps; m++) {
    F[m] = .5f * (Fl[m] + Fr[m] - c_max * (Ur[m] - Ul[m]));
  }
}

// ----------------------------------------------------------------------
// fluxes_hll

static void
fluxes_hll(mrc_fld_data_t F[], mrc_fld_data_t Ul[], mrc_fld_data_t Ur[],
	   mrc_fld_data_t Wl[], mrc_fld_data_t Wr[], int i)
{
  mrc_fld_data_t Fl[s_n_comps], Fr[s_n_comps];
  mrc_fld_data_t cf;

  cf = wavespeed(Ul, Wl, i);
  mrc_fld_data_t cp_l = Wl[VX] + cf;
  mrc_fld_data_t cm_l = Wl[VX] - cf;

  cf = wavespeed(Ur, Wr, i);
  mrc_fld_data_t cp_r = Wr[VX] + cf;
  mrc_fld_data_t cm_r = Wr[VX] - cf;

  mrc_fld_data_t c_l =  mrc_fld_min(mrc_fld_min(cm_l, cm_r), 0.);
  mrc_fld_data_t c_r =  mrc_fld_max(mrc_fld_max(cp_l, cp_r), 0.);

  fluxes(Fl, Ul, Wl, i);
  fluxes(Fr, Ur, Wr, i);

  for (int m = 0; m < s_n_comps; m++) {
    F[m] = ((c_r * Fl[m] - c_l * Fr[m]) + (c_r * c_l * (Ur[m] - Ul[m]))) / (c_r - c_l);
  }
}

// ----------------------------------------------------------------------
// fluxes_hllc

static void
fluxes_hllc(mrc_fld_data_t F[5], mrc_fld_data_t Ul[5], mrc_fld_data_t Ur[5],
	    mrc_fld_data_t Wl[5], mrc_fld_data_t Wr[5], int i)
{
  assert(s_opt_eqn == OPT_EQN_MHD_SCONS || s_opt_eqn == OPT_EQN_HD); 

  mrc_fld_data_t Fl[5], Fr[5];

  fluxes(Fl, Ul, Wl, i);
  fluxes(Fr, Ur, Wr, i);

  mrc_fld_data_t cs2;

  cs2 = s_gamma * Wl[PP] / Wl[RR];
  mrc_fld_data_t cpv_l = Wl[VX] + sqrtf(cs2);
  mrc_fld_data_t cmv_l = Wl[VX] - sqrtf(cs2); 

  cs2 = s_gamma * Wr[PP] / Wr[RR];
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
// hlld_calc_state_s

static void
hlld_calc_state_s(mrc_fld_data_t Ws[], mrc_fld_data_t U[], mrc_fld_data_t W[], 
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
fluxes_hlld(mrc_fld_data_t F[], mrc_fld_data_t Ul[], mrc_fld_data_t Ur[],
	    mrc_fld_data_t Wl[], mrc_fld_data_t Wr[], int i)
{
  assert(s_opt_eqn == OPT_EQN_MHD_FCONS);

    mrc_fld_data_t Fl[s_n_comps], Fr[s_n_comps];
    mrc_fld_data_t bb, cs2, as2, cf;
    
    bb = sqr(Wl[BX]) + sqr(Wl[BY]) + sqr(Wl[BZ]);
    cs2 = s_gamma * Wl[PP] / Wl[RR];
    as2 = bb / Wl[RR]; 
    cf = sqrtf(.5 * (cs2 + as2 + sqrtf(sqr(as2 + cs2)
				       - (4. * sqr(sqrt(cs2) * Wl[BX]) / Wl[RR]))));       
    
    fluxes(Fl, Ul, Wl, i);
    mrc_fld_data_t cpv_l = Wl[VX] + cf;
    mrc_fld_data_t cmv_l = Wl[VX] - cf; 
    
    bb = sqr(Wr[BX]) + sqr(Wr[BY]) + sqr(Wr[BZ]);
    cs2 = s_gamma * Wr[PP] / Wr[RR];
    as2 = bb / Wr[RR]; 
    cf = sqrtf(.5 * (cs2 + as2 + sqrtf(sqr(as2 + cs2)
				       - (4. * sqr(sqrt(cs2) * Wr[BX]) / Wr[RR]))));     
    
    fluxes(Fr, Ur, Wr, i);
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
    
    mrc_fld_data_t Urs[s_n_comps], Uls[s_n_comps], Wls[s_n_comps], Wrs[s_n_comps]; 
    mrc_fld_data_t Urss[s_n_comps], Ulss[s_n_comps], Wlss[s_n_comps], Wrss[s_n_comps];
    
    hlld_calc_state_s(Wls, Ul, Wl, SL, SM, sPt, SLmUL);
    hlld_calc_state_s(Wrs, Ur, Wr, SR, SM, sPt, SRmUR);
    
    // MK eq. 49 
    Wlss[RR] = Wls[RR];
    Wrss[RR] = Wrs[RR];
    
    // MK eq. 50 
    mrc_fld_data_t ssPt _mrc_unused = sPt; // FIXME!!!
    
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
        
    for (int m = 0; m < s_n_comps; m++) {
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
// mhd_riemann

static void _mrc_unused
mhd_riemann(fld1d_state_t F, fld1d_state_t U_l, fld1d_state_t U_r,
	    fld1d_state_t W_l, fld1d_state_t W_r, int ib, int ie)
{
  // if applicable, solve GLM part of equation first, which will then be used
  // in fluxes()
  mhd_divb_glm_riemann(U_l, U_r, W_l, W_r, ib, ie);

  if (s_opt_riemann == OPT_RIEMANN_RUSANOV) {
    for (int i = ib; i < ie; i++) {
      fluxes_rusanov(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
		     &F1S(W_l, 0, i), &F1S(W_r, 0, i), i);
    }
  } else if (s_opt_riemann == OPT_RIEMANN_HLL) {
    for (int i = ib; i < ie; i++) {
      fluxes_hll(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
		 &F1S(W_l, 0, i), &F1S(W_r, 0, i), i);
    }
  } else if (s_opt_riemann == OPT_RIEMANN_HLLC) {
    for (int i = ib; i < ie; i++) {
      fluxes_hllc(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
		  &F1S(W_l, 0, i), &F1S(W_r, 0, i), i);
    }
  } else if (s_opt_riemann == OPT_RIEMANN_HLLD) {
    for (int i = ib; i < ie; i++) {
      fluxes_hlld(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
		  &F1S(W_l, 0, i), &F1S(W_r, 0, i), i);
    }
  } else {
    assert(0);
  }
}

#endif
