
// FIXME, at least uppercase
#define sign(x) (( x > 0. ) - ( x < 0. ))

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
// fluxes_hydro

// FIXME, sc vs hydro is a kinda arbitrary distinction, but the former does not
// contain pressure, because that's what's needed for Jimmy-MHD ("c3")

static void
fluxes_hydro(mrc_fld_data_t F[5], mrc_fld_data_t U[5], mrc_fld_data_t W[5])
{
  F[RR]  = W[RR] * W[VX];
  F[RVX] = W[RR] * W[VX] * W[VX] + W[PP];
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
  mrc_fld_data_t cs2 = s_gamma * W[PP] / W[RR];
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
  cs2 = s_gamma * Wl[PP] / Wl[RR];
  mrc_fld_data_t cmsv_l = sqrtf(vv) + sqrtf(cs2);

  vv = sqr(Wr[VX]) + sqr(Wr[VY]) + sqr(Wr[VZ]);
  cs2 = s_gamma * Wr[PP] / Wr[RR];
  mrc_fld_data_t cmsv_r = sqrtf(vv) + sqrtf(cs2);

  mrc_fld_data_t lambda = .5 * (cmsv_l + cmsv_r);
  
  for (int m = 0; m < 5; m++) {
    F[m] = .5f * ((Fr[m] + Fl[m]) - lambda * (Ur[m] - Ul[m]));
  }
}

// ----------------------------------------------------------------------
// fluxes_rusanov_hydro

static void
fluxes_rusanov_hydro(mrc_fld_data_t F[5], mrc_fld_data_t Ul[5], mrc_fld_data_t Ur[5],
		     mrc_fld_data_t Wl[5], mrc_fld_data_t Wr[5])
{
  mrc_fld_data_t Fl[5], Fr[5];
  
  fluxes_hydro(Fl, Ul, Wl);
  fluxes_hydro(Fr, Ur, Wr);

  mrc_fld_data_t vv, cs2;
  vv = sqr(Wl[VX]) + sqr(Wl[VY]) + sqr(Wl[VZ]);
  cs2 = s_gamma * Wl[PP] / Wl[RR];
  mrc_fld_data_t cmsv_l = sqrtf(vv) + sqrtf(cs2);

  vv = sqr(Wr[VX]) + sqr(Wr[VY]) + sqr(Wr[VZ]);
  cs2 = s_gamma * Wr[PP] / Wr[RR];
  mrc_fld_data_t cmsv_r = sqrtf(vv) + sqrtf(cs2);

  mrc_fld_data_t lambda = .5 * (cmsv_l + cmsv_r);
  
  for (int m = 0; m < 5; m++) {
    F[m] = .5f * ((Fr[m] + Fl[m]) - lambda * (Ur[m] - Ul[m]));
  }
}

// ----------------------------------------------------------------------
// fluxes_hll_fc

static void
fluxes_hll_fc(mrc_fld_data_t F[8], mrc_fld_data_t Ul[8], mrc_fld_data_t Ur[8],
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

  mrc_fld_data_t c_l =  mrc_fld_min(mrc_fld_min(cm_l, cm_r), 0.); 
  mrc_fld_data_t c_r =  mrc_fld_max(mrc_fld_max(cp_l, cp_r), 0.); 

  for (int m = 0; m < 8; m++) {
    F[m] = ((c_r * Fl[m] - c_l * Fr[m]) + (c_r * c_l * (Ur[m] - Ul[m]))) / (c_r - c_l);
  }
}

// ----------------------------------------------------------------------
// fluxes_hll_hydro

static void
fluxes_hll_hydro(mrc_fld_data_t F[5], mrc_fld_data_t Ul[5], mrc_fld_data_t Ur[5],
		 mrc_fld_data_t Wl[5], mrc_fld_data_t Wr[5])
{
  mrc_fld_data_t Fl[5], Fr[5];

  fluxes_hydro(Fl, Ul, Wl);
  fluxes_hydro(Fr, Ur, Wr);

  mrc_fld_data_t cs2;

  cs2 = s_gamma * Wl[PP] / Wl[RR];
  mrc_fld_data_t cpv_l = Wl[VX] + sqrtf(cs2);
  mrc_fld_data_t cmv_l = Wl[VX] - sqrtf(cs2); 

  cs2 = s_gamma * Wr[PP] / Wr[RR];
  mrc_fld_data_t cpv_r = Wr[VX] + sqrtf(cs2);
  mrc_fld_data_t cmv_r = Wr[VX] - sqrtf(cs2); 

  mrc_fld_data_t SR =  fmaxf(fmaxf(cpv_l, cpv_r), 0.); 
  mrc_fld_data_t SL =  fminf(fminf(cmv_l, cmv_r), 0.); 

  //  mrc_fld_data_t lambda = .5 * (cmsv_l + cmsv_r);  
  for (int m = 0; m < 5; m++) {
    F[m] = ((SR * Fl[m] - SL * Fr[m]) + (SR * SL * (Ur[m] - Ul[m]))) / (SR - SL);
  }
}

// ----------------------------------------------------------------------
// fluxes_hllc_hydro

static void
fluxes_hllc_hydro(mrc_fld_data_t F[5], mrc_fld_data_t Ul[5], mrc_fld_data_t Ur[5],
		  mrc_fld_data_t Wl[5], mrc_fld_data_t Wr[5])
{
  mrc_fld_data_t Fl[5], Fr[5];

  fluxes_hydro(Fl, Ul, Wl);
  fluxes_hydro(Fr, Ur, Wr);

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
hlld_calc_state_s(mrc_fld_data_t Ws[8], mrc_fld_data_t U[8], mrc_fld_data_t W[8], 
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
// fluxes_hlld_fc
// 
// Miyoshi & Kusano (2005)

static void
fluxes_hlld_fc(mrc_fld_data_t F[8], mrc_fld_data_t Ul[8], mrc_fld_data_t Ur[8],
	       mrc_fld_data_t Wl[8], mrc_fld_data_t Wr[8])
{
    mrc_fld_data_t Fl[8], Fr[8];
    mrc_fld_data_t bb, cs2, as2, cf;
    
    bb = sqr(Wl[BX]) + sqr(Wl[BY]) + sqr(Wl[BZ]);
    cs2 = s_gamma * Wl[PP] / Wl[RR];
    as2 = bb / Wl[RR]; 
    cf = sqrtf(.5 * (cs2 + as2 + sqrtf(sqr(as2 + cs2)
				       - (4. * sqr(sqrt(cs2) * Wl[BX]) / Wl[RR]))));       
    
    fluxes_fc(Fl, Ul, Wl);
    mrc_fld_data_t cpv_l = Wl[VX] + cf;
    mrc_fld_data_t cmv_l = Wl[VX] - cf; 
    
    bb = sqr(Wr[BX]) + sqr(Wr[BY]) + sqr(Wr[BZ]);
    cs2 = s_gamma * Wr[PP] / Wr[RR];
    as2 = bb / Wr[RR]; 
    cf = sqrtf(.5 * (cs2 + as2 + sqrtf(sqr(as2 + cs2)
				       - (4. * sqr(sqrt(cs2) * Wr[BX]) / Wr[RR]))));     
    
    fluxes_fc(Fr, Ur, Wr);
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
// mhd_riemann_rusanov_run_hydro

static void _mrc_unused
mhd_riemann_rusanov_run_hydro(fld1d_state_t F,
			      fld1d_state_t U_l, fld1d_state_t U_r,
			      fld1d_state_t W_l, fld1d_state_t W_r,
			      int ldim, int l, int r, int dim)
{
  for (int i = -l; i < ldim + r; i++) {
    fluxes_rusanov_hydro(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
			 &F1S(W_l, 0, i), &F1S(W_r, 0, i));
  }
}

// ----------------------------------------------------------------------
// mhd_riemann_hll_run_hydro

static void _mrc_unused
mhd_riemann_hll_run_hydro(fld1d_state_t F,
			  fld1d_state_t U_l, fld1d_state_t U_r,
			  fld1d_state_t W_l, fld1d_state_t W_r,
			  int ldim, int l, int r, int dim)
{
  for (int i = -l; i < ldim + r; i++) {
    fluxes_hll_hydro(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
		     &F1S(W_l, 0, i), &F1S(W_r, 0, i));
  }
}

// ----------------------------------------------------------------------
// mhd_riemann_hllc_run_hydro

static void _mrc_unused
mhd_riemann_hllc_run_hydro(fld1d_state_t F,
			   fld1d_state_t U_l, fld1d_state_t U_r,
			   fld1d_state_t W_l, fld1d_state_t W_r,
			   int ldim, int l, int r, int dim)
{
  for (int i = -l; i < ldim + r; i++) {
    fluxes_hllc_hydro(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
		      &F1S(W_l, 0, i), &F1S(W_r, 0, i));
  }
}

// ----------------------------------------------------------------------
// mhd_riemann_hlld_run_fc

static void _mrc_unused
mhd_riemann_hlld_run_fc(fld1d_state_t F,
			fld1d_state_t U_l, fld1d_state_t U_r,
			fld1d_state_t W_l, fld1d_state_t W_r,
			int ldim, int l, int r, int dim)
{
  for (int i = -l; i < ldim + r; i++) {
    fluxes_hlld_fc(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
		   &F1S(W_l, 0, i), &F1S(W_r, 0, i));
  }
}

// ----------------------------------------------------------------------
// mhd_riemann_run_fc

static void _mrc_unused
mhd_riemann_run_fc(fld1d_state_t F, fld1d_state_t U_l, fld1d_state_t U_r,
		   fld1d_state_t W_l, fld1d_state_t W_r,
		   int ldim, int l, int r, int dim)
{
  if (s_opt_riemann == OPT_RIEMANN_RUSANOV) {
    for (int i = -l; i < ldim + r; i++) {
      fluxes_rusanov_fc(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
			&F1S(W_l, 0, i), &F1S(W_r, 0, i));
    }
  } else if (s_opt_riemann == OPT_RIEMANN_HLL) {
    for (int i = -l; i < ldim + r; i++) {
      fluxes_hll_fc(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
		    &F1S(W_l, 0, i), &F1S(W_r, 0, i));
    }
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mhd_riemann_run_sc

static void _mrc_unused
mhd_riemann_run_sc(fld1d_state_t F, fld1d_state_t U_l, fld1d_state_t U_r,
		   fld1d_state_t W_l, fld1d_state_t W_r,
		   int ldim, int l, int r, int dim)
{
  if (s_opt_riemann == OPT_RIEMANN_RUSANOV) {
    for (int i = -l; i < ldim + r; i++) {
      fluxes_rusanov_sc(&F1S(F, 0, i), &F1S(U_l, 0, i), &F1S(U_r, 0, i),
			&F1S(W_l, 0, i), &F1S(W_r, 0, i));
    }
  } else {
    assert(0);
  }
}
