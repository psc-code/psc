#include "ggcm_mhd_step_cweno_private.h" 

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>

#include <assert.h>
#include <stdio.h>
#include <math.h>

#define U_M(var, m) MRC_F3(u_m[m], var, ix, iy, iz)
#define U_P(var, m) MRC_F3(u_p[m], var, \
		     (ix) - ((m) == 0), \
		     (iy) - ((m) == 1),	\
		     (iz) - ((m) == 2))

#if SEMICONSV 

#define P_M(ndx) (gamma - 1.f) * ( U_M( _UU1, ndx) - ((.5f / U_M( _RR1, ndx) ) * \
		 (sqr(U_M( _RV1X, ndx)) + sqr(U_M( _RV1Y, ndx)) + sqr(U_M( _RV1Z, ndx)) )))

#define P_P(ndx) (gamma - 1.f) * ( U_P( _UU1, ndx) - ((.5f / U_P( _RR1, ndx) ) * \
		 (sqr(U_P( _RV1X, ndx)) + sqr(U_P( _RV1Y, ndx)) + sqr(U_P( _RV1Z, ndx)) )))
#else 

#define P_M(ndx) (gamma - 1.f) * ( U_M( _UU1, ndx) - ((.5f / U_M( _RR1, ndx) ) * \
		   (sqr(U_M( _RV1X, ndx)) + sqr(U_M( _RV1Y, ndx)) + sqr(U_M( _RV1Z, ndx)) )) - \
	    (.5f * (sqr(U_M( _B1X, ndx)) + sqr(U_M(_B1Y, ndx)) + sqr(U_M(_B1Z, ndx)))))

#define P_P(ndx) (gamma - 1.f) * ( U_P( _UU1, ndx) - ((.5f / U_P( _RR1, ndx) ) * \
		   (sqr(U_P( _RV1X, ndx)) + sqr(U_P( _RV1Y, ndx)) + sqr(U_P( _RV1Z, ndx)) )) - \
	    (.5f * (sqr(U_P( _B1X, ndx))  + sqr(U_P( _B1Y, ndx))  + sqr(U_P(_B1Z, ndx)))))
#endif 



// ----------------------------------------------------------------------
// calc_KNP_fluxes 
// A. Kurganov, S. Noelle, G. Petrova, SIAM J. Sci. Comput. 23 (2001) 707.
// (Ziegler 2004 section 3.1) 

void
calc_KNP_fluxes(struct ggcm_mhd *mhd, struct mrc_fld *_flux[3],
	    struct mrc_fld *_flux_p[3], struct mrc_fld *_flux_m[3],
	    struct mrc_fld *_u,
	    struct mrc_fld *_u_p[3], struct mrc_fld *_u_m[3],
	    struct mrc_fld *_E_p[3], struct mrc_fld *_E_m[3])
{
  float gamma = mhd->par.gamm;
  float d_i = mhd->par.d_i;
  float mpermi = 1.f;

  struct mrc_fld *flux[3], *flux_p[3], *flux_m[3];
  struct mrc_fld *u_p[3], *u_m[3], *E_p[3], *E_m[3];
  struct mrc_fld *u = mrc_fld_get_as(_u, "mhd_fc_float");
  for (int f = 0; f < 3; f++) {
    flux[f] = mrc_fld_get_as(_flux[f], "float");
    flux_p[f] = mrc_fld_get_as(_flux_p[f], "float");
    flux_m[f] = mrc_fld_get_as(_flux_m[f], "float");
    u_p[f] = mrc_fld_get_as(_u_p[f], "float");
    u_m[f] = mrc_fld_get_as(_u_m[f], "float");
    E_p[f] = mrc_fld_get_as(_E_p[f], "float");
    E_m[f] = mrc_fld_get_as(_E_m[f], "float");
  }

  float *bdx1 = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD1);
  float *bdy1 = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD1);
  float *bdz1 = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD1);

  mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
    
   // Coeffiecents ap and am   
 
    float rhoip = 1.f / MRC_F3(u_p[0], _RR1, ix-1,iy,iz);	
    float rhoim = 1.f / MRC_F3(u_m[0], _RR1, ix,iy,iz);	
    
    float ppp = P_P(0);
    float ppm = P_M(0);
      
    float csp = sqrtf((gamma * ppp) / (MRC_F3(u_p[0], _RR1, ix-1,iy,iz)));
    float csm = sqrtf((gamma * ppm) / (MRC_F3(u_m[0], _RR1, ix,iy,iz)));   
       
    float cAp = sqrtf( (sqr(MRC_F3(u_p[0], _B1X, ix-1,iy,iz))+
			sqr(MRC_F3(u_p[0], _B1Y, ix-1,iy,iz))+
			sqr(MRC_F3(u_p[0], _B1Z, ix-1,iy,iz)))/ MRC_F3( u_p[0], _RR1, ix-1,iy,iz) );

    float cAm = sqrtf( (sqr(MRC_F3(u_m[0], _B1X, ix,iy,iz))+
			sqr(MRC_F3(u_m[0], _B1Y, ix,iy,iz))+
			sqr(MRC_F3(u_m[0], _B1Z, ix,iy,iz)))/ MRC_F3(u_m[0], _RR1, ix,iy,iz) );

#if 1
    float tmpp = sqr(csp) + sqr(cAp);
    float cfp = sqrtf( 0.5 * ( tmpp + sqrtf( sqr( sqr(cAp) + sqr(csp) )
					     - (4. * mpermi * sqr(csp * MRC_F3(u, _B1X, ix-1,iy,iz)) /  
						MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) )) );      
    float tmpm = sqr(csm) + sqr(cAm);
    float cfm = sqrtf( 0.5 * ( tmpm + sqrtf( sqr( sqr(cAm) + sqr(csm) )
					     - (4.* mpermi * sqr(csm * MRC_F3(u, _B1X, ix,iy,iz)) /  
						MRC_F3(u_m[0], _RR1, ix,iy,iz))  )) );
    
#else 

    float tmpp = sqr(csp) + sqr(cAp);
    float cfp = sqrtf( 0.5 * ( tmpp + sqrtf( sqr( sqr(cAp) + sqr(csp) )
					     - (4. * mpermi * sqr(csp * MRC_F3(u_p[0], _B1X, ix-1,iy,iz)) /  
						MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) )) );      
    float tmpm = sqr(csm) + sqr(cAm);
    float cfm = sqrtf( 0.5 * ( tmpm + sqrtf( sqr( sqr(cAm) + sqr(csm) )
					     - (4.* mpermi * sqr(csm * MRC_F3(u_p[0], _B1X, ix,iy,iz)) /  
						MRC_F3(u_m[0], _RR1, ix,iy,iz))  )) );
    
#endif 

    float cwp = d_i * cAp * sqrtf(1./MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) * bdx1[ix+2]  ;
    float cwm = d_i * cAm * sqrtf(1./MRC_F3(u_m[0], _RR1, ix,iy,iz)) * bdx1[ix+2]  ; 

#if incws == 0
    cwp = 0;
    cwm = 0;
#endif
    float acfm = cfm; 
    float acfp = cfp; 
    float acAp = cAp; 
    float acAm = cAm;
    float acsp = csp; 
    float acsm = csm;
    float acpp = U_P( _UU1, 0) ; 
    float acpm = U_M( _UU1, 0) ; 
    //printf("cAp %f cwp %f %f \n", cAp, d_i*(cAp*M_PI*gmhd->crdx[BD1][ix+2]) 
    // sqrtf(MRC_F3(u_p[0], _RR1, ix-1,iy,iz)),gmhd->crdx[BD1][ix+2]  ); 

    float ap = fmaxf(fmaxf((MRC_F3(u_p[0], _RV1X, ix-1,iy,iz) / MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) + cfp + cwp,
			   (MRC_F3(u_m[0], _RV1X, ix,iy,iz) / MRC_F3(u_m[0], _RR1, ix,iy,iz)) + cfm + cwm ), 0.f);

    float am = fminf(fminf( (MRC_F3(u_p[0], _RV1X, ix-1,iy,iz) / MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) - cfp - cwp,
			    (MRC_F3(u_m[0], _RV1X, ix,iy,iz) / MRC_F3(u_m[0], _RR1,ix,iy,iz)) - cfm - cwm), 0.f);

#if KT == 1
    ap = fmaxf(ap,-am);
    am=-ap;
#endif    

    // Coeffiecents bp and bm 
    rhoip = 1.f / MRC_F3(u_p[1], _RR1, ix,iy-1,iz);
    rhoim = 1.f / MRC_F3(u_m[1], _RR1, ix,iy,iz); 
  
    ppp = P_P(1); 
    ppm = P_M(1); 

    csp = sqrtf((gamma * ppp) / (MRC_F3(u_p[1], _RR1, ix,iy-1,iz)));
    csm = sqrtf((gamma * ppm) / (MRC_F3(u_m[1], _RR1, ix,iy,iz)));

    cAp = sqrtf( (sqr(MRC_F3(u, _B1X, ix,iy-1,iz))+
		  sqr(MRC_F3(u, _B1Y, ix,iy-1,iz))+
		  sqr(MRC_F3(u, _B1Z, ix,iy-1,iz)))/ MRC_F3(u_p[1], _RR1, ix,iy-1,iz) );

    cAm = sqrtf( (sqr(MRC_F3(u, _B1X, ix,iy,iz))+
		  sqr(MRC_F3(u, _B1Y, ix,iy,iz))+
		  sqr(MRC_F3(u, _B1Z, ix,iy,iz)))/ MRC_F3(u_m[1], _RR1, ix,iy,iz) );
    			  	  
    tmpp = sqr(csp) + sqr(cAp);
    cfp = sqrtf( 0.5 * (tmpp + sqrtf( sqr( sqr(cAp) + sqr(csp) ) - 
				      (4.f * mpermi * sqr(csp * MRC_F3(u, _B1Y, ix,iy-1,iz)) * rhoip))));      
    tmpm = sqr(csm) + sqr(cAm);
    cfm = sqrtf( 0.5 * (tmpm + sqrtf( sqr( sqr(cAm) + sqr(csm) ) -  
				      (4.f * mpermi * sqr(csm * MRC_F3(u, _B1Y, ix,iy,iz)) * rhoim)))); 

    cwp =  d_i * cAp * sqrtf(1./MRC_F3(u_p[1], _RR1, ix,iy-1,iz)) * bdy1[iy+2];
    cwm =  d_i * cAm * sqrtf(1./MRC_F3(u_m[1], _RR1, ix,iy,iz)) * bdy1[iy+2]; 
    
#if incws == 0
    cwp = 0;
    cwm = 0;
#endif
    
    float bp = fmaxf(fmaxf( (MRC_F3(u_p[1], _RV1Y, ix,iy-1,iz) / MRC_F3(u_p[1], _RR1, ix,iy-1,iz)) + cfp + cwp,
			    (MRC_F3(u_m[1], _RV1Y, ix,iy,iz) / MRC_F3(u_m[1], _RR1, ix,iy,iz)) + cfm + cwm), 0.f);
    float bm = fminf(fminf( (MRC_F3(u_p[1], _RV1Y, ix,iy-1,iz) / MRC_F3(u_p[1], _RR1, ix,iy-1,iz)) - cfp - cwp,
			    (MRC_F3(u_m[1], _RV1Y, ix,iy,iz) / MRC_F3(u_m[1], _RR1, ix,iy,iz)) - cfm - cwm), 0.f);
    
#if KT == 1
    bp = fmaxf(bp,-bm);
    bm=-bp;
#endif    

   // Coeffiecents cp and cm 
    rhoip = 1.f / MRC_F3(u_p[2], _RR1, ix,iy,iz-1);
    rhoim = 1.f / MRC_F3(u_m[2], _RR1, ix,iy,iz); 

    ppp = P_P(2); 
    ppm = P_M(2);


    csp = sqrtf((gamma * ppp) * rhoip);
    csm = sqrtf((gamma * ppm) * rhoim);
   			
    cAp = sqrtf( (sqr(MRC_F3(u_p[2], _B1X, ix,iy,iz-1))+
		  sqr(MRC_F3(u_p[2], _B1Y, ix,iy,iz-1))+
		  sqr(MRC_F3(u_p[2], _B1Z, ix,iy,iz-1))) * rhoip);// MRC_F3(u_p[2], _RR1, ix,iy,iz-1) );
    
    cAm = sqrtf( (sqr(MRC_F3(u_m[2], _B1X, ix,iy,iz))+
		  sqr(MRC_F3(u_m[2], _B1Y, ix,iy,iz))+
		  sqr(MRC_F3(u_m[2], _B1Z, ix,iy,iz))) * rhoim); // MRC_F3(u_m[2], _RR1, ix,iy,iz) );

    tmpp = sqr(csp) + sqr(cAp);
    cfp = sqrtf( 0.5f * (tmpp + sqrtf( sqr( sqr(cAp) + sqr(csp) ) - 
				       (4. * mpermi * sqr(csp * MRC_F3(u_p[2], _B1Z, ix,iy,iz-1)) *  
					rhoip )) ));      
    if (!isfinite(cfp)) {
      //mprintf("ix %d %d %d cfp = csp %g cAp %g rr %g ppp %g ppm %g tmpbe %g tmpke %g tmpee %g \n", ix,iy,iz, csp, cAp, MRC_F3(u_p[2], _RR1, ix,iy,iz-1), ppp, ppm, tmpbe, tmpke, tmpee);
    }

    tmpm = sqr(csm) + sqr(cAm);
    cfm = sqrtf( 0.5f * (tmpm + sqrtf( sqr( sqr(cAm) + sqr(csm) ) -
				       (4.* mpermi * sqr(csm * MRC_F3(u_m[2], _B1Z, ix,iy,iz)) *  
					rhoim)) ));

    cwp = d_i * cAp * sqrtf(rhoip) * bdz1[iz+2] ;
    cwm = d_i * cAm * sqrtf(rhoim) * bdz1[iz+2] ; 

    //printf("cAp %f cwp %f %f \n", cAp,cAp * M_PI * sqrtf(MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) * gmhd->crdz[BD1][iz+2] , ); 

#if incws == 0
    cwp = 0;
    cwm = 0;
#endif

    float cp = fmaxf(fmaxf( (MRC_F3(u_p[2], _RV1Z, ix,iy,iz-1) / MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) + cfp + cwp,
			    (MRC_F3(u_m[2], _RV1Z, ix,iy,iz) / MRC_F3(u_m[2], _RR1, ix,iy,iz)) + cfm + cwm), 0.f);
    float cm = fminf(fminf( (MRC_F3(u_p[2], _RV1Z, ix,iy,iz-1) / MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) - cfp - cwp,
			    (MRC_F3(u_m[2], _RV1Z, ix,iy,iz) / MRC_F3(u_m[2], _RR1, ix,iy,iz)) - cfm - cwm), 0.f);
    if (cp == 0. && cm == 0.) {
      //mprintf("ix %d %d %d cfp %g cwp %g\n", ix,iy,iz, cfp, cwp);
    }

#if KT == 1
    cp = fmaxf(cp,-cm);
    cm=-cp;
#endif    

    assert(isfinite(ap));
    assert(isfinite(am));


#if 0

    //    ap = 1e-3; am = -1e-3;
    // Flux of _EX,_EY,_EZ through the x faces
    FLUX(flux, 0, _EX, ix,iy,iz) = 
      (1.f/(ap - am)) * ( (ap*am) * ( MRC_F3(u_m[0], _B1X, ix,iy,iz) - MRC_F3(u_p[0], _B1X, ix-1,iy,iz)));
    FLUX(flux, 0, _EY, ix,iy,iz) = 
      (1.f/(ap - am)) * (   - ap * MRC_F3(E_p[0], 2, ix-1, iy,iz) + am *  MRC_F3(E_m[0], 2, ix,iy,iz) + 
			 (ap*am)*   ( MRC_F3(u_m[0], _B1Y, ix,iy,iz) -  MRC_F3(u_p[0], _B1Y, ix-1,iy,iz) ));
    FLUX(flux, 0, _EZ, ix,iy,iz) = 
      (1.f/(ap - am)) * (     ap * MRC_F3(E_p[0], 1, ix-1,iy,iz) - am *  MRC_F3(E_m[0], 1, ix,iy,iz) + 
			 (ap*am) * ( MRC_F3(u_m[0], _B1Z, ix,iy,iz) - MRC_F3(u_p[0], _B1Z, ix-1,iy,iz) ));  
    
    assert(isfinite(bp));
    assert(isfinite(bm));
    //bp = 1e-3; bm = -1e-3;
    // flux of _EX,_EY,_EZ through the y faces    
    FLUX(flux, 1, _EX, ix,iy,iz) =
       (1.f/(bp - bm)) * (     bp  *  MRC_F3(E_p[1], 2, ix,iy-1,iz) - bm * MRC_F3(E_m[1], 2, ix,iy,iz) + 
		         (bp * bm) *  ( MRC_F3(u_m[1], _B1X, ix,iy,iz) -    MRC_F3(u_p[1], _B1X, ix,iy-1,iz) ));
    FLUX(flux, 1, _EY, ix,iy,iz) =
       (1.f/(bp - bm)) * ( (bp * bm)* ( MRC_F3(u_m[1], _B1Y, ix,iy,iz) - MRC_F3(u_p[1], _B1Y, ix,iy-1,iz) ));    
    FLUX(flux, 1, _EZ, ix,iy,iz) =
       (1.f/(bp - bm)) * (    - bp * MRC_F3(E_p[1], 0, ix,iy-1,iz) + bm *  MRC_F3(E_m[1], 0, ix,iy,iz) + 
		         (bp * bm) *  ( MRC_F3(u_m[1], _B1Z, ix,iy,iz) -     MRC_F3(u_p[1], _B1Z, ix,iy-1,iz) ));  
    
    assert(isfinite(cp));
    assert(isfinite(cm));
    //cp = 1e-3; cm = -1e-3;
    // flux of _EX,_EY,_EZ through the z faces
    FLUX(flux, 2, _EX, ix,iy,iz) = 
      (1.f/(cp - cm))*( - cp * MRC_F3(E_p[2], 1, ix,iy,iz-1) + cm * MRC_F3(E_m[2], 1, ix,iy,iz) + 
			 (cp * cm) * ( MRC_F3(u_m[2], _B1X, ix,iy,iz) - MRC_F3(u_p[2], _B1X, ix,iy,iz-1) ));
    FLUX(flux, 2, _EY, ix,iy,iz) =
      (1.f/(cp - cm)) *(  cp * MRC_F3(E_p[2], 0, ix,iy,iz-1) - cm * MRC_F3(E_m[2], 0, ix,iy,iz) +
		         (cp * cm) * ( MRC_F3(u_m[2], _B1Y, ix,iy,iz) -    MRC_F3(u_p[2], _B1Y, ix,iy,iz-1) ));
    FLUX(flux, 2, _EZ, ix,iy,iz) = 
      (1.f/(cp - cm))*(  (cp * cm) * ( MRC_F3(u_m[2], _B1Z, ix,iy,iz) - MRC_F3(u_p[2], _B1Z, ix,iy,iz-1) ));    

#else 
 //    ap = 1e-3; am = -1e-3;
    // Flux of _EX,_EY,_EZ through the x faces
    FLUX(flux, 0, _EX, ix,iy,iz) = 
      (1.f/(ap - am)) * ( (ap*am) * ( MRC_F3(u_m[0], _B1X, ix,iy,iz) - MRC_F3(u_p[0], _B1X, ix-1,iy,iz)));
    FLUX(flux, 0, _EY, ix,iy,iz) = 
      (1.f/(ap - am)) * (ap*am)*   ( MRC_F3(u_m[0], _B1Y, ix,iy,iz) -  MRC_F3(u_p[0], _B1Y, ix-1,iy,iz));
    FLUX(flux, 0, _EZ, ix,iy,iz) = 
      (1.f/(ap - am)) * (ap*am) * ( MRC_F3(u_m[0], _B1Z, ix,iy,iz) - MRC_F3(u_p[0], _B1Z, ix-1,iy,iz));  
     
    assert(isfinite(bp));
    assert(isfinite(bm));
    //bp = 1e-3; bm = -1e-3;
    // flux of _EX,_EY,_EZ through the y faces    
    FLUX(flux, 1, _EX, ix,iy,iz) =
      (1.f/(bp - bm)) * (bp * bm) *  ( MRC_F3(u_m[1], _B1X, ix,iy,iz) -    MRC_F3(u_p[1], _B1X, ix,iy-1,iz));
    FLUX(flux, 1, _EY, ix,iy,iz) =
      (1.f/(bp - bm)) * ( (bp * bm)* ( MRC_F3(u_m[1], _B1Y, ix,iy,iz) - MRC_F3(u_p[1], _B1Y, ix,iy-1,iz)));    
    FLUX(flux, 1, _EZ, ix,iy,iz) =
      (1.f/(bp - bm)) * (bp * bm) *  ( MRC_F3(u_m[1], _B1Z, ix,iy,iz) -     MRC_F3(u_p[1], _B1Z, ix,iy-1,iz));  
    
    assert(isfinite(cp));
    assert(isfinite(cm));
    //cp = 1e-3; cm = -1e-3;
    // flux of _EX,_EY,_EZ through the z faces
    FLUX(flux, 2, _EX, ix,iy,iz) = 
      (1.f/(cp - cm)) *  (cp * cm) * ( MRC_F3(u_m[2], _B1X, ix,iy,iz) - MRC_F3(u_p[2], _B1X, ix,iy,iz-1));
    FLUX(flux, 2, _EY, ix,iy,iz) =
      (1.f/(cp - cm)) *  (cp * cm) * ( MRC_F3(u_m[2], _B1Y, ix,iy,iz) -    MRC_F3(u_p[2], _B1Y, ix,iy,iz-1));;
    FLUX(flux, 2, _EZ, ix,iy,iz) = 
      (1.f/(cp - cm))*(  (cp * cm) * ( MRC_F3(u_m[2], _B1Z, ix,iy,iz) - MRC_F3(u_p[2], _B1Z, ix,iy,iz-1)));    


#endif






    for (int m = 0; m <= _UU1; m++) {
#if 0
      FLUX(flux, 0, m, ix,iy,iz) =
	(ap * FLUX(flux_p, 0, m, ix-1,iy,iz) - am * FLUX(flux_m, 0, m, ix,iy,iz)) / (ap - am) +
	(ap * am) / (ap - am) * (MRC_F3(u_m[0], m, ix ,iy,iz) - MRC_F3(u_p[0], m, ix-1,iy,iz));
      FLUX(flux, 1, m, ix,iy,iz) = 
	(bp * FLUX(flux_p, 1, m, ix,iy-1,iz) - bm * FLUX(flux_m, 1, m, ix,iy,iz)) / (bp - bm) +
	(bp * bm) / (bp - bm) * (MRC_F3(u_m[1], m, ix,iy ,iz) - MRC_F3(u_p[1], m, ix,iy-1 ,iz));
      FLUX(flux, 2, m, ix,iy,iz) =   
	(cp * FLUX(flux_p, 2, m, ix,iy,iz-1) - cm * FLUX(flux_m, 2, m, ix,iy,iz)) / (cp - cm) +
	(cp * cm) / (cp - cm) * (MRC_F3(u_m[2], m, ix,iy,iz ) - MRC_F3(u_p[2], m, ix,iy,iz-1));
#else
      FLUX(flux, 0, m, ix,iy,iz) = .5f * (FLUX(flux_p, 0, m, ix-1,iy,iz) + FLUX(flux_m, 0, m, ix,iy,iz));
      FLUX(flux, 1, m, ix,iy,iz) = .5f * (FLUX(flux_p, 1, m, ix-1,iy,iz) + FLUX(flux_m, 1, m, ix,iy,iz));
      FLUX(flux, 2, m, ix,iy,iz) = .5f * (FLUX(flux_p, 2, m, ix-1,iy,iz) + FLUX(flux_m, 2, m, ix,iy,iz));
#endif

      
      if (ix > 1 && ix < 64 &&
	  iy > 1 && iy < 64 &&
	  iz > 1 && iz < 64) {
	//assert(isfinite(FLUX(flux, 0, m, ix,iy,iz)));
	//	assert(isfinite(FLUX(flux, 1, m, ix,iy,iz)));

	//assert(isfinite(FLUX(flux, 0, m, ix,iy,iz)));
      }
      
    } 
  } mrc_fld_foreach_end;
 
  mrc_fld_put_as(u, _u);
  for (int f = 0; f < 3; f++) {
    mrc_fld_put_as(flux[f], _flux[f]);
    mrc_fld_put_as(flux_p[f], _flux_p[f]);
    mrc_fld_put_as(flux_m[f], _flux_m[f]);
    mrc_fld_put_as(u_p[f], _u_p[f]);
    mrc_fld_put_as(u_m[f], _u_m[f]);
    mrc_fld_put_as(E_p[f], _E_p[f]);
    mrc_fld_put_as(E_m[f], _E_m[f]);
  }
}
