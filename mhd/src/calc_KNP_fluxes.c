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

#define CAP(m)  sqrtf( ( \
	     sqr(MRC_F3(u_p[m], _B1X,  (ix) - ((m) == 0), (iy) - ((m) == 1), (iz) - ((m) == 2))) + \
	     sqr(MRC_F3(u_p[m], _B1Y,  (ix) - ((m) == 0), (iy) - ((m) == 1), (iz) - ((m) == 2))) + \
	     sqr(MRC_F3(u_p[m], _B1Z,  (ix) - ((m) == 0), (iy) - ((m) == 1), (iz) - ((m) == 2)))) * rhoip)

#define CAM(m) sqrtf( (  sqr(MRC_F3(u_m[m], _B1X, ix,iy,iz))+ \
	     sqr(MRC_F3(u_m[m], _B1Y, ix,iy,iz)) + sqr(MRC_F3(u_m[m], _B1Z, ix,iy,iz)))* rhoim)

#define CFP(m)  sqrtf( .5f * ( sqr(csp) + sqr(cAp) + sqrtf( sqr( sqr(cAp) + sqr(csp) ) - \
              (4.f * sqr(csp * MRC_F3(u_p[m], _B1X+m, \
              (ix) - ((m) == 0), (iy) - ((m) == 1), (iz) - ((m) == 2))) * rhoip))));

#define CFM(m)  sqrtf( .5f * ( sqr(csm) + sqr(cAm) + sqrtf( sqr( sqr(cAm) + sqr(csm) ) - \
	      (4.f * sqr(csm * MRC_F3(u_m[m], _B1X+m, ix,iy,iz)) * rhoim))));


// ----------------------------------------------------------------------
// calc_KNP_fluxes 
// A. Kurganov, S. Noelle, G. Petrova, SIAM J. Sci. Comput. 23 (2001) 707.
// (Ziegler 2004 section 3.1) 

void
calc_KNP_fluxes(struct ggcm_mhd *mhd, struct mrc_fld *flux[3],
	    struct mrc_fld *flux_p[3], struct mrc_fld *flux_m[3],
	    struct mrc_fld *_u,
	    struct mrc_fld *u_p[3], struct mrc_fld *u_m[3],
	    struct mrc_fld *E_p[3], struct mrc_fld *E_m[3])
{
  float gamma = mhd->par.gamm;
  float d_i = mhd->par.d_i;
  float mpermi = 1.f;

  /*
  struct mrc_fld *flux[3], *flux_p[3], *flux_m[3];
  struct mrc_fld *u_p[3], *u_m[3], *E_p[3], *E_m[3];
  */
  struct mrc_fld *u = mrc_fld_get_as(_u, "float");
  
  /*
  for (int f = 0; f < 3; f++) {
    flux[f] = mrc_fld_get_as(_flux[f], "float");
    flux_p[f] = mrc_fld_get_as(_flux_p[f], "float");
    flux_m[f] = mrc_fld_get_as(_flux_m[f], "float");
    u_p[f] = mrc_fld_get_as(_u_p[f], "float");
    u_m[f] = mrc_fld_get_as(_u_m[f], "float");
    E_p[f] = mrc_fld_get_as(_E_p[f], "float");
    E_m[f] = mrc_fld_get_as(_E_m[f], "float");
  }
  */
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
       
    float cAp = CAP(0); float cAm = CAM(0);         
    float cfp = CFP(0); float cfm = CFM(0);


    float cwp = d_i * cAp * sqrtf(1./MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) * bdx1[ix+2]  ;
    float cwm = d_i * cAm * sqrtf(1./MRC_F3(u_m[0], _RR1, ix,iy,iz)) * bdx1[ix+2]  ; 

#if incws == 0
    cwp = 0;
    cwm = 0;
#endif
 
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

    cAp = CAP(1); cAm = CAM(1);         
    cfp = CFP(1); cfm =  CFM(1);

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
   			
    cAp = CAP(2); cAm = CAM(2);
    cfp = CFP(2); cfm = CFM(2);
 
    cwp = d_i * cAp * sqrtf(rhoip) * bdz1[iz+2] ;
    cwm = d_i * cAm * sqrtf(rhoim) * bdz1[iz+2] ; 

#if incws == 0
    cwp = 0;
    cwm = 0;
#endif

    float cp = fmaxf(fmaxf( (MRC_F3(u_p[2], _RV1Z, ix,iy,iz-1) / MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) + cfp + cwp,
			    (MRC_F3(u_m[2], _RV1Z, ix,iy,iz) / MRC_F3(u_m[2], _RR1, ix,iy,iz)) + cfm + cwm), 0.f);
    float cm = fminf(fminf( (MRC_F3(u_p[2], _RV1Z, ix,iy,iz-1) / MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) - cfp - cwp,
			    (MRC_F3(u_m[2], _RV1Z, ix,iy,iz) / MRC_F3(u_m[2], _RR1, ix,iy,iz)) - cfm - cwm), 0.f);
  

#if KT == 1
    cp = fmaxf(cp,-cm);
    cm=-cp;
#endif    

    assert(isfinite(ap));
    assert(isfinite(am));


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

    for (int m = 0; m <= _UU1; m++) {
      FLUX(flux, 0, m, ix,iy,iz) =
	(ap * FLUX(flux_p, 0, m, ix-1,iy,iz) - am * FLUX(flux_m, 0, m, ix,iy,iz)) / (ap - am) +
	(ap * am) / (ap - am) * (MRC_F3(u_m[0], m, ix ,iy,iz) - MRC_F3(u_p[0], m, ix-1,iy,iz));
      FLUX(flux, 1, m, ix,iy,iz) = 
	(bp * FLUX(flux_p, 1, m, ix,iy-1,iz) - bm * FLUX(flux_m, 1, m, ix,iy,iz)) / (bp - bm) +
	(bp * bm) / (bp - bm) * (MRC_F3(u_m[1], m, ix,iy ,iz) - MRC_F3(u_p[1], m, ix,iy-1 ,iz));
      FLUX(flux, 2, m, ix,iy,iz) =   
	(cp * FLUX(flux_p, 2, m, ix,iy,iz-1) - cm * FLUX(flux_m, 2, m, ix,iy,iz)) / (cp - cm) +
	(cp * cm) / (cp - cm) * (MRC_F3(u_m[2], m, ix,iy,iz ) - MRC_F3(u_p[2], m, ix,iy,iz-1));
    }  
     
  } mrc_fld_foreach_end;
 
  mrc_fld_put_as(u, _u);
  /*
  for (int f = 0; f < 3; f++) {
    mrc_fld_put_as(flux[f], _flux[f]);
    mrc_fld_put_as(flux_p[f], _flux_p[f]);
    mrc_fld_put_as(flux_m[f], _flux_m[f]);
    mrc_fld_put_as(u_p[f], _u_p[f]);
    mrc_fld_put_as(u_m[f], _u_m[f]);
    mrc_fld_put_as(E_p[f], _E_p[f]);
    mrc_fld_put_as(E_m[f], _E_m[f]);
  }
  */
}
