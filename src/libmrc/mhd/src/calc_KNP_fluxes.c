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

#define P_M(ndx) (gamma - 1.f) * ( U_M( UU, ndx) - ((.5f / U_M( RR, ndx) ) * \
		 (sqr(U_M( RVX, ndx)) + sqr(U_M( RVY, ndx)) + sqr(U_M( RVZ, ndx)) )))

#define P_P(ndx) (gamma - 1.f) * ( U_P( UU, ndx) - ((.5f / U_P( RR, ndx) ) * \
		 (sqr(U_P( RVX, ndx)) + sqr(U_P( RVY, ndx)) + sqr(U_P( RVZ, ndx)) )))
#else 

#define P_M(ndx) (gamma - 1.f) * ( U_M( UU, ndx) - ((.5f / U_M( RR, ndx) ) * \
		   (sqr(U_M( RVX, ndx)) + sqr(U_M( RVY, ndx)) + sqr(U_M( RVZ, ndx)) )) - \
	    (.5f * (sqr(U_M( BX, ndx)) + sqr(U_M(BY, ndx)) + sqr(U_M(BZ, ndx)))))

#define P_P(ndx) (gamma - 1.f) * ( U_P( UU, ndx) - ((.5f / U_P( RR, ndx) ) * \
		   (sqr(U_P( RVX, ndx)) + sqr(U_P( RVY, ndx)) + sqr(U_P( RVZ, ndx)) )) - \
	    (.5f * (sqr(U_P( BX, ndx))  + sqr(U_P( BY, ndx))  + sqr(U_P(BZ, ndx)))))
#endif 

#define CAP(m)  sqrtf( ( \
	     sqr(MRC_F3(u_p[m], BX,  (ix) - ((m) == 0), (iy) - ((m) == 1), (iz) - ((m) == 2))) + \
	     sqr(MRC_F3(u_p[m], BY,  (ix) - ((m) == 0), (iy) - ((m) == 1), (iz) - ((m) == 2))) + \
	     sqr(MRC_F3(u_p[m], BZ,  (ix) - ((m) == 0), (iy) - ((m) == 1), (iz) - ((m) == 2)))) * rhoip)

#define CAM(m) sqrtf( (  sqr(MRC_F3(u_m[m], BX, ix,iy,iz))+ \
	     sqr(MRC_F3(u_m[m], BY, ix,iy,iz)) + sqr(MRC_F3(u_m[m], BZ, ix,iy,iz)))* rhoim)

#define CFP(m)  sqrtf( .5f * ( sqr(csp) + sqr(cAp) + sqrtf( sqr( sqr(cAp) + sqr(csp) ) - \
              (4.f * sqr(csp * MRC_F3(u_p[m], BX+m, \
              (ix) - ((m) == 0), (iy) - ((m) == 1), (iz) - ((m) == 2))) * rhoip))));

#define CFM(m)  sqrtf( .5f * ( sqr(csm) + sqr(cAm) + sqrtf( sqr( sqr(cAm) + sqr(csm) ) - \
	      (4.f * sqr(csm * MRC_F3(u_m[m], BX+m, ix,iy,iz)) * rhoim))));


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

  struct mrc_fld *u = mrc_fld_get_as(_u, "float");
  
  float *bdx1 = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD1);
  float *bdy1 = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD1);
  float *bdz1 = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD1);

  mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
    
   // Coeffiecents ap and am   
 
    float rhoip = 1.f / MRC_F3(u_p[0], RR, ix-1,iy,iz);	
    float rhoim = 1.f / MRC_F3(u_m[0], RR, ix,iy,iz);	
    
    float ppp = P_P(0);
    float ppm = P_M(0);
      
    float csp = sqrtf((gamma * ppp) / (MRC_F3(u_p[0], RR, ix-1,iy,iz)));
    float csm = sqrtf((gamma * ppm) / (MRC_F3(u_m[0], RR, ix,iy,iz)));   
       
    float cAp = CAP(0); float cAm = CAM(0);         
    float cfp = CFP(0); float cfm = CFM(0);


    float cwp = d_i * cAp * sqrtf(1./MRC_F3(u_p[0], RR, ix-1,iy,iz)) * bdx1[ix+2]  ;
    float cwm = d_i * cAm * sqrtf(1./MRC_F3(u_m[0], RR, ix,iy,iz)) * bdx1[ix+2]  ; 

#if incws == 0
    cwp = 0;
    cwm = 0;
#endif
 
    float ap = fmaxf(fmaxf((MRC_F3(u_p[0], RVX, ix-1,iy,iz) / MRC_F3(u_p[0], RR, ix-1,iy,iz)) + cfp + cwp,
			   (MRC_F3(u_m[0], RVX, ix,iy,iz) / MRC_F3(u_m[0], RR, ix,iy,iz)) + cfm + cwm ), 0.f);

    float am = fminf(fminf( (MRC_F3(u_p[0], RVX, ix-1,iy,iz) / MRC_F3(u_p[0], RR, ix-1,iy,iz)) - cfp - cwp,
			    (MRC_F3(u_m[0], RVX, ix,iy,iz) / MRC_F3(u_m[0], RR,ix,iy,iz)) - cfm - cwm), 0.f);

#if KT == 1
    ap = fmaxf(ap,-am);
    am=-ap;
#endif    

    // Coeffiecents bp and bm 
    rhoip = 1.f / MRC_F3(u_p[1], RR, ix,iy-1,iz);
    rhoim = 1.f / MRC_F3(u_m[1], RR, ix,iy,iz); 
  
    ppp = P_P(1); 
    ppm = P_M(1); 

    csp = sqrtf((gamma * ppp) / (MRC_F3(u_p[1], RR, ix,iy-1,iz)));
    csm = sqrtf((gamma * ppm) / (MRC_F3(u_m[1], RR, ix,iy,iz)));

    cAp = CAP(1); cAm = CAM(1);         
    cfp = CFP(1); cfm =  CFM(1);

    cwp =  d_i * cAp * sqrtf(1./MRC_F3(u_p[1], RR, ix,iy-1,iz)) * bdy1[iy+2];
    cwm =  d_i * cAm * sqrtf(1./MRC_F3(u_m[1], RR, ix,iy,iz)) * bdy1[iy+2]; 
    
#if incws == 0
    cwp = 0;
    cwm = 0;
#endif
    
    float bp = fmaxf(fmaxf( (MRC_F3(u_p[1], RVY, ix,iy-1,iz) / MRC_F3(u_p[1], RR, ix,iy-1,iz)) + cfp + cwp,
			    (MRC_F3(u_m[1], RVY, ix,iy,iz) / MRC_F3(u_m[1], RR, ix,iy,iz)) + cfm + cwm), 0.f);
    float bm = fminf(fminf( (MRC_F3(u_p[1], RVY, ix,iy-1,iz) / MRC_F3(u_p[1], RR, ix,iy-1,iz)) - cfp - cwp,
			    (MRC_F3(u_m[1], RVY, ix,iy,iz) / MRC_F3(u_m[1], RR, ix,iy,iz)) - cfm - cwm), 0.f);
    
#if KT == 1
    bp = fmaxf(bp,-bm);
    bm=-bp;
#endif    

   // Coeffiecents cp and cm 
    rhoip = 1.f / MRC_F3(u_p[2], RR, ix,iy,iz-1);
    rhoim = 1.f / MRC_F3(u_m[2], RR, ix,iy,iz); 

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

    float cp = fmaxf(fmaxf( (MRC_F3(u_p[2], RVZ, ix,iy,iz-1) / MRC_F3(u_p[2], RR, ix,iy,iz-1)) + cfp + cwp,
			    (MRC_F3(u_m[2], RVZ, ix,iy,iz) / MRC_F3(u_m[2], RR, ix,iy,iz)) + cfm + cwm), 0.f);
    float cm = fminf(fminf( (MRC_F3(u_p[2], RVZ, ix,iy,iz-1) / MRC_F3(u_p[2], RR, ix,iy,iz-1)) - cfp - cwp,
			    (MRC_F3(u_m[2], RVZ, ix,iy,iz) / MRC_F3(u_m[2], RR, ix,iy,iz)) - cfm - cwm), 0.f);
  

#if KT == 1
    cp = fmaxf(cp,-cm);
    cm=-cp;
#endif    

    assert(isfinite(ap));
    assert(isfinite(am));


    //    ap = 1e-3; am = -1e-3;
    // Flux of _EX,_EY,_EZ through the x faces
    FLUX(flux, 0, _EX, ix,iy,iz) = 
      (1.f/(ap - am)) * ( (ap*am) * ( MRC_F3(u_m[0], BX, ix,iy,iz) - MRC_F3(u_p[0], BX, ix-1,iy,iz)));
    FLUX(flux, 0, _EY, ix,iy,iz) = 
      (1.f/(ap - am)) * (   - ap * MRC_F3(E_p[0], 2, ix-1, iy,iz) + am *  MRC_F3(E_m[0], 2, ix,iy,iz) + 
			 (ap*am)*   ( MRC_F3(u_m[0], BY, ix,iy,iz) -  MRC_F3(u_p[0], BY, ix-1,iy,iz) ));
    FLUX(flux, 0, _EZ, ix,iy,iz) = 
      (1.f/(ap - am)) * (     ap * MRC_F3(E_p[0], 1, ix-1,iy,iz) - am *  MRC_F3(E_m[0], 1, ix,iy,iz) + 
			 (ap*am) * ( MRC_F3(u_m[0], BZ, ix,iy,iz) - MRC_F3(u_p[0], BZ, ix-1,iy,iz) ));  
    
    assert(isfinite(bp));
    assert(isfinite(bm));
    //bp = 1e-3; bm = -1e-3;
    // flux of _EX,_EY,_EZ through the y faces    
    FLUX(flux, 1, _EX, ix,iy,iz) =
       (1.f/(bp - bm)) * (     bp  *  MRC_F3(E_p[1], 2, ix,iy-1,iz) - bm * MRC_F3(E_m[1], 2, ix,iy,iz) + 
		         (bp * bm) *  ( MRC_F3(u_m[1], BX, ix,iy,iz) -    MRC_F3(u_p[1], BX, ix,iy-1,iz) ));
    FLUX(flux, 1, _EY, ix,iy,iz) =
       (1.f/(bp - bm)) * ( (bp * bm)* ( MRC_F3(u_m[1], BY, ix,iy,iz) - MRC_F3(u_p[1], BY, ix,iy-1,iz) ));    
    FLUX(flux, 1, _EZ, ix,iy,iz) =
       (1.f/(bp - bm)) * (    - bp * MRC_F3(E_p[1], 0, ix,iy-1,iz) + bm *  MRC_F3(E_m[1], 0, ix,iy,iz) + 
		         (bp * bm) *  ( MRC_F3(u_m[1], BZ, ix,iy,iz) -     MRC_F3(u_p[1], BZ, ix,iy-1,iz) ));  
    
    assert(isfinite(cp));
    assert(isfinite(cm));
    //cp = 1e-3; cm = -1e-3;
    // flux of _EX,_EY,_EZ through the z faces
    FLUX(flux, 2, _EX, ix,iy,iz) = 
      (1.f/(cp - cm))*( - cp * MRC_F3(E_p[2], 1, ix,iy,iz-1) + cm * MRC_F3(E_m[2], 1, ix,iy,iz) + 
			 (cp * cm) * ( MRC_F3(u_m[2], BX, ix,iy,iz) - MRC_F3(u_p[2], BX, ix,iy,iz-1) ));
    FLUX(flux, 2, _EY, ix,iy,iz) =
      (1.f/(cp - cm)) *(  cp * MRC_F3(E_p[2], 0, ix,iy,iz-1) - cm * MRC_F3(E_m[2], 0, ix,iy,iz) +
		         (cp * cm) * ( MRC_F3(u_m[2], BY, ix,iy,iz) -    MRC_F3(u_p[2], BY, ix,iy,iz-1) ));
    FLUX(flux, 2, _EZ, ix,iy,iz) = 
      (1.f/(cp - cm))*(  (cp * cm) * ( MRC_F3(u_m[2], BZ, ix,iy,iz) - MRC_F3(u_p[2], BZ, ix,iy,iz-1) ));    

    for (int m = 0; m <= UU; m++) {
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
  
}
