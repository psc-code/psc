
#include "ggcm_mhd_step_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"

#include <mrc_domain.h>
#include <mrc_ddc.h>

#include <assert.h>
#include <math.h>

// ----------------------------------------------------------------------
// ggcm_mhd_get_fields

static struct mrc_f3 *
ggcm_mhd_get_fields(struct ggcm_mhd *mhd, const char *name, int nr_comps)
{ 
  struct mrc_f3 *f3 = mrc_domain_f3_create(mhd->domain, SW_2, NULL);
  mrc_f3_set_name(f3, name);
  mrc_f3_set_nr_comps(f3, nr_comps);
  mrc_f3_setup(f3);
  return f3;
}

// ======================================================================
// ggcm_mhd_step subclass "cweno"

// Define limiter: [1] van Leer(1977) , [2] minmod (Roe 1976) [3] moncen [4] genminmod

#define LMTR 1
#define sign(x) (( x > 0 ) - ( x < 0 ))
// KNP[0] or KT[1]? 
#define KT 0

enum {
  _EX = _B1Z + 1,
  _EY,
  _EZ,
  _JX,
  _JY,
  _JZ,
  __NR_FLDS,
};

static void
calc_neg_divg(struct mrc_f3 *rhs, int m, struct mrc_f3 *flux, struct mrc_crds *crds)
{
  mrc_f3_foreach(rhs, ix, iy, iz, 0, 0) {
    int ind[3] = { ix, iy, iz };

    MRC_F3(rhs, m, ix, iy, iz) = 0.;
    for(int i=0; i<3; i++) {
      int dind[3] = {0, 0, 0};
      dind[i] = 1;      

      float cw =(MRC_CRD(crds, i, ind[i]+1) - MRC_CRD(crds, i, ind[i]));
      MRC_F3(rhs, m, ix, iy, iz) -=( MRC_F3(flux, i, ix+dind[0],iy+dind[1],iz+dind[2]) -
      				     MRC_F3(flux, i, ix,iy,iz))/ cw;
    }
  } mrc_f3_foreach_end; 

  // make code stop if NaN or inf encountered, use 1 grid value in density as test   
  assert(isfinite(  MRC_F3(rhs, 0, 0, 0, 0)));
}

static void
calc_fct_rhs(struct ggcm_mhd *mhd, struct mrc_f3 *rhs, struct mrc_f3 *fld, struct mrc_f3 *flux_E[3])
{  
  // compute edge centered electric fields by interpolation (iv)  here tmp_fld are edge centered 
  // so that e.g.  MRC_F3(tmp_fld, 0, 0, 0, 0) is E_x 0,-1/2,-1/2  
  //         i.e.  MRC_F3(tmp_fld, 0, ix,iy,iz) is E_x ix,iy-1/2,iz-1/2   
  //           and MRC_F3(tmp_fld, 1, 0, 0, 0) is E_y -1/2,0,-1/2  etc etc 
                       

  struct mrc_ddc *ddc = mrc_domain_get_ddc(mhd->domain);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct mrc_f3 *E_ec = ggcm_mhd_get_fields(mhd, "E_ec", 3);

  struct mrc_f3 *fex = flux_E[0], *fey = flux_E[1], *fez = flux_E[2];

  //initialize cell edge center Electric field structure      
  mrc_f3_foreach(rhs, ix, iy,  iz, 1, 1) { 
    MRC_F3(E_ec, 0, ix, iy, iz) = 0.25*(-MRC_F3(fez, 1, ix, iy, iz) -
					 MRC_F3(fez, 1, ix, iy, iz-1) +
					 MRC_F3(fey, 2, ix, iy, iz) +
					 MRC_F3(fey, 2, ix, iy-1, iz));    
    //MRC_F3(fld,_EX, ix, iy, iz) = MRC_F3(E_ec, 0, ix, iy, iz);    
    MRC_F3(E_ec, 1, ix, iy, iz) = 0.25*( MRC_F3(fez, 0, ix, iy, iz) +
					 MRC_F3(fez, 0, ix, iy, iz-1) -
					 MRC_F3(fex, 2, ix, iy, iz) -
					 MRC_F3(fex, 2, ix-1, iy, iz));    
    //MRC_F3(fld,_EY, ix, iy, iz) = MRC_F3(E_ec, 1, ix, iy, iz);        
    MRC_F3(E_ec, 2, ix, iy, iz) = 0.25*(-MRC_F3(fey, 0, ix, iy, iz) -
					 MRC_F3(fey, 0, ix, iy-1, iz) +
					 MRC_F3(fex, 1, ix, iy, iz) +
					 MRC_F3(fex, 1, ix-1, iy, iz));    
    //MRC_F3(fld,_EZ, ix, iy, iz) = MRC_F3(E_ec, 2, ix, iy, iz);
  } mrc_f3_foreach_end;
  

  mrc_ddc_fill_ghosts(ddc, 0, 3, E_ec);
  mrc_f3_foreach(rhs, ix, iy,  iz, 1, 1) {
   
    B1X(rhs, ix, iy, iz) =  
      (-(( MRC_F3(E_ec, 2, ix, iy+1, iz) - MRC_F3(E_ec, 2, ix, iy, iz) ) /
	 (.5f*( MRC_CRDY(crds, iy+1) - MRC_CRDY(crds, iy-1) ))) 
       +(( MRC_F3(E_ec, 1, ix, iy, iz+1) - MRC_F3(E_ec, 1, ix, iy, iz) ) /
	 (.5f*( MRC_CRDZ(crds, iz+1) - MRC_CRDZ(crds, iz-1) )))); 
    
    B1Y(rhs, ix, iy, iz) = 
      (-(( MRC_F3(E_ec,  0, ix, iy, iz+1) - MRC_F3(E_ec, 0, ix, iy, iz) ) /
	 (.5f*( MRC_CRD(crds, 2, iz+1) - MRC_CRD(crds, 2, iz-1) )))
       +(( MRC_F3(E_ec, 2, ix+1, iy, iz) - MRC_F3(E_ec, 2, ix, iy, iz) ) /
	 (.5f*(MRC_CRD(crds, 0, ix+1) - MRC_CRD(crds, 0, ix-1) )))); 
    
    B1Z(rhs, ix, iy, iz) = 
      (-(( MRC_F3(E_ec, 1, ix+1, iy, iz) - MRC_F3(E_ec, 1, ix, iy, iz) ) /
	 (.5f*( MRC_CRD(crds, 0, ix+1) - MRC_CRD(crds, 0, ix-1) )))  
       +(( MRC_F3(E_ec, 0, ix, iy+1, iz) - MRC_F3(E_ec, 0, ix, iy, iz) ) /
	 (.5f*( MRC_CRD(crds, 1, iy+1) - MRC_CRD(crds, 1, iy-1) )))); 
   
  } mrc_f3_foreach_end;
    
  mrc_f3_destroy(E_ec); 
  //  MHERE;
}

// ----------------------------------------------------------------------
// calc_fluxes_per_faces
//
// this calculates fluxes on the face i using reconstructed variables (fld) that are
// given on the respective face

static void
calc_fluxes_per_face(struct mrc_f3 **flux, struct ggcm_mhd *mhd, struct mrc_f3 *fld, int i)
{
  float mpermi = 1.f;
  float gamma = mhd->par.gamm;
  float d_i = mhd->par.d_i;

  mrc_f3_foreach(fld, ix, iy, iz, 1, 1) {
    float rhoi = 1.f / MRC_F3(fld, _RR1, ix,iy,iz);
      
    float BB = (0.5f) *mpermi * (sqr(MRC_F3(fld, _B1X, ix,iy,iz)) +
					  sqr(MRC_F3(fld, _B1Y, ix,iy,iz)) +
					  sqr(MRC_F3(fld, _B1Z, ix,iy,iz)));
    
    float mB = (MRC_F3(fld, _B1X, ix,iy,iz)*MRC_F3(fld, _RV1X, ix,iy,iz)) + 
               (MRC_F3(fld, _B1Y, ix,iy,iz)*MRC_F3(fld, _RV1Y, ix,iy,iz)) + 
               (MRC_F3(fld, _B1Z, ix,iy,iz)*MRC_F3(fld, _RV1Z, ix,iy,iz)) ; 

    
    float JB = -(MRC_F3(fld, _B1X, ix,iy,iz)*MRC_F3(fld, _JX, ix,iy,iz))  
                -(MRC_F3(fld, _B1Y, ix,iy,iz)*MRC_F3(fld, _JY, ix,iy,iz))  
                -(MRC_F3(fld, _B1Z, ix,iy,iz)*MRC_F3(fld, _JZ, ix,iy,iz)) ; 



    float pp = (gamma - 1.f) *
      (MRC_F3(fld, _UU1, ix,iy,iz) - .5f * rhoi * (sqr(MRC_F3(fld, _RV1X, ix,iy,iz)) +
						 sqr(MRC_F3(fld, _RV1Y, ix,iy,iz)) +
						 sqr(MRC_F3(fld, _RV1Z, ix,iy,iz)))- 
       (.5f * mpermi * (sqr(MRC_F3(fld, _B1X, ix,iy,iz)) +
				 sqr(MRC_F3(fld, _B1Y, ix,iy,iz)) +
				 sqr(MRC_F3(fld, _B1Z, ix,iy,iz)))));
    
    //float BB = (0.5f) *mpermi * (sqr(MRC_F3(fld, _B1X+i, ix,iy,iz)));
    //float pp = (mhd->T)*MRC_F3(fld, _RR1, ix,iy,iz);
    // mass consv. 
    MRC_F3(flux[_RR1], i, ix,iy,iz) = MRC_F3(fld, _RV1X+i, ix,iy,iz);
      //MRC_F3(fld, _RV1X+i, ix,iy,iz);
    
    // momentum eq. 
    for (int j = 0; j < 3; j++) {
      
      //BB = (0.5f) *mpermi * (sqr(MRC_F3(fld, _B1X+j, ix,iy,iz)));
      MRC_F3(flux[_RV1X+i], j, ix,iy,iz) = 
	rhoi * MRC_F3(fld, _RV1X+j, ix,iy,iz) * MRC_F3(fld, _RV1X+i, ix,iy,iz) +
	((j == i) ? pp : 0.) + 
	((j == i) ? BB : 0.) - mpermi * (MRC_F3(fld, _B1X+i, ix,iy,iz) * MRC_F3(fld, _B1X+j, ix,iy,iz));
      //mpermi * MRC_F3(fld, _B1X+j, ix,iy,iz) * MRC_F3(fld, _B1X+i, ix,iy,iz);
       
    }
    
    // energy eq. 
    MRC_F3(flux[_UU1], i, ix,iy,iz) =
      ( ((MRC_F3(fld, _UU1, ix,iy,iz) + pp + BB)*MRC_F3(fld, _RV1X+i, ix,iy,iz))-
	(mpermi * mB * MRC_F3(fld, _B1X+i, ix,iy,iz)) + 
	(d_i * ( -0.5*MRC_F3(fld, _JX+i, ix,iy,iz)*BB - MRC_F3(fld, _B1X+i, ix,iy,iz)*JB)) ) * rhoi;

    /*
    MRC_F3(flux[_UU1], i, ix,iy,iz) =
      (MRC_F3(fld, _UU1, ix,iy,iz) + pp+BB) * rhoi * MRC_F3(fld, _RV1X+i, ix,iy,iz) +
      (mpermi * rhoi* MRC_F3(fld, _RV1X+i, ix,iy,iz) * MRC_F3(fld, _B1X+i, ix,iy,iz)) * 
      (MRC_F3(fld, _B1X+i, ix,iy,iz));
    */


  } mrc_f3_foreach_end;
 
}

// ----------------------------------------------------------------------
// calc_cweno_fluxes
//
// calculates CWENO fluxes on faces in flux, flux_E, from the original state
// vector u (which is cell centered / on the Yee grid)

static void
calc_cweno_fluxes(struct mrc_f3 **flux, struct mrc_f3 *flux_E[3], struct ggcm_mhd *mhd,
		  struct mrc_f3 *u, struct mrc_crds *crds)
{
  struct mrc_ddc *ddc = mrc_domain_get_ddc(mhd->domain);

  float mpermi = 1.f;
  float gamma = mhd->par.gamm;
  float d_i = mhd->par.d_i;
  float eta = mhd->par.diffco;

  mrc_ddc_fill_ghosts(ddc, 0, _B1Z + 1, u);

  // initialize deltas for reconstruction
  struct mrc_f3 *u_delta[3];
  for (int f = 0; f < 3; f++) {
    u_delta[f] = ggcm_mhd_get_fields(mhd, "u_delta", _B1Z + 1);
  }

#if LMTR == 1
  // TVD slope limiter  (van Leer 1977) geometric mean
  // (dxu_)ijk = max{ (u_i+1jk-u_ijk)*(u_ijk-u_i-1jk)  / (u_i+1jk - u_i-1jk)}
  mrc_f3_foreach(u, ix,iy,iz, 1, 1) {
    for (int m = 0; m < _B1Z+1; m++) {
      MRC_F3(u_delta[0], m, ix,iy,iz) = 
	0.5*fmaxf((MRC_F3(u, m, ix+1,iy,iz) - MRC_F3(u, m, ix  ,iy,iz)) *
	      (MRC_F3(u, m, ix  ,iy,iz) - MRC_F3(u, m, ix-1,iy,iz)) , 0.f) /
	(MRC_F3(u, m, ix+1,iy,iz) - MRC_F3(u, m, ix-1,iy,iz));
      MRC_F3(u_delta[1], m, ix,iy,iz) = 
	0.5*fmaxf((MRC_F3(u, m, ix,iy+1,iz) - MRC_F3(u, m, ix,iy  ,iz)) *
	      (MRC_F3(u, m, ix,iy  ,iz) - MRC_F3(u, m, ix,iy-1,iz)), 0.f) / 
	(MRC_F3(u, m, ix,iy+1,iz) - MRC_F3(u, m, ix,iy-1,iz));
      MRC_F3(u_delta[2], m, ix,iy,iz) = 
	0.5*fmaxf((MRC_F3(u, m, ix,iy,iz+1) - MRC_F3(u, m, ix,iy,iz  )) *
	      (MRC_F3(u, m, ix,iy,iz  ) - MRC_F3(u, m, ix,iy,iz-1)), 0.f) /
	(MRC_F3(u, m, ix,iy,iz+1) - MRC_F3(u, m, ix,iy,iz-1));
      // FIXME, need to make sure NaN -> 0
      if (!isfinite(MRC_F3(u_delta[0], m, ix,iy,iz))) MRC_F3(u_delta[0], m, ix,iy,iz) = 0.f;
      if (!isfinite(MRC_F3(u_delta[1], m, ix,iy,iz))) MRC_F3(u_delta[1], m, ix,iy,iz) = 0.f;
      if (!isfinite(MRC_F3(u_delta[2], m, ix,iy,iz))) MRC_F3(u_delta[2], m, ix,iy,iz) = 0.f;
    }   
    } mrc_f3_foreach_end;
#endif

#if LMTR == 2
   // MinMod limiter (Roe, 1986)
  // return the smallest if all are positive 
  // return the largest if all are negative
  // return zero if they are not all postivie or all negative  
  // minmod(a,b) = 0.5[sign(a)+sign(b)]*min(|a|,|b|) eq (41) 
  // 
  int dind[3] = {0, 0, 0};  
  for (int i = 0; i < 3; i++) {    
    dind[i] = 1; 
    for (int m = 0; m < _B1Z+1; m++) {    
      mrc_f3_foreach(u, ix,iy,iz, 1, 1) {
	int ind[3] = { ix, iy, iz };           
	float rl =  ( MRC_F3(u, m, ix+dind[0], iy+dind[1], iz+dind[2])  - MRC_F3(u, m, ix, iy, iz));
	float rr =  ( MRC_F3(u, m, ix, iy, iz)  - MRC_F3(u, m, ix-dind[0], iy-dind[1], iz-dind[2]));
	if ( rl*rr > 0 ) {
	  if ( fabsf(rr)<fabsf(rl) ) {
	    MRC_F3(u_delta[i], m, ix,iy,iz) = 0.5*rr ;            
	    MHERE;
	  } else { 
	    MRC_F3(u_delta[i], m, ix,iy,iz) = 0.5*rl ;    
	    MHERE;
	  }
	} else { 
	  MRC_F3(u_delta[i], m, ix,iy,iz) = 0.0; 
	  MHERE;
	}
      } mrc_f3_foreach_end;      
    }
    dind[i] = 0;
  }
#endif
 
#if LMTR == 3
   // MONCEN limiter (Van Leer)
  int dind[3] = {0, 0, 0};
  
  for (int i = 0; i < 3; i++) {    
    dind[i] = 1; 
    for (int m = 0; m < _B1Z+1; m++) {    
      mrc_f3_foreach(u, ix,iy,iz, 1, 1) {
	int ind[3] = { ix, iy, iz };           
	float rl =  2.0*( MRC_F3(u, m, ix+dind[0], iy+dind[1], iz+dind[2])  - MRC_F3(u, m, ix, iy, iz));
	float rr =  2.0*( MRC_F3(u, m, ix, iy, iz)  - MRC_F3(u, m, ix-dind[0], iy-dind[1], iz-dind[2]));
	float cen = 0.5*( MRC_F3(u, m, ix+dind[0], iy+dind[1], iz+dind[2])  - 
		      MRC_F3(u, m, ix-dind[0], iy-dind[1], iz-dind[2]) ) ; 
	
	if ( (cen*rl*rr) > 0) {
	  MRC_F3(u_delta[i], m, ix,iy,iz) = 0.5*sign(cen)*fmin(fabsf(cen),fmin(fabsf(rl),fabsf(rr))) ; 
	} else {
	  MRC_F3(u_delta[i], m, ix,iy,iz) = 0.0 ; 
	}
	
	//MRC_F3(u_delta[i], m, ix,iy,iz) = 0.0 ; 
      if (!isfinite(MRC_F3(u_delta[i], m, ix,iy,iz))) MRC_F3(u_delta[i], m, ix,iy,iz) = 0.f;
      } mrc_f3_foreach_end;      
    }
    dind[i] = 0;
  }

#endif
 
#if LMTR == 4
   // generalised minmod limiter with parameter(Van Leer 1979)
  int dind[3] = {0, 0, 0};
  float theta = 1.0; 
  for (int i = 0; i < 3; i++) {    
    dind[i] = 1; 
    for (int m = 0; m < _B1Z+1; m++) {    
      mrc_f3_foreach(u, ix,iy,iz, 1, 1) {
	int ind[3] = { ix, iy, iz };           
	float rl =  theta*( MRC_F3(u, m, ix+dind[0], iy+dind[1], iz+dind[2])  - MRC_F3(u, m, ix, iy, iz));
	float rr =  theta*( MRC_F3(u, m, ix, iy, iz)  - MRC_F3(u, m, ix-dind[0], iy-dind[1], iz-dind[2]));
	float cen = 0.5*( MRC_F3(u, m, ix+dind[0], iy+dind[1], iz+dind[2])  - 
		      MRC_F3(u, m, ix-dind[0], iy-dind[1], iz-dind[2]) ) ; 
	
	if ( (cen*rl*rr) > 0) {
	  MRC_F3(u_delta[i], m, ix,iy,iz) = 0.5*sign(cen)*fmin(fabsf(cen),fmin(fabsf(rl),fabsf(rr))) ; 
	} else {
	  MRC_F3(u_delta[i], m, ix,iy,iz) = 0.0 ; 
	}
	
	//MRC_F3(u_delta[i], m, ix,iy,iz) = 0.0 ; 
      if (!isfinite(MRC_F3(u_delta[i], m, ix,iy,iz))) MRC_F3(u_delta[i], m, ix,iy,iz) = 0.f;
      } mrc_f3_foreach_end;      
    }
    dind[i] = 0;
  }

#endif

  mrc_ddc_fill_ghosts(ddc, 0, _B1Z + 1, u_delta[0]);
  mrc_ddc_fill_ghosts(ddc, 0, _B1Z + 1, u_delta[1]);
  mrc_ddc_fill_ghosts(ddc, 0, _B1Z + 1, u_delta[2]);
                       
  //initialize cell surface center variables			    
  struct mrc_f3 *u_p[3], *u_m[3];
  for (int f = 0; f < 3; f++) {
    u_p[f] = ggcm_mhd_get_fields(mhd, "u_p", _JZ + 1);
    u_m[f] = ggcm_mhd_get_fields(mhd, "u_m", _JZ + 1);
  }

  // Reonstruction    UijkE  = u_ijk + (dxu_)ijk    UijkW = u_ijk - (dxu_)ijk
  mrc_f3_foreach(u, ix,iy,iz, 1, 1) {
    for (int m = 0; m < _UU1+1; m++) {
      // defined like this, both u_p and u_m are coplaner when indices are the same

      MRC_F3(u_p[0], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) + MRC_F3(u_delta[0], m, ix,iy,iz);
      MRC_F3(u_p[1], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) + MRC_F3(u_delta[1], m, ix,iy,iz);
      MRC_F3(u_p[2], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) + MRC_F3(u_delta[2], m, ix,iy,iz);
      MRC_F3(u_m[0], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) - MRC_F3(u_delta[0], m, ix,iy,iz);
      MRC_F3(u_m[1], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) - MRC_F3(u_delta[1], m, ix,iy,iz);
      MRC_F3(u_m[2], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) - MRC_F3(u_delta[2], m, ix,iy,iz);
    }
  } mrc_f3_foreach_end;
  
  // calculation of cell surface center averages for B 
  // note that B is staggered here.. so that index 1 for B is 1-1/2 i.e. cell surface for V
  mrc_f3_foreach(u, ix,iy,iz, 1, 1) {
    for (int i = 0; i < 3; i++) {

      int cycl[5]={0,1,2,0,1};
      int dind[3]={0,0,0};
      int ip1= cycl[i+1];
      int ip2= cycl[i+2];

      // reconstructio for B compoents. e.g. i E-W location, Bx is already in correct locatio 
      // but transverse components are first reconstructed to cell edges then averaged to cell surfaces
   

      // _p  i-> 0:E 1:   
      dind[i]=1; 
      MRC_F3(u_p[i], _B1X+i, ix,iy,iz) = B1XYZ(u, i, ix+dind[0],iy+dind[1],iz+dind[2]);
      // +1 --> Bxi+1/2y,z because staggered grid 
      dind[i]=0;

      dind[ip1]=1;
      MRC_F3(u_p[i], _B1X+ip1, ix,iy,iz) =
	(0.5*(B1XYZ(u, ip1, ix+dind[0],iy+dind[1],iz+dind[2]) +
	      MRC_F3(u_delta[i], _B1X+ip1, ix+dind[0],iy+dind[1],iz+dind[2]) +
	      B1XYZ(u, ip1, ix,iy,iz) + MRC_F3(u_delta[i], _B1X+ip1, ix,iy,iz)));      
      dind[ip1]=0;
       
      dind[ip2]=1;
      MRC_F3(u_p[i], _B1X+ip2, ix,iy,iz) =
	(0.5*(B1XYZ(u, ip2, ix+dind[0],iy+dind[1],iz+dind[2] ) + 
	      MRC_F3(u_delta[i], _B1X+ip2, ix+dind[0],iy+dind[1],iz+dind[2]) +
	      B1XYZ(u, ip2, ix,iy,iz) + MRC_F3(u_delta[i], _B1X+ip2, ix,iy,iz)));
      dind[ip2]=0;
      
      // _m
      dind[i]=0; 
      MRC_F3(u_m[i], _B1X+i, ix,iy,iz) = B1XYZ(u, i, ix,iy,iz);
      //  +0 --> Bxi-1/2y,z because staggered grid
      
      dind[ip1]=1;
      MRC_F3(u_m[i], _B1X+ip1, ix,iy,iz) = 
	(0.5*(B1XYZ(u, ip1, ix+dind[0],iy+dind[1],iz+dind[2]) -
	      MRC_F3(u_delta[i], _B1X+ip1, ix+dind[0],iy+dind[1],iz+dind[2]) +
	      B1XYZ(u, ip1, ix,iy,iz) - MRC_F3(u_delta[i], _B1X+ip1, ix,iy,iz)));
      dind[ip1]=0;

      dind[ip2]=1;
      MRC_F3(u_m[i], _B1X+ip2, ix,iy,iz) =
	(0.5*(B1XYZ(u, ip2, ix+dind[0],iy+dind[1],iz+dind[2] ) -
	      MRC_F3(u_delta[i], _B1X+ip2, ix+dind[0],iy+dind[1],iz+dind[2]) +
	      B1XYZ(u, ip2, ix,iy,iz) - MRC_F3(u_delta[i], _B1X+ip2, ix,iy,iz)));
      dind[ip2]=0;
    }
  } mrc_f3_foreach_end;
  
  float *bdx3 = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD3);
  float *bdy3 = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD3);
  float *bdz3 = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD3);

  for (int i = 0; i < 3; i++) {    
    mrc_f3_foreach(u, ix,iy,iz, 2, 2) {	
      // _p
      MRC_F3(u_p[i],_JX, ix,iy,iz) = 
	0.5*((MRC_F3(u,_B1Z, ix,iy+1,iz) - MRC_F3(u,_B1Z, ix,iy-1,iz)) * bdy3[iy] - 
	     (MRC_F3(u,_B1Y, ix,iy,iz+1) - MRC_F3(u,_B1Y, ix,iy,iz-1)) * bdz3[iz]);     

      MRC_F3(u_p[i],_JY, ix,iy,iz) =
	0.5*((MRC_F3(u,_B1X, ix,iy,iz+1) - MRC_F3(u,_B1X, ix,iy,iz-1)) * bdz3[iz] -
	     (MRC_F3(u,_B1Z, ix+1,iy,iz) - MRC_F3(u,_B1Z, ix-1,iy,iz)) * bdx3[ix]);       

      MRC_F3(u_p[i],_JZ, ix,iy,iz) = 
	0.5*((MRC_F3(u,_B1Y, ix+1,iy,iz) - MRC_F3(u,_B1Y, ix-1,iy,iz)) * bdx3[ix] - 
	     (MRC_F3(u,_B1X, ix,iy+1,iz) - MRC_F3(u,_B1X, ix,iy-1,iz)) * bdy3[iy]); 

     
      // _m 
      MRC_F3(u_m[i], _JX, ix,iy,iz) = MRC_F3(u_p[i], _JX, ix-1,iy,iz) ; 
      MRC_F3(u_m[i], _JY, ix,iy,iz) = MRC_F3(u_p[i], _JY, ix,iy-1,iz) ;
      MRC_F3(u_m[i], _JZ, ix,iy,iz) = MRC_F3(u_p[i], _JZ, ix,iy,iz-1) ;      
      
    } mrc_f3_foreach_end;
  }
  
  // Calculation of cell surface center  averages for E using v and B just calculated 
  //  E^N(SWETB) = -V^N(SWETB) X B^N(SWETB)      + eta*J  [  + di J x B ] 
  mrc_f3_foreach(u, ix, iy, iz, 1, 1) {  
    for (int i = 0; i < 3; i++) {     
      // _p 
      MRC_F3(u_p[i], _EX, ix,iy,iz) = 
	-((MRC_F3(u_p[i], _RV1Y, ix,iy,iz) * MRC_F3(u_p[i], _B1Z, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz))- 	
	  (MRC_F3(u_p[i], _B1Y, ix,iy,iz) * MRC_F3(u_p[i], _RV1Z, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz)))+
	  eta *  MRC_F3(u_p[i],_JX, ix,iy,iz)
	+ d_i * ((MRC_F3(u_p[i], _JY, ix,iy,iz) * MRC_F3(u_p[i], _B1Z, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz))- 	
		 (MRC_F3(u_p[i], _B1Y, ix,iy,iz) * MRC_F3(u_p[i], _JZ, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz)));  
  
      MRC_F3(u_p[i], _EY, ix,iy,iz) = 
	((MRC_F3(u_p[i], _RV1X, ix,iy,iz) * MRC_F3(u_p[i], _B1Z, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz))- 	
	 (MRC_F3(u_p[i], _B1X, ix,iy,iz) * MRC_F3(u_p[i], _RV1Z, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz)))+
	 eta * MRC_F3(u_p[i],_JY, ix,iy,iz)
	- d_i * ((MRC_F3(u_p[i], _JX, ix,iy,iz) * MRC_F3(u_p[i], _B1Z, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz))- 	
		 (MRC_F3(u_p[i], _B1X, ix,iy,iz) * MRC_F3(u_p[i], _JZ, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz)));
      
      MRC_F3(u_p[i], _EZ, ix,iy,iz) = 
	-((MRC_F3(u_p[i], _RV1X, ix,iy,iz) * MRC_F3(u_p[i], _B1Y, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz))- 	
	  (MRC_F3(u_p[i], _B1X, ix,iy,iz) * MRC_F3(u_p[i], _RV1Y, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz)))+
	  eta * MRC_F3(u_p[i],_JZ, ix,iy,iz)
        + d_i * ((MRC_F3(u_p[i], _JX, ix,iy,iz) * MRC_F3(u_p[i], _B1Y, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz))- 	
		     (MRC_F3(u_p[i], _B1X, ix,iy,iz) * MRC_F3(u_p[i], _JY, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz)));

	
      // _m 
      MRC_F3(u_m[i], _EX, ix,iy,iz) =  
	-((MRC_F3(u_m[i], _RV1Y, ix,iy,iz) * MRC_F3(u_m[i], _B1Z, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz))- 	
	  (MRC_F3(u_m[i], _B1Y, ix,iy,iz) * MRC_F3(u_m[i], _RV1Z, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz)))+
   	  eta * MRC_F3(u_m[i],_JX, ix,iy,iz)
	+ d_i * ((MRC_F3(u_m[i], _JY, ix,iy,iz) * MRC_F3(u_m[i], _B1Z, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz))- 	
		 (MRC_F3(u_m[i], _B1Y, ix,iy,iz) * MRC_F3(u_m[i], _JZ, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz)));
      
      MRC_F3(u_m[i], _EY, ix,iy,iz) = 
	((MRC_F3(u_m[i], _RV1X, ix,iy,iz) * MRC_F3(u_m[i], _B1Z, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz))- 	
	 (MRC_F3(u_m[i], _B1X, ix,iy,iz) * MRC_F3(u_m[i], _RV1Z, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz)))+
	 eta * MRC_F3(u_m[i],_JY, ix,iy,iz)
	- d_i * ((MRC_F3(u_m[i], _JX, ix,iy,iz) * MRC_F3(u_m[i], _B1Z, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz))- 	
		 (MRC_F3(u_m[i], _B1X, ix,iy,iz) * MRC_F3(u_m[i], _JZ, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz)));
      
      MRC_F3(u_m[i], _EZ, ix,iy,iz) =  
	-((MRC_F3(u_m[i], _RV1X, ix,iy,iz) * MRC_F3(u_m[i], _B1Y, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz))- 	
	  (MRC_F3(u_m[i], _B1X, ix,iy,iz) * MRC_F3(u_m[i], _RV1Y, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz)))+
	  eta * MRC_F3(u_m[i],_JZ, ix,iy,iz) 
        + d_i * ((MRC_F3(u_m[i], _JX, ix,iy,iz) * MRC_F3(u_m[i], _B1Y, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz))- 	
		 (MRC_F3(u_m[i], _B1X, ix,iy,iz) * MRC_F3(u_m[i], _JY, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz)));   
    }    
  } mrc_f3_foreach_end;
  
  //mrc_f3_destroy(j);

  // initialize flux functions
  struct mrc_f3 *flux_p[_B1Z + 1], *flux_m[_B1Z + 1];
  for (int m = 0; m < _B1Z + 1; m++) {
    flux_p[m] = ggcm_mhd_get_fields(mhd, "flux_p", 3);
    flux_m[m] = ggcm_mhd_get_fields(mhd, "flux_m", 3);
  }
  
  // calculate fluxes per face (small f's) using reconstructed 
  // variables U^N(SWETB) and B^N(SWETB) = (Bx,By,Bz)^N(SEWTB)
  for (int f = 0; f < 3; f++) {
    calc_fluxes_per_face(flux_p, mhd, u_p[f], f);
    calc_fluxes_per_face(flux_m, mhd, u_m[f], f);
  }
  
  //  float cs;
  float csp; 
  float csm; 
  float tmp;
  float cAp;
  float cAm;
  float cfp; 
  float cfm; 
  float ap; 
  float am; 
  float bp; 
  float bm; 
  float cp;
  float cm;
  float ppp;
  float ppm; 
  //  float rhoi;
  float rhoip;
  float rhoim; 

  mrc_f3_foreach(u, ix,iy,iz, 1, 1) {
    
   // Coeffiecents ap and am   
 
    rhoip = 1.f / MRC_F3(u_p[0], _RR1, ix-1,iy,iz);	
    rhoim = 1.f / MRC_F3(u_m[0], _RR1, ix,iy,iz);	
    ppp = (gamma - 1.f) *
      (MRC_F3(u_p[0], _UU1, ix-1,iy,iz) - .5f * rhoip * (sqr(MRC_F3(u_p[0], _RV1X, ix-1,iy,iz)) +
						     sqr(MRC_F3(u_p[0], _RV1Y, ix-1,iy,iz)) +
						     sqr(MRC_F3(u_p[0], _RV1Z, ix-1,iy,iz)))- 
       .5f * (sqr(MRC_F3(u_p[0], _B1X, ix-1,iy,iz)) +
	      sqr(MRC_F3(u_p[0], _B1Y, ix-1,iy,iz)) +
	      sqr(MRC_F3(u_p[0], _B1Z, ix-1,iy,iz))));
       
    ppm = (gamma - 1.f) *
      (MRC_F3(u_m[0], _UU1, ix,iy,iz) - .5f * rhoim * (sqr(MRC_F3(u_m[0], _RV1X, ix,iy,iz)) +
						       sqr(MRC_F3(u_m[0], _RV1Y, ix,iy,iz)) +
						       sqr(MRC_F3(u_m[0], _RV1Z, ix,iy,iz)))-       
       .5f * (sqr(MRC_F3(u_m[0], _B1X, ix,iy,iz)) +
	      sqr(MRC_F3(u_m[0], _B1Y, ix,iy,iz)) +
	      sqr(MRC_F3(u_m[0], _B1Z, ix,iy,iz))));
    
    csp = sqrtf((gamma * ppp) / (MRC_F3(u_p[0], _RR1, ix-1,iy,iz)));
    csm = sqrtf((gamma * ppm) / (MRC_F3(u_m[0], _RR1, ix,iy,iz)));   
       
    cAp = sqrtf( (sqr(MRC_F3(u_p[0], _B1X, ix-1,iy,iz))+
		  sqr(MRC_F3(u_p[0], _B1Y, ix-1,iy,iz))+
		  sqr(MRC_F3(u_p[0], _B1Z, ix-1,iy,iz)))/ MRC_F3( u_p[0], _RR1, ix-1,iy,iz) );

    cAm = sqrtf( (sqr(MRC_F3(u_m[0], _B1X, ix,iy,iz))+
		  sqr(MRC_F3(u_m[0], _B1Y, ix,iy,iz))+
		  sqr(MRC_F3(u_m[0], _B1Z, ix,iy,iz)))/ MRC_F3(u_m[0], _RR1, ix,iy,iz) );
    //cAp = sqrtf( (sqr(MRC_F3(u_p[0], _B1X, ix-1,iy,iz)))/ MRC_F3( u_p[0], _RR1, ix-1,iy,iz) );
    //cAm = sqrtf( (sqr(MRC_F3(u_m[0], _B1X, ix,iy,iz)))/ MRC_F3(u_m[0], _RR1, ix,iy,iz) );		
 						      
    tmp = sqr(csp) + sqr(cAp);
    cfp = sqrtf( 0.5 * ( tmp + sqrtf( sqr( sqr(cAp) + sqr(csp) )
				     - (4. * mpermi * sqr(csp * MRC_F3(u, _B1X, ix-1,iy,iz)) /  
					MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) )) );      
    //cfp = sqrt(tmp);//sqrt(2.0)*sqrt(sqr(csp) + sqr(cAp));
    tmp = sqr(csm) + sqr(cAm);
    cfm = sqrtf( 0.5 * ( tmp + sqrtf( sqr( sqr(cAm) + sqr(csm) )
				     - (4.* mpermi * sqr(csm * MRC_F3(u, _B1X, ix,iy,iz)) /  
					MRC_F3(u_m[0], _RR1, ix,iy,iz))  )) );

    //cfm = sqrt(tmp);//sqrt(2.0)*sqrt(sqr(csm) + sqr(cAm));
    ap = fmaxf(fmaxf((MRC_F3(u_p[0], _RV1X, ix-1,iy,iz) / MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) + cfp,
		     (MRC_F3(u_m[0], _RV1X, ix,iy,iz) / MRC_F3(u_m[0], _RR1, ix,iy,iz)) + cfm), 0.f);

    am = fminf(fminf( (MRC_F3(u_p[0], _RV1X, ix-1,iy,iz) / MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) - cfp,
		     (MRC_F3(u_m[0], _RV1X, ix,iy,iz) / MRC_F3(u_m[0], _RR1,ix,iy,iz)) - cfm), 0.f);
    
   
    
#if KT == 1
    ap = fmaxf(ap,-am);
    am=-ap;
#endif    


    // Coeffiecents bp and bm 
    rhoip = 1.f / MRC_F3(u_p[1], _RR1, ix,iy-1,iz);
    rhoim = 1.f / MRC_F3(u_m[1], _RR1, ix,iy,iz); 
    
    ppp = (gamma - 1.f) *
      (MRC_F3(u_p[1], _UU1, ix,iy-1,iz) - .5f * rhoip * (sqr(MRC_F3(u_p[1], _RV1X, ix,iy-1,iz)) +
						    sqr(MRC_F3(u_p[1], _RV1Y, ix,iy-1,iz)) +
						     sqr(MRC_F3(u_p[1], _RV1Z, ix,iy-1,iz)))-       
                                             .5f * (sqr(MRC_F3(u_p[1], _B1X, ix,iy-1,iz)) +
						    sqr(MRC_F3(u_p[1], _B1Y, ix,iy-1,iz)) +
	 					    sqr(MRC_F3(u_p[1], _B1Z, ix,iy-1,iz))));    
   
    ppm = (gamma - 1.f) *
      (MRC_F3(u_m[1], _UU1, ix,iy,iz) - .5f * rhoim * (sqr(MRC_F3(u_m[1], _RV1X, ix,iy,iz)) +
						       sqr(MRC_F3(u_m[1], _RV1Y, ix,iy,iz)) +
						       sqr(MRC_F3(u_m[1], _RV1Z, ix,iy,iz)))-	
       .5f * (sqr(MRC_F3(u_m[1], _B1X, ix,iy,iz)) +
	      sqr(MRC_F3(u_m[1], _B1Y, ix,iy,iz)) +
	      sqr(MRC_F3(u_m[1], _B1Z, ix,iy,iz))));
    
    csp = sqrtf((gamma * ppp) / (MRC_F3(u_p[1], _RR1, ix,iy-1,iz)));
    csm = sqrtf((gamma * ppm) / (MRC_F3(u_m[1], _RR1, ix,iy,iz)));

    cAp = sqrtf( (sqr(MRC_F3(u, _B1X, ix,iy-1,iz))+
		  sqr(MRC_F3(u, _B1Y, ix,iy-1,iz))+
		  sqr(MRC_F3(u, _B1Z, ix,iy-1,iz)))/ MRC_F3(u_p[1], _RR1, ix,iy-1,iz) );

    cAm = sqrtf( (sqr(MRC_F3(u, _B1X, ix,iy,iz))+
		  sqr(MRC_F3(u, _B1Y, ix,iy,iz))+
		  sqr(MRC_F3(u, _B1Z, ix,iy,iz)))/ MRC_F3(u_m[1], _RR1, ix,iy,iz) );
		   	       
    			  	  
    tmp = sqr(csp) + sqr(cAp);
    cfp = sqrtf( 0.5 * (tmp + sqrtf( sqr( sqr(cAp) + sqr(csp) ) - 
				     (4.f * mpermi * sqr(csp * MRC_F3(u, _B1Y, ix,iy-1,iz)) /  
				      MRC_F3(u_p[1], _RR1, ix,iy-1,iz))) ));      
    //cfp = sqrt(tmp);//sqrt(2.0)*sqrt( sqr(csp) + sqr(cAp));
    tmp = sqr(csm) + sqr(cAm);
    cfm = sqrtf( 0.5 * (tmp + sqrtf( sqr( sqr(cAm) + sqr(csm) ) -  
				     (4.f * mpermi * sqr(csm * MRC_F3(u, _B1Y, ix,iy,iz)) /  
				      MRC_F3(u_m[1], _RR1, ix,iy,iz))) ));
    //cfm = sqrt(tmp);////sqrt(2.0)*sqrt( sqr(csm) + sqr(cAm));
   
    bp = fmaxf(fmaxf( (MRC_F3(u_p[1], _RV1Y, ix,iy-1,iz) / MRC_F3(u_p[1], _RR1, ix,iy-1,iz)) + cfp,
		      (MRC_F3(u_m[1], _RV1Y, ix,iy,iz) / MRC_F3(u_m[1], _RR1, ix,iy,iz)) + cfm), 0.f);
    bm = fminf(fminf( (MRC_F3(u_p[1], _RV1Y, ix,iy-1,iz) / MRC_F3(u_p[1], _RR1, ix,iy-1,iz)) - cfp,
		      (MRC_F3(u_m[1], _RV1Y, ix,iy,iz) / MRC_F3(u_m[1], _RR1, ix,iy,iz)) - cfm), 0.f);

  
    
#if KT == 1
    bp = fmaxf(bp,-bm);
    bm=-bp;
#endif    

   // Coeffiecents cp and cm 

    rhoip = 1.f / MRC_F3(u_p[2], _RR1, ix,iy,iz-1);
    rhoim = 1.f / MRC_F3(u_m[2], _RR1, ix,iy,iz); 


    ppp = (gamma - 1.f) *
      (MRC_F3(u_p[2], _UU1, ix,iy,iz-1) - .5f * rhoip * (sqr(MRC_F3(u_p[2], _RV1X, ix,iy,iz-1)) +
						     sqr(MRC_F3(u_p[2], _RV1Y, ix,iy,iz-1)) +
						     sqr(MRC_F3(u_p[2], _RV1Z, ix,iy,iz-1)))-      
       .5f * (sqr(MRC_F3(u_p[2], _B1X, ix,iy,iz-1)) +
	      sqr(MRC_F3(u_p[2], _B1Y, ix,iy,iz-1)) +
	      sqr(MRC_F3(u_p[2], _B1Z, ix,iy,iz-1))));   
    
    ppm = (gamma - 1.f) *
      (MRC_F3(u_m[2], _UU1, ix,iy,iz) - .5f * rhoim * (sqr(MRC_F3(u_m[2], _RV1X, ix,iy,iz)) +
						       sqr(MRC_F3(u_m[2], _RV1Y, ix,iy,iz)) +
						       sqr(MRC_F3(u_m[2], _RV1Z, ix,iy,iz)))-        
       .5f * (sqr(MRC_F3(u_m[2], _B1X, ix,iy,iz)) +
	      sqr(MRC_F3(u_m[2], _B1Y, ix,iy,iz)) +
	      sqr(MRC_F3(u_m[2], _B1Z, ix,iy,iz))));
       

  
    csp = sqrtf((gamma * ppp) / (MRC_F3(u_p[2], _RR1, ix,iy,iz-1)));
    csm = sqrtf((gamma * ppm) / (MRC_F3(u_m[2], _RR1, ix,iy,iz)));
   			
    cAp = sqrtf( (sqr(MRC_F3(u_p[2], _B1X, ix,iy,iz-1))+
		  sqr(MRC_F3(u_p[2], _B1Y, ix,iy,iz-1))+
		  sqr(MRC_F3(u_p[2], _B1Z, ix,iy,iz-1)))/ MRC_F3(u_p[2], _RR1, ix,iy,iz-1) );
    
    cAm = sqrtf( (sqr(MRC_F3(u_m[2], _B1X, ix,iy,iz))+
		  sqr(MRC_F3(u_m[2], _B1Y, ix,iy,iz))+
		  sqr(MRC_F3(u_m[2], _B1Z, ix,iy,iz)))/ MRC_F3(u_m[2], _RR1, ix,iy,iz) );


    tmp = sqr(csp) + sqr(cAp);
    cfp = sqrtf( 0.5f * (tmp + sqrtf( sqr( sqr(cAp) + sqr(csp) ) - 
				      (4. * mpermi * sqr(csp * MRC_F3(u_p[2], _B1Z, ix,iy,iz-1)) /  
				       MRC_F3(u_p[2], _RR1, ix,iy,iz-1))) ));      
    //cfp = tmp ;//sqrt(2.0)*sqrt(tmp);
    tmp = sqr(csm) + sqr(cAm);
    cfm = sqrtf( 0.5f * (tmp + sqrtf( sqr( sqr(cAm) + sqr(csm) ) -
				      (4.* mpermi * sqr(csm * MRC_F3(u_m[2], _B1Z, ix,iy,iz)) /  
				       MRC_F3(u_m[2], _RR1, ix,iy,iz))) ));
    //cfm = tmp ;//sqrt(2.0)*sqrt(tmp);

    cp = fmaxf(fmaxf( (MRC_F3(u_p[2], _RV1Z, ix,iy,iz-1) / MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) + cfp,
 		      (MRC_F3(u_m[2], _RV1Z, ix,iy,iz) / MRC_F3(u_m[2], _RR1, ix,iy,iz)) + cfm), 0.f);
    cm = fminf(fminf( (MRC_F3(u_p[2], _RV1Z, ix,iy,iz-1) / MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) - cfp,
		      (MRC_F3(u_m[2], _RV1Z, ix,iy,iz) / MRC_F3(u_m[2], _RR1, ix,iy,iz)) - cfm), 0.f);


    
#if KT == 1
    cp = fmaxf(cp,-cm);
    cm=-cp;
#endif    

    // Flux of _EX,_EY,_EZ through the x faces
    MRC_F3(flux_E[0], 0, ix,iy,iz) = 
      (1.f/(ap - am)) * ( (ap*am) * ( MRC_F3(u_m[0], _B1X, ix,iy,iz) - MRC_F3(u_p[0], _B1X, ix-1,iy,iz)));
    MRC_F3(flux_E[1], 0, ix,iy,iz) = 
      (1.f/(ap - am)) * (   - ap *     MRC_F3(u_p[0], _EZ, ix-1, iy,iz) + am *  MRC_F3(u_m[0], _EZ, ix,iy,iz) + 
			 (ap*am)*   ( MRC_F3(u_m[0], _B1Y, ix,iy,iz) -  MRC_F3(u_p[0], _B1Y, ix-1,iy,iz) ));
    MRC_F3(flux_E[2], 0, ix,iy,iz) = 
      (1.f/(ap - am)) * (     ap  * MRC_F3(u_p[0], _EY, ix-1,iy,iz) - am *  MRC_F3(u_m[0], _EY, ix,iy,iz) + 
			 (ap*am) * ( MRC_F3(u_m[0], _B1Z, ix,iy,iz) -     MRC_F3(u_p[0], _B1Z, ix-1,iy,iz) ));  
    
    // flux of _EX,_EY,_EZ through the y faces    
    // CHECK
    MRC_F3(flux_E[0], 1, ix,iy,iz) =
       (1.f/(bp - bm)) * (     bp  *    MRC_F3(u_p[1], _EZ, ix,iy-1,iz) - bm * MRC_F3(u_m[1], _EZ, ix,iy,iz) + 
		         (bp * bm) *  ( MRC_F3(u_m[1], _B1X, ix,iy,iz) -    MRC_F3(u_p[1], _B1X, ix,iy-1,iz) ));
    MRC_F3(flux_E[1], 1, ix,iy,iz) =
       (1.f/(bp - bm)) * ( (bp * bm)* ( MRC_F3(u_m[1], _B1Y, ix,iy,iz) - MRC_F3(u_p[1], _B1Y, ix,iy-1,iz) ));    
    MRC_F3(flux_E[2], 1, ix,iy,iz) =
       (1.f/(bp - bm)) * (    - bp *  MRC_F3(u_p[1], _EX, ix,iy-1,iz) + bm *  MRC_F3(u_m[1], _EX, ix,iy,iz) + 
		         (bp * bm) *  ( MRC_F3(u_m[1], _B1Z, ix,iy,iz) -     MRC_F3(u_p[1], _B1Z, ix,iy-1,iz) ));  
    
    // flux of _EX,_EY,_EZ through the z faces
    MRC_F3(flux_E[0], 2, ix,iy,iz) = 
      (1.f/(cp - cm))*( - cp  *  MRC_F3(u_p[2], _EY, ix,iy,iz-1)  +  cm * MRC_F3(u_m[2], _EY, ix,iy,iz) + 
			 (cp * cm) * ( MRC_F3(u_m[2], _B1X, ix,iy,iz) - MRC_F3(u_p[2], _B1X, ix,iy,iz-1) ));
    MRC_F3(flux_E[1], 2, ix,iy,iz) =
      (1.f/(cp - cm)) *( cp * MRC_F3(u_p[2], _EX, ix,iy,iz-1) - cm * MRC_F3(u_m[2], _EX, ix,iy,iz) +
		         (cp * cm) * ( MRC_F3(u_m[2], _B1Y, ix,iy,iz) -    MRC_F3(u_p[2], _B1Y, ix,iy,iz-1) ));
    MRC_F3(flux_E[2], 2, ix,iy,iz) = 
      (1.f/(cp - cm))*(  (cp * cm) * ( MRC_F3(u_m[2], _B1Z, ix,iy,iz) - MRC_F3(u_p[2], _B1Z, ix,iy,iz-1) ));    

      for (int m = 0; m < _UU1+1; m++) {
	MRC_F3(flux[m], 0, ix,iy,iz) =
	  (ap * MRC_F3(flux_p[m], 0, ix-1,iy,iz) - am * MRC_F3(flux_m[m], 0, ix,iy,iz)) / (ap - am) +
	  (ap * am) / (ap - am) * (MRC_F3(u_m[0], m, ix ,iy,iz) - MRC_F3(u_p[0], m, ix-1,iy,iz));
	MRC_F3(flux[m], 1, ix,iy,iz) = 
	  (bp  * MRC_F3(flux_p[m], 1, ix,iy-1,iz) - bm * MRC_F3(flux_m[m], 1, ix,iy ,iz)) / (bp - bm) +
	  (bp * bm) / (bp - bm) * (MRC_F3(u_m[1], m, ix,iy ,iz) - MRC_F3(u_p[1], m, ix,iy-1 ,iz));
	MRC_F3(flux[m], 2, ix,iy,iz) =   
	  (cp  * MRC_F3(flux_p[m], 2, ix,iy,iz-1) - cm * MRC_F3(flux_m[m], 2, ix,iy,iz )) / (cp - cm) +
	  (cp * cm) / (cp - cm) * (MRC_F3(u_m[2], m, ix,iy,iz ) - MRC_F3(u_p[2], m, ix,iy,iz-1));
      } 
     } mrc_f3_foreach_end;
 
  for (int f = 0; f < 3; f++) {
    mrc_f3_destroy(u_delta[f]);
    mrc_f3_destroy(u_p[f]);
    mrc_f3_destroy(u_m[f]);
  }
  for (int m = 0; m < _B1Z + 1; m++) {
    mrc_f3_destroy(flux_p[m]);
    mrc_f3_destroy(flux_m[m]);
  }
}

static void
ggcm_mhd_step_cweno_calc_rhs(struct ggcm_mhd_step *step, struct mrc_f3 *rhs,
			     struct mrc_f3 *fld)
{
  struct ggcm_mhd *mhd = step->mhd;
  
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct mrc_f3 *flux[_UU1 + 1];
  for (int m = 0; m < _UU1 + 1; m++) {
    flux[m] = ggcm_mhd_get_fields(mhd, "flux", 3);
  }
  struct mrc_f3 *flux_E[3];
  for (int m = 0; m < 3; m++) {
    flux_E[m] = ggcm_mhd_get_fields(mhd, "flux_E", 3);
  }
   
  calc_cweno_fluxes(flux, flux_E, mhd, fld, crds);
  //zero_flux(flux, mhd, fld); 

  for (int m = 0; m < _UU1+1; m++) {
    calc_neg_divg(rhs, m, flux[m], crds);
  }
  calc_fct_rhs(mhd, rhs, fld, flux_E); 
  
  for (int m = 0; m < _UU1 + 1; m++) {
    mrc_f3_destroy(flux[m]);
  }
  for (int m = 0; m < 3; m++) {
    mrc_f3_destroy(flux_E[m]);
  }

  /*
  struct mrc_io *io=mrc_io_create(MPI_COMM_WORLD); 
  mrc_io_set_param_string(io,"basename","rhs");
  mrc_io_setup(io);
  mrc_io_open(io, "w",0,0);
  for (int m = 0; m < __NR_FLDS; m++){
    char s[10]; 
    sprintf(s,"m%d",m);
    mrc_f3_set_comp_name(rhs,m,s);
  }
  mrc_f3_write(rhs,io);
  mrc_io_close(io);
  Abort2();  
  mrc_io_destroy(io); 
  */
}

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "cweno"

struct ggcm_mhd_step_ops ggcm_mhd_step_cweno_ops = {
  .name        = "cweno",
  .calc_rhs    = ggcm_mhd_step_cweno_calc_rhs,
};

