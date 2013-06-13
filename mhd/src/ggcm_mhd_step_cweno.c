
#include "ggcm_mhd_step_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"

#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_io.h>

#include <assert.h>
#include <math.h>

//#define DEBUG

// ----------------------------------------------------------------------
// ggcm_mhd_get_fields

static struct mrc_fld *
ggcm_mhd_get_fields(struct ggcm_mhd *mhd, const char *name, int nr_comp)
{   
  struct mrc_fld *f3 = mrc_domain_fld_create(mhd->domain, SW_2, NULL);
  mrc_fld_set_name(f3, name);
  mrc_fld_set_nr_comps(f3, nr_comp);
  mrc_fld_setup(f3);
  return f3;
}



// ======================================================================
// ggcm_mhd_step subclass "cweno"

// Define limiter: [1] van Leer(1977) , [2] minmod (Roe 1976) [3] moncen [4] genminmod

#define LMTR 1
#define sign(x) (( x > 0 ) - ( x < 0 ))
// KNP[0] or KT[1]? 
#define KT 0
#define incws 1

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
calc_neg_divg(struct mrc_fld *_rhs, int m, struct mrc_fld *_flux, struct mrc_crds *crds)
{
  struct mrc_fld *rhs = mrc_fld_get_as(_rhs, "float");
  struct mrc_fld *flux = mrc_fld_get_as(_flux, "float");

  mrc_fld_foreach(rhs, ix, iy, iz, 0, 0) {
    int ind[3] = { ix, iy, iz };

    MRC_F3(rhs, m, ix, iy, iz) = 0.;
    for(int i=0; i<3; i++) {
      int dind[3] = {0, 0, 0};
      dind[i] = 1;      

      float cw =(MRC_CRD(crds, i, ind[i]+1) - MRC_CRD(crds, i, ind[i]));
      MRC_F3(rhs, m, ix, iy, iz) -=( MRC_F3(flux, i, ix+dind[0],iy+dind[1],iz+dind[2]) -
      				     MRC_F3(flux, i, ix,iy,iz))/ cw;
      //      assert(isfinite(MRC_F3(rhs, m, ix,iy,iz)));
    }
  } mrc_fld_foreach_end; 

  mrc_fld_put_as(rhs, _rhs);
  mrc_fld_put_as(flux, _flux);
}

// ----------------------------------------------------------------------
// calc_fluxes_per_face
//
// this calculates fluxes on the face i using reconstructed variables (fld) that are
// given on the respective face

static void
calc_fluxes_per_face(struct mrc_fld **_flux, struct ggcm_mhd *mhd, struct mrc_fld *_fld, int i)
{
  float mpermi = 1.f;
  float gamma = mhd->par.gamm;
  float d_i = mhd->par.d_i;

  struct mrc_fld *fld = mrc_fld_get_as(_fld, "float");
  struct mrc_fld *flux[5];
  for (int m = 0; m < 5; m++) {
    flux[m] = mrc_fld_get_as(_flux[m], "float");
  }
  mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {
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
    
    // mass consv. 
    MRC_F3(flux[_RR1], i, ix,iy,iz) = MRC_F3(fld, _RV1X+i, ix,iy,iz);
    
    // momentum eq. 
    for (int j = 0; j < 3; j++) {
      
      MRC_F3(flux[_RV1X+i], j, ix,iy,iz) = 
	rhoi * MRC_F3(fld, _RV1X+j, ix,iy,iz) * MRC_F3(fld, _RV1X+i, ix,iy,iz) +
	((j == i) ? pp : 0.) + 
	((j == i) ? BB : 0.) - mpermi * (MRC_F3(fld, _B1X+i, ix,iy,iz) * MRC_F3(fld, _B1X+j, ix,iy,iz));
    }
    
    // energy eq. 
    MRC_F3(flux[_UU1], i, ix,iy,iz) =
      ( ((MRC_F3(fld, _UU1, ix,iy,iz) + pp + BB)*MRC_F3(fld, _RV1X+i, ix,iy,iz))-
	(mpermi * mB * MRC_F3(fld, _B1X+i, ix,iy,iz)) + 
	(d_i * ( -0.5*MRC_F3(fld, _JX+i, ix,iy,iz)*BB - MRC_F3(fld, _B1X+i, ix,iy,iz)*JB)) ) * rhoi;

  } mrc_fld_foreach_end;

  mrc_fld_put_as(fld, _fld);
  for (int m = 0; m < 5; m++) {
    mrc_fld_put_as(flux[m], _flux[m]);
  }
}

// ----------------------------------------------------------------------
// calc_cweno_fluxes
//
// calculates CWENO fluxes on faces in flux, flux_E, from the original state
// vector u (which is cell centered / on the Yee grid)

static void
calc_cweno_fluxes(struct mrc_fld **_flux, struct mrc_fld *_flux_E[3], struct ggcm_mhd *mhd,
		  struct mrc_fld *_u, struct mrc_crds *crds)
{
  float mpermi = 1.f;
  float gamma = mhd->par.gamm;
  float d_i = mhd->par.d_i;
  float eta = mhd->par.diffco;

  ggcm_mhd_fill_ghosts(mhd, _u, 0, mhd->time);

  struct mrc_fld *u = mrc_fld_get_as(_u, "float");

  // initialize deltas for reconstruction
  struct mrc_fld *_u_delta[3], *u_delta[3];
  for (int f = 0; f < 3; f++) {
    _u_delta[f] = ggcm_mhd_get_fields(mhd, "u_delta", _B1Z + 1);
    u_delta[f] = mrc_fld_get_as(_u_delta[f], "float");
  }

#if LMTR == 1
  // TVD slope limiter  (van Leer 1977) geometric mean
  // (dxu_)ijk = max{ (u_i+1jk-u_ijk)*(u_ijk-u_i-1jk)  / (u_i+1jk - u_i-1jk)}
  mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
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
    } mrc_fld_foreach_end;
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
      mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
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
      } mrc_fld_foreach_end;      
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
      mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
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
      } mrc_fld_foreach_end;      
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
      mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
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
      } mrc_fld_foreach_end;      
    }
    dind[i] = 0;
  }

#endif

#ifdef DEBUG
  mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
    assert(MRC_F3(u, 0, ix,iy,iz) > 0.);
  } mrc_fld_foreach_end;
#endif

  for (int f = 0; f < 3; f++) {
    mrc_fld_put_as(u_delta[f], _u_delta[f]);
  }

  ggcm_mhd_fill_ghosts(mhd, _u_delta[0], 0, mhd->time);
  ggcm_mhd_fill_ghosts(mhd, _u_delta[1], 0, mhd->time);
  ggcm_mhd_fill_ghosts(mhd, _u_delta[2], 0, mhd->time);

#ifdef DEBUG
  mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
    assert(MRC_F3(u, 0, ix,iy,iz) > 0.);
  } mrc_fld_foreach_end;
#endif

  //initialize cell surface center variables			    
  struct mrc_fld *_u_p[3], *_u_m[3], *u_p[3], *u_m[3];
  for (int f = 0; f < 3; f++) {
    _u_p[f] = ggcm_mhd_get_fields(mhd, "u_p", _JZ + 1);
    _u_m[f] = ggcm_mhd_get_fields(mhd, "u_m", _JZ + 1);
  }

  for (int f = 0; f < 3; f++) {
    u_delta[f] = mrc_fld_get_as(_u_delta[f], "float");
    u_p[f] = mrc_fld_get_as(_u_p[f], "float");
    u_m[f] = mrc_fld_get_as(_u_m[f], "float");
  }

  // Reonstruction    UijkE  = u_ijk + (dxu_)ijk    UijkW = u_ijk - (dxu_)ijk
  mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
    for (int m = 0; m < _UU1+1; m++) {
      // defined like this, both u_p and u_m are coplaner when indices are the same

      MRC_F3(u_p[0], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) + MRC_F3(u_delta[0], m, ix,iy,iz);
      MRC_F3(u_p[1], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) + MRC_F3(u_delta[1], m, ix,iy,iz);
      MRC_F3(u_p[2], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) + MRC_F3(u_delta[2], m, ix,iy,iz);
      MRC_F3(u_m[0], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) - MRC_F3(u_delta[0], m, ix,iy,iz);
      MRC_F3(u_m[1], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) - MRC_F3(u_delta[1], m, ix,iy,iz);
      MRC_F3(u_m[2], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) - MRC_F3(u_delta[2], m, ix,iy,iz);
    }
  } mrc_fld_foreach_end;
  
  // calculation of cell surface center averages for B 
  // note that B is staggered here.. so that index 1 for B is 1-1/2 i.e. cell surface for V
  mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
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
  } mrc_fld_foreach_end;
  
  float *bdx3 = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD3);
  float *bdy3 = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD3);
  float *bdz3 = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD3);

  for (int i = 0; i < 3; i++) {    
    mrc_fld_foreach(u, ix,iy,iz, 2, 2) {	
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
      
    } mrc_fld_foreach_end;
  }
  
  // Calculation of cell surface center  averages for E using v and B just calculated 
  //  E^N(SWETB) = -V^N(SWETB) X B^N(SWETB)      + eta*J  [  + di J x B ] 
  mrc_fld_foreach(u, ix, iy, iz, 1, 1) {  
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
  } mrc_fld_foreach_end;
  
  mrc_fld_put_as(u, _u);
  for (int f = 0; f < 3; f++) {
    mrc_fld_put_as(u_delta[f], _u_delta[f]);
    mrc_fld_put_as(u_p[f], _u_p[f]);
    mrc_fld_put_as(u_m[f], _u_m[f]);
  }
  //mrc_fld_destroy(j);

  // initialize flux functions
  struct mrc_fld *_flux_p[_B1Z + 1], *_flux_m[_B1Z + 1];
  for (int m = 0; m < _B1Z + 1; m++) {
    _flux_p[m] = ggcm_mhd_get_fields(mhd, "flux_p", 3);
    _flux_m[m] = ggcm_mhd_get_fields(mhd, "flux_m", 3);
  }
  
  // calculate fluxes per face (small f's) using reconstructed 
  // variables U^N(SWETB) and B^N(SWETB) = (Bx,By,Bz)^N(SEWTB)
  for (int f = 0; f < 3; f++) {
    calc_fluxes_per_face(_flux_p, mhd, u_p[f], f);
    calc_fluxes_per_face(_flux_m, mhd, u_m[f], f);
  }
  
  struct mrc_fld *flux_p[_B1Z + 1], *flux_m[_B1Z + 1];
  for (int m = 0; m < _B1Z + 1; m++) {
    flux_p[m] = mrc_fld_get_as(_flux_p[m], "float");
    flux_m[m] = mrc_fld_get_as(_flux_m[m], "float");
  }

  //  float cs;
  float csp; 
  float csm; 
  float tmp;
  float cAp;
  float cAm;
  float cfp; 
  float cfm; 
  float cwp; // whistler speed estimate
  float cwm; 
  float ap; 
  float am; 
  float bp; 
  float bm; 
  float cp;
  float cm;
  float ppp;
  float ppm; 
  float rhoip;
  float rhoim; 

  float *bdx1 = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD1);
  float *bdy1 = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD1);
  float *bdz1 = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD1);

  u = mrc_fld_get_as(_u, "float");
  for (int f = 0; f < 3; f++) {
    u_p[f] = mrc_fld_get_as(_u_p[f], "float");
    u_m[f] = mrc_fld_get_as(_u_m[f], "float");
  }

  struct mrc_fld *flux_E[3];
  for (int m = 0; m < 3; m++) {
    flux_E[m] = mrc_fld_get_as(_flux_E[m], "float");
  }
  struct mrc_fld *flux[5];
  for (int m = 0; m < 5; m++) {
    flux[m] = mrc_fld_get_as(_flux[m], "float");
  }
  
  mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
    
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
 						      
    tmp = sqr(csp) + sqr(cAp);
    cfp = sqrtf( 0.5 * ( tmp + sqrtf( sqr( sqr(cAp) + sqr(csp) )
				     - (4. * mpermi * sqr(csp * MRC_F3(u, _B1X, ix-1,iy,iz)) /  
					MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) )) );      
    tmp = sqr(csm) + sqr(cAm);
    cfm = sqrtf( 0.5 * ( tmp + sqrtf( sqr( sqr(cAm) + sqr(csm) )
				     - (4.* mpermi * sqr(csm * MRC_F3(u, _B1X, ix,iy,iz)) /  
					MRC_F3(u_m[0], _RR1, ix,iy,iz))  )) );

    /*
    ap = fmaxf(fmaxf((MRC_F3(u_p[0], _RV1X, ix-1,iy,iz) / MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) + cfp,
		     (MRC_F3(u_m[0], _RV1X, ix,iy,iz) / MRC_F3(u_m[0], _RR1, ix,iy,iz)) + cfm), 0.f);

    am = fminf(fminf( (MRC_F3(u_p[0], _RV1X, ix-1,iy,iz) / MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) - cfp,
		     (MRC_F3(u_m[0], _RV1X, ix,iy,iz) / MRC_F3(u_m[0], _RR1,ix,iy,iz)) - cfm), 0.f);
    */

    // cw =  ( |B|*pi )  /  ( e * rho * delta x ) 
    cwp = d_i * cAp * sqrtf(1./MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) * bdx1[ix+2]  ;
    cwm = d_i * cAm * sqrtf(1./MRC_F3(u_m[0], _RR1, ix,iy,iz)) * bdx1[ix+2]  ; 

#if incws == 0
    cwp = 0;
    cwm = 0;
#endif

    //printf("cAp %f cwp %f %f \n", cAp, d_i*(cAp*M_PI*gmhd->crdx[BD1][ix+2]) 
    // sqrtf(MRC_F3(u_p[0], _RR1, ix-1,iy,iz)),gmhd->crdx[BD1][ix+2]  ); 

    ap = fmaxf(fmaxf((MRC_F3(u_p[0], _RV1X, ix-1,iy,iz) / MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) + cfp + cwp,
		     (MRC_F3(u_m[0], _RV1X, ix,iy,iz) / MRC_F3(u_m[0], _RR1, ix,iy,iz)) + cfm + cwm ), 0.f);

    am = fminf(fminf( (MRC_F3(u_p[0], _RV1X, ix-1,iy,iz) / MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) - cfp - cwp,
		     (MRC_F3(u_m[0], _RV1X, ix,iy,iz) / MRC_F3(u_m[0], _RR1,ix,iy,iz)) - cfm - cwm), 0.f);
    
   
    
    
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
    tmp = sqr(csm) + sqr(cAm);
    cfm = sqrtf( 0.5 * (tmp + sqrtf( sqr( sqr(cAm) + sqr(csm) ) -  
				     (4.f * mpermi * sqr(csm * MRC_F3(u, _B1Y, ix,iy,iz)) /  
				      MRC_F3(u_m[1], _RR1, ix,iy,iz))) ));

    /*
    bp = fmaxf(fmaxf( (MRC_F3(u_p[1], _RV1Y, ix,iy-1,iz) / MRC_F3(u_p[1], _RR1, ix,iy-1,iz)) + cfp,
		      (MRC_F3(u_m[1], _RV1Y, ix,iy,iz) / MRC_F3(u_m[1], _RR1, ix,iy,iz)) + cfm), 0.f);
    bm = fminf(fminf( (MRC_F3(u_p[1], _RV1Y, ix,iy-1,iz) / MRC_F3(u_p[1], _RR1, ix,iy-1,iz)) - cfp,
		      (MRC_F3(u_m[1], _RV1Y, ix,iy,iz) / MRC_F3(u_m[1], _RR1, ix,iy,iz)) - cfm), 0.f);
    */
    cwp =  d_i * cAp * sqrtf(1./MRC_F3(u_p[1], _RR1, ix,iy-1,iz)) * bdy1[iy+2];
    cwm =  d_i * cAm * sqrtf(1./MRC_F3(u_m[1], _RR1, ix,iy,iz)) * bdy1[iy+2]; 
    
   
#if incws == 0
    cwp = 0;
    cwm = 0;
#endif
    
    bp = fmaxf(fmaxf( (MRC_F3(u_p[1], _RV1Y, ix,iy-1,iz) / MRC_F3(u_p[1], _RR1, ix,iy-1,iz)) + cfp + cwp,
		      (MRC_F3(u_m[1], _RV1Y, ix,iy,iz) / MRC_F3(u_m[1], _RR1, ix,iy,iz)) + cfm + cwm), 0.f);
    bm = fminf(fminf( (MRC_F3(u_p[1], _RV1Y, ix,iy-1,iz) / MRC_F3(u_p[1], _RR1, ix,iy-1,iz)) - cfp - cwp,
		      (MRC_F3(u_m[1], _RV1Y, ix,iy,iz) / MRC_F3(u_m[1], _RR1, ix,iy,iz)) - cfm - cwm), 0.f);

  

   
    
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

    tmp = sqr(csm) + sqr(cAm);
    cfm = sqrtf( 0.5f * (tmp + sqrtf( sqr( sqr(cAm) + sqr(csm) ) -
				      (4.* mpermi * sqr(csm * MRC_F3(u_m[2], _B1Z, ix,iy,iz)) /  
				       MRC_F3(u_m[2], _RR1, ix,iy,iz))) ));

    /*
    cp = fmaxf(fmaxf( (MRC_F3(u_p[2], _RV1Z, ix,iy,iz-1) / MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) + cfp,
 		      (MRC_F3(u_m[2], _RV1Z, ix,iy,iz) / MRC_F3(u_m[2], _RR1, ix,iy,iz)) + cfm), 0.f);
    cm = fminf(fminf( (MRC_F3(u_p[2], _RV1Z, ix,iy,iz-1) / MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) - cfp,
		      (MRC_F3(u_m[2], _RV1Z, ix,iy,iz) / MRC_F3(u_m[2], _RR1, ix,iy,iz)) - cfm), 0.f);
    */

    
    //cwp = d_i * cAp * M_PI * sqrtf(MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) * bdz1[iz+2]  ;
    //cwm = d_i * cAm * M_PI * sqrtf(MRC_F3(u_m[2], _RR1, ix,iy,iz)) * bdz1[iz+2]  ; 

  
    cwp = d_i * cAp * sqrtf(1./MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) * bdz1[iz+2] ;
    cwm = d_i * cAm * sqrtf(1./MRC_F3(u_m[2], _RR1, ix,iy,iz)) * bdz1[iz+2] ; 


    
    //printf("cAp %f cwp %f %f \n", cAp,cAp * M_PI * sqrtf(MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) * gmhd->crdz[BD1][iz+2] , ); 

#if incws == 0
    cwp = 0;
    cwm = 0;
#endif

    cp = fmaxf(fmaxf( (MRC_F3(u_p[2], _RV1Z, ix,iy,iz-1) / MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) + cfp + cwp,
 		      (MRC_F3(u_m[2], _RV1Z, ix,iy,iz) / MRC_F3(u_m[2], _RR1, ix,iy,iz)) + cfm + cwm), 0.f);
    cm = fminf(fminf( (MRC_F3(u_p[2], _RV1Z, ix,iy,iz-1) / MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) - cfp - cwp,
		      (MRC_F3(u_m[2], _RV1Z, ix,iy,iz) / MRC_F3(u_m[2], _RR1, ix,iy,iz)) - cfm - cwm), 0.f);

    
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
			 (ap*am) * ( MRC_F3(u_m[0], _B1Z, ix,iy,iz) - MRC_F3(u_p[0], _B1Z, ix-1,iy,iz) ));  
    
    // flux of _EX,_EY,_EZ through the y faces    
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
     } mrc_fld_foreach_end;
 
  for (int m = 0; m < 3; m++) {
    mrc_fld_put_as(flux_E[m], _flux_E[m]);
  }
  for (int m = 0; m < 5; m++) {
    mrc_fld_put_as(flux[m], _flux[m]);
  }
  mrc_fld_put_as(u, _u);
  for (int m = 0; m < 3; m++) {
    mrc_fld_put_as(u_p[m], _u_p[m]);
    mrc_fld_put_as(u_m[m], _u_m[m]);
  }
  for (int m = 0; m < _B1Z + 1; m++) {
    mrc_fld_put_as(flux_p[m], _flux_p[m]);
    mrc_fld_put_as(flux_m[m], _flux_m[m]);
  }
  for (int f = 0; f < 3; f++) {
    mrc_fld_destroy(_u_delta[f]);
    mrc_fld_destroy(_u_p[f]);
    mrc_fld_destroy(_u_m[f]);
  }
  for (int m = 0; m < _B1Z + 1; m++) {
    mrc_fld_destroy(flux_p[m]);
    mrc_fld_destroy(flux_m[m]);
  }
}

#if 0

// ----------------------------------------------------------------------
// fill_ghost_fld 
// 
// This fills ghost cells for fld objects that have been duplicated in the 
// time-stepping routine (c.f. mrc_ts_rk2.c). Without this, zero-values at 
// boundaries will give inf and nan values for non-periodic boudnaries.

static void
fill_ghost_fld(struct ggcm_mhd *mhd, struct mrc_fld *_fld)
{
  struct mrc_fld *f3 = mrc_fld_get_as(mhd->fld, "float");
  struct mrc_fld *fld = mrc_fld_get_as(_fld, "float");
  int sw = 2; 
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, 0, &info);

  int bc[3];
  mrc_domain_get_param_int(mhd->domain, "bcx", &bc[0]); // FIXME in libmrc
  mrc_domain_get_param_int(mhd->domain, "bcy", &bc[1]);
  mrc_domain_get_param_int(mhd->domain, "bcz", &bc[2]);
  const int *dims = mrc_fld_dims(f3);
  int nx = dims[0], ny = dims[1], nz = dims[2];

  if ( bc[0] != BC_PERIODIC ) {    
    if (info.off[0] == 0) { 
      for (int m =0; m < __NR_FLDS+1; m++) { 
	for (int ix = -sw; ix < 0; ix++) {  
	  for (int iy = -sw; iy < ny + sw; iy++) {
	    for (int iz = -sw; iz < nz + sw; iz++) {
	      MRC_F3(fld,m, ix,iy,iz) = MRC_F3(f3,m, ix,iy,iz); 
	    }	     
	  }
	}
      }
    }    
    if (info.off[0] + info.ldims[0] == gdims[0]) { 
      for ( int m = 0; m < __NR_FLDS+1; m++ ) {
	for (int ix = nx; ix < nx + sw; ix++) {  
	  for (int iy = -sw; iy < ny + sw; iy++) {
	    for (int iz = -sw; iz < nz + sw; iz++) {
	      MRC_F3(fld,m, ix,iy,iz) = MRC_F3(f3,m, ix,iy,iz);
		}		
	  }
	}
      }
    }    
  }  

  if ( bc[1] != BC_PERIODIC ) {    
    if (info.off[1] == 0) { 
      for (int m =0; m < __NR_FLDS+1; m++) { 
	for (int iy = -sw; iy < 0; iy++) {  
	  for (int ix = -sw; ix < nx + sw; ix++) {
	    for (int iz = -sw; iz < nz + sw; iz++) {
	      MRC_F3(fld,m, ix,iy,iz) = MRC_F3(f3,m, ix,iy,iz); 
	    }
	  }
	}
      }
    }    
    if (info.off[1] + info.ldims[1] == gdims[1]) { 
      for ( int m = 0; m < __NR_FLDS+1; m++ ) {
	for (int iy = ny; iy < ny + sw; iy++) {  
	  for (int ix = -sw; ix < nx + sw; ix++) {
	    for (int iz = -sw; iz < nz + sw; iz++) {
	      MRC_F3(fld,m, ix,iy,iz) = MRC_F3(f3,m, ix,iy,iz);
		}		
	  }
	}
      }
    }    
  }  

  if ( bc[2] != BC_PERIODIC ) {    
    if (info.off[2] == 0) { 
      for (int m =0; m < __NR_FLDS+1; m++) { 
	for (int iz = -sw; iz < 0; iz++) {  
	  for (int iy = -sw; iy < ny + sw; iy++) {
	    for (int ix = -sw; ix < nx + sw; ix++) {
	      MRC_F3(fld,m, ix,iy,iz) = MRC_F3(f3,m, ix,iy,iz); 
	    }
	  }
	}
      }
    }    
    if (info.off[2] + info.ldims[2] == gdims[2]) { 
      for ( int m = 0; m < __NR_FLDS+1; m++ ) {
	for (int iz = nz; iz < nz + sw; iz++) {  
	  for (int iy = -sw; iy < ny + sw; iy++) {
	    for (int ix = -sw; ix < nx + sw; ix++) {
	      MRC_F3(fld,m, ix,iy,iz) = MRC_F3(f3,m, ix,iy,iz);
		}		
	  }
	}
       }
    }    
  }  

  mrc_fld_put_as(f3, mhd->fld);
  mrc_fld_put_as(fld, _fld);
}

#endif

static void
calc_fct_rhs(struct ggcm_mhd *mhd, struct mrc_fld *_rhs, struct mrc_fld *_fld, struct mrc_fld *_flux_E[3])
{  
  // compute edge centered electric fields by interpolation (iv)  here tmp_fld are edge centered 
  // so that e.g.  MRC_F3(tmp_fld, 0, 0, 0, 0) is E_x 0,-1/2,-1/2  
  //         i.e.  MRC_F3(tmp_fld, 0, ix,iy,iz) is E_x ix,iy-1/2,iz-1/2   
  //           and MRC_F3(tmp_fld, 1, 0, 0, 0) is E_y -1/2,0,-1/2  etc etc 
                       

  struct mrc_ddc *ddc = mrc_domain_get_ddc(mhd->domain);
  mrc_ddc_set_param_int(ddc, "max_n_fields", 3);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  /* struct mrc_fld *E_ec = ggcm_mhd_get_fields(mhd, "E_ec", 3); */

  struct mrc_fld *fld = mrc_fld_get_as(_fld, "float");
  struct mrc_fld *fex = mrc_fld_get_as(_flux_E[0], "float");
  struct mrc_fld *fey = mrc_fld_get_as(_flux_E[1], "float");
  struct mrc_fld *fez = mrc_fld_get_as(_flux_E[2], "float");

  //initialize cell edge center Electric field structure      
  mrc_fld_foreach(fld, ix, iy,  iz, 1, 1) { 

    /*
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
    */

    MRC_F3(fld, _EX, ix, iy, iz) = 0.25*(-MRC_F3(fez, 1, ix, iy, iz) -
					 MRC_F3(fez, 1, ix, iy, iz-1) +
					 MRC_F3(fey, 2, ix, iy, iz) +
					 MRC_F3(fey, 2, ix, iy-1, iz));    
    //MRC_F3(fld,_EX, ix, iy, iz) = MRC_F3(E_ec, 0, ix, iy, iz);    
    MRC_F3(fld, _EY, ix, iy, iz) = 0.25*( MRC_F3(fez, 0, ix, iy, iz) +
					 MRC_F3(fez, 0, ix, iy, iz-1) -
					 MRC_F3(fex, 2, ix, iy, iz) -
					 MRC_F3(fex, 2, ix-1, iy, iz));    
    //MRC_F3(fld,_EY, ix, iy, iz) = MRC_F3(E_ec, 1, ix, iy, iz);        
    MRC_F3(fld, _EZ, ix, iy, iz) = 0.25*(-MRC_F3(fey, 0, ix, iy, iz) -
					 MRC_F3(fey, 0, ix, iy-1, iz) +
					 MRC_F3(fex, 1, ix, iy, iz) +
					 MRC_F3(fex, 1, ix-1, iy, iz));    
    //MRC_F3(fld,_EZ, ix, iy, iz) = MRC_F3(E_ec, 2, ix, iy, iz);

  } mrc_fld_foreach_end;
  
  /////////////////////////////////////////////////
  //  Use field here rather than E_ec hold the electric field so that 
  //  we can (correctly apply the boundary condition). 
  /////////////////////////////////////////////////

  mrc_fld_put_as(fld, _fld);
  mrc_fld_put_as(fex, _flux_E[0]);
  mrc_fld_put_as(fey, _flux_E[1]);
  mrc_fld_put_as(fez, _flux_E[2]);

  //  fill_ghost_fld(mhd, fld);

  struct mrc_fld *rhs = mrc_fld_get_as(_rhs, "float");
  fld = mrc_fld_get_as(_fld, "float");

  mrc_fld_foreach(rhs, ix, iy,  iz, 1, 1) {
    
    B1X(rhs, ix, iy, iz) =  
      (-(( MRC_F3(fld, _EZ, ix, iy+1, iz) - MRC_F3(fld, _EZ, ix, iy, iz) ) /
	 (.5f*( MRC_CRDY(crds, iy+1) - MRC_CRDY(crds, iy-1) ))) 
       +(( MRC_F3(fld, _EY, ix, iy, iz+1) - MRC_F3(fld, _EY, ix, iy, iz) ) /
	 (.5f*( MRC_CRDZ(crds, iz+1) - MRC_CRDZ(crds, iz-1) )))); 
    
    B1Y(rhs, ix, iy, iz) = 
      (-(( MRC_F3(fld, _EX, ix, iy, iz+1) - MRC_F3(fld, _EX, ix, iy, iz) ) /
	 (.5f*( MRC_CRD(crds, 2, iz+1) - MRC_CRD(crds, 2, iz-1) )))
       +(( MRC_F3(fld, _EZ, ix+1, iy, iz) - MRC_F3(fld, _EZ, ix, iy, iz) ) /
	 (.5f*(MRC_CRD(crds, 0, ix+1) - MRC_CRD(crds, 0, ix-1) )))); 
    
    B1Z(rhs, ix, iy, iz) = 
      (-(( MRC_F3( fld, _EY, ix+1, iy, iz) - MRC_F3(fld, _EY, ix, iy, iz) ) /
	 (.5f*( MRC_CRD(crds, 0, ix+1) - MRC_CRD(crds, 0, ix-1) )))  
       +(( MRC_F3( fld, _EX, ix, iy+1, iz) - MRC_F3(fld, _EX, ix, iy, iz) ) /
	 (.5f*( MRC_CRD(crds, 1, iy+1) - MRC_CRD(crds, 1, iy-1) )))); 

    /*
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
    */
  } mrc_fld_foreach_end;
    
  mrc_fld_put_as(fld, _fld);
  mrc_fld_put_as(rhs, _rhs);

  /* mrc_fld_destroy(E_ec);  */
}

static void
ggcm_mhd_step_cweno_calc_rhs(struct ggcm_mhd_step *step, struct mrc_fld *rhs,
			     struct mrc_fld *fld)
{
  struct ggcm_mhd *mhd = step->mhd;  
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  
  //fill_ghost_fld(mhd, fld);

#ifdef DEBUG
  {
    static struct ggcm_mhd_diag *diag;
    static int cnt;
    if (!diag) {
      diag = ggcm_mhd_diag_create(ggcm_mhd_comm(mhd));
      ggcm_mhd_diag_set_param_obj(diag, "mhd", mhd);
      ggcm_mhd_diag_set_param_string(diag, "run", "fld");
      ggcm_mhd_diag_setup(diag);
      ggcm_mhd_diag_view(diag);
    }
    ggcm_mhd_diag_run_now(diag, mhd->fld, DIAG_TYPE_3D, cnt++);
  }
#endif

  struct mrc_fld *flux[_UU1 + 1];
  for (int m = 0; m < _UU1 + 1; m++) {
    flux[m] = ggcm_mhd_get_fields(mhd, "flux", 3);
  }
  struct mrc_fld *flux_E[3];
  for (int m = 0; m < 3; m++) {
    flux_E[m] = ggcm_mhd_get_fields(mhd, "flux_E", 3);
  }
   
  calc_cweno_fluxes(flux, flux_E, mhd, fld, crds);

  for (int m = 0; m < _UU1+1; m++) {
    calc_neg_divg(rhs, m, flux[m], crds);
  }
  //fill_ghost_fld(mhd, fld);
  calc_fct_rhs(mhd, rhs, fld, flux_E); 

#ifdef DEBUG
  {
    static struct ggcm_mhd_diag *diag;
    static int cnt;
    if (!diag) {
      diag = ggcm_mhd_diag_create(ggcm_mhd_comm(mhd));
      ggcm_mhd_diag_set_param_obj(diag, "mhd", mhd);
      ggcm_mhd_diag_set_param_string(diag, "run", "rhs");
      ggcm_mhd_diag_setup(diag);
      ggcm_mhd_diag_view(diag);
    }
    ggcm_mhd_diag_run_now(diag, mhd->fld, DIAG_TYPE_3D, cnt++);
  }
#endif

  for (int m = 0; m < _UU1 + 1; m++) {
    mrc_fld_destroy(flux[m]);
  }
  for (int m = 0; m < 3; m++) {
    mrc_fld_destroy(flux_E[m]);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "cweno"

struct ggcm_mhd_step_ops ggcm_mhd_step_cweno_ops = {
  .name        = "cweno",
  .calc_rhs    = ggcm_mhd_step_cweno_calc_rhs,
};

