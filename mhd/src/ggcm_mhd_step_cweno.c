
#include "ggcm_mhd_step_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"

#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_io.h>

#include <assert.h>
#include <stdio.h>
#include <math.h>

#define DEBUG

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
#define incws 0

enum {
  // reuse B in the _fluxes_ (only) to store E field
  _EX = _B1X,
  _EY,
  _EZ,

  _JX = 11,
  _JY,
  _JZ,

  __NR_FLDS,
};

#define FLUX(flux, f, m, ix,iy,iz) MRC_F3((flux)[f], m, ix,iy,iz)

static void
calc_neg_divg(struct ggcm_mhd *mhd, struct mrc_fld *_rhs, struct mrc_fld *_flux[3])
{
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  struct mrc_fld *rhs = mrc_fld_get_as(_rhs, "mhd_fc_float");

  struct mrc_fld *flux[3];
  for (int f = 0; f < 3; f++) {
    flux[f] = mrc_fld_get_as(_flux[f], "float");
  }

  for (int m = 0; m <= _UU1; m++) {
    mrc_fld_foreach(rhs, ix, iy, iz, 0, 0) {
      int ind[3] = { ix, iy, iz };
      
      MRC_F3(rhs, m, ix, iy, iz) = 0.;
      for(int i=0; i<3; i++) {
	int dind[3] = {0, 0, 0};
	dind[i] = 1;      
	
	MRC_F3(rhs, m, ix, iy, iz) -=
	  (FLUX(flux, i, m, ix+dind[0],iy+dind[1],iz+dind[2]) - 
	   FLUX(flux, i, m, ix,iy,iz))
	  / (MRC_CRD(crds, i, ind[i]+1) - MRC_CRD(crds, i, ind[i]));
	assert(isfinite(MRC_F3(rhs, m, ix, iy, iz)));
      }
    } mrc_fld_foreach_end; 
  }

  for (int f = 0; f < 3; f++) {
    mrc_fld_put_as(flux[f], _flux[f]);
  }
  mrc_fld_put_as(rhs, _rhs);
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
  struct mrc_fld *flux[3];
  for (int f = 0; f < 3; f++) {
    flux[f] = mrc_fld_get_as(_flux[f], "float");
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
    FLUX(flux, i, _RR1, ix,iy,iz) = MRC_F3(fld, _RV1X+i, ix,iy,iz);
    
    // momentum eq. 
    for (int j = 0; j < 3; j++) {
      FLUX(flux, j, _RV1X+i, ix,iy,iz) = 
	rhoi * MRC_F3(fld, _RV1X+j, ix,iy,iz) * MRC_F3(fld, _RV1X+i, ix,iy,iz) +
	((j == i) ? pp : 0.) + 
	((j == i) ? BB : 0.) - mpermi * (MRC_F3(fld, _B1X+i, ix,iy,iz) * MRC_F3(fld, _B1X+j, ix,iy,iz));
    }
    
    // energy eq. 
    FLUX(flux, i, _UU1, ix,iy,iz) =
      ( ((MRC_F3(fld, _UU1, ix,iy,iz) + pp + BB)*MRC_F3(fld, _RV1X+i, ix,iy,iz))-
	(mpermi * mB * MRC_F3(fld, _B1X+i, ix,iy,iz)) + 
	(d_i * ( -0.5*MRC_F3(fld, _JX+i, ix,iy,iz)*BB - MRC_F3(fld, _B1X+i, ix,iy,iz)*JB)) ) * rhoi;

  } mrc_fld_foreach_end;

  mrc_fld_put_as(fld, _fld);
  for (int f = 0; f < 3; f++) {
    mrc_fld_put_as(flux[f], _flux[f]);
  }
}

// ----------------------------------------------------------------------
// limit_1
//
// TVD slope limiter  (van Leer 1977) harmonic mean
// (dxu_)ijk = max{ (u_i+1jk-u_ijk)*(u_ijk-u_i-1jk)  / (u_i+1jk - u_i-1jk)}

static void __unused
limit_1(struct mrc_fld *u_delta[3], struct mrc_fld *u)
{
  mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
    for (int m = 0; m <= _UU1; m++) {
      MRC_F3(u_delta[0], m, ix,iy,iz) = 
	fmaxf((MRC_F3(u, m, ix+1,iy,iz) - MRC_F3(u, m, ix  ,iy,iz)) *
	      (MRC_F3(u, m, ix  ,iy,iz) - MRC_F3(u, m, ix-1,iy,iz)) , 0.f) /
	(MRC_F3(u, m, ix+1,iy,iz) - MRC_F3(u, m, ix-1,iy,iz));
      MRC_F3(u_delta[1], m, ix,iy,iz) = 
	fmaxf((MRC_F3(u, m, ix,iy+1,iz) - MRC_F3(u, m, ix,iy  ,iz)) *
	      (MRC_F3(u, m, ix,iy  ,iz) - MRC_F3(u, m, ix,iy-1,iz)), 0.f) / 
	(MRC_F3(u, m, ix,iy+1,iz) - MRC_F3(u, m, ix,iy-1,iz));
      MRC_F3(u_delta[2], m, ix,iy,iz) = 
	fmaxf((MRC_F3(u, m, ix,iy,iz+1) - MRC_F3(u, m, ix,iy,iz  )) *
	      (MRC_F3(u, m, ix,iy,iz  ) - MRC_F3(u, m, ix,iy,iz-1)), 0.f) /
	(MRC_F3(u, m, ix,iy,iz+1) - MRC_F3(u, m, ix,iy,iz-1));
      // FIXME, need to make sure NaN -> 0
      if (!isfinite(MRC_F3(u_delta[0], m, ix,iy,iz))) MRC_F3(u_delta[0], m, ix,iy,iz) = 0.f;
      if (!isfinite(MRC_F3(u_delta[1], m, ix,iy,iz))) MRC_F3(u_delta[1], m, ix,iy,iz) = 0.f;
      if (!isfinite(MRC_F3(u_delta[2], m, ix,iy,iz))) MRC_F3(u_delta[2], m, ix,iy,iz) = 0.f;
    }
    for (int m = 0; m < 3; m++) {
      MRC_F3(u_delta[0], _B1X + m, ix,iy,iz) = 
	fmaxf((B1XYZ(u, m, ix+1,iy,iz) - B1XYZ(u, m, ix  ,iy,iz)) *
	      (B1XYZ(u, m, ix  ,iy,iz) - B1XYZ(u, m, ix-1,iy,iz)), 0.f) /
	(B1XYZ(u, m, ix+1,iy,iz) - B1XYZ(u, m, ix-1,iy,iz));
      MRC_F3(u_delta[1], _B1X + m, ix,iy,iz) = 
	fmaxf((B1XYZ(u, m, ix,iy+1,iz) - B1XYZ(u, m, ix,iy  ,iz)) *
	      (B1XYZ(u, m, ix,iy  ,iz) - B1XYZ(u, m, ix,iy-1,iz)), 0.f) / 
	(B1XYZ(u, m, ix,iy+1,iz) - B1XYZ(u, m, ix,iy-1,iz));
      MRC_F3(u_delta[2], _B1X + m, ix,iy,iz) = 
	fmaxf((B1XYZ(u, m, ix,iy,iz+1) - B1XYZ(u, m, ix,iy,iz  )) *
	      (B1XYZ(u, m, ix,iy,iz  ) - B1XYZ(u, m, ix,iy,iz-1)), 0.f) /
	(B1XYZ(u, m, ix,iy,iz+1) - B1XYZ(u, m, ix,iy,iz-1));
    // FIXME, need to make sure NaN -> 0
      if (!isfinite(MRC_F3(u_delta[0], _B1X + m, ix,iy,iz))) MRC_F3(u_delta[0], _B1X + m, ix,iy,iz) = 0.f;
      if (!isfinite(MRC_F3(u_delta[1], _B1X + m, ix,iy,iz))) MRC_F3(u_delta[1], _B1X + m, ix,iy,iz) = 0.f;
      if (!isfinite(MRC_F3(u_delta[2], _B1X + m, ix,iy,iz))) MRC_F3(u_delta[2], _B1X + m, ix,iy,iz) = 0.f;
    }
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// limit_2
//
// MinMod limiter (Roe, 1986)
// return the smallest if all are positive 
// return the largest if all are negative
// return zero if they are not all postivie or all negative  
// minmod(a,b) = 0.5[sign(a)+sign(b)]*min(|a|,|b|) eq (41) 

static void __unused
limit_2(struct mrc_fld *u_delta[3], struct mrc_fld *u)
{
  int dind[3] = {0, 0, 0};  
  for (int i = 0; i < 3; i++) {    
    dind[i] = 1; 
    for (int m = 0; m <= _B1Z; m++) {    
      mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
	float rl =  ( MRC_F3(u, m, ix+dind[0], iy+dind[1], iz+dind[2])  - MRC_F3(u, m, ix, iy, iz));
	float rr =  ( MRC_F3(u, m, ix, iy, iz)  - MRC_F3(u, m, ix-dind[0], iy-dind[1], iz-dind[2]));
	if ( rl*rr > 0 ) {
	  if ( fabsf(rr)<fabsf(rl) ) {
	    MRC_F3(u_delta[i], m, ix,iy,iz) = 0.5*rr ;            
	  } else { 
	    MRC_F3(u_delta[i], m, ix,iy,iz) = 0.5*rl ;    
	  }
	} else { 
	  MRC_F3(u_delta[i], m, ix,iy,iz) = 0.0; 
	}
      } mrc_fld_foreach_end;      
    }
    dind[i] = 0;
  }
}

// ----------------------------------------------------------------------
// limit_3
//
// MONCEN limiter (Van Leer)

static void __unused
limit_3(struct mrc_fld *u_delta[3], struct mrc_fld *u)
{
  int dind[3] = {0, 0, 0};
  
  for (int i = 0; i < 3; i++) {    
    dind[i] = 1; 
    for (int m = 0; m <= _B1Z; m++) {    
      mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
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
}

// ----------------------------------------------------------------------
// limit_4
//

static void __unused
limit_4(struct mrc_fld *u_delta[3], struct mrc_fld *u)
{
   // generalised minmod limiter with parameter(Van Leer 1979)
  int dind[3] = {0, 0, 0};
  float theta = 1.0; 
  for (int i = 0; i < 3; i++) {    
    dind[i] = 1; 
    for (int m = 0; m <= _B1Z; m++) {    
      mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
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
}

// ----------------------------------------------------------------------
// calc_u_delta

static void
calc_u_delta(struct mrc_fld *_u_delta[3], struct mrc_fld *_u)
{
  struct mrc_fld *u = mrc_fld_get_as(_u, "mhd_fc_float");
  struct mrc_fld *u_delta[3];
  for (int f = 0; f < 3; f++) {
    u_delta[f] = mrc_fld_get_as(_u_delta[f], "float");
  }

  // --- find u_delta
#if LMTR == 1
  limit_1(u_delta, u);
#elif LMTR == 2
  limit_2(u_delta, u);
#elif LMTR == 3
  limit_3(u_delta, u);
#elif LMTR == 4
  limit_4(u_delta, u);
#endif

  mrc_fld_put_as(u, _u);
  for (int f = 0; f < 3; f++) {
    mrc_fld_put_as(u_delta[f], _u_delta[f]);
  }
}

// ----------------------------------------------------------------------
// calc_u_pm

static void
calc_u_pm(struct ggcm_mhd *mhd, struct mrc_fld *_u_p[3], struct mrc_fld *_u_m[3],
	  struct mrc_fld *_E_p[3], struct mrc_fld *_E_m[3],
	  struct mrc_fld *_u, struct mrc_fld *_u_delta[3])
{
  float d_i = mhd->par.d_i;
  float eta = mhd->par.diffco;

  struct mrc_fld *u = mrc_fld_get_as(_u, "mhd_fc_float");
  struct mrc_fld *u_delta[3], *u_p[3], *u_m[3], *E_p[3], *E_m[3];
  for (int f = 0; f < 3; f++) {
    u_delta[f] = mrc_fld_get_as(_u_delta[f], "float");
    u_p[f] = mrc_fld_get_as(_u_p[f], "float");
    u_m[f] = mrc_fld_get_as(_u_m[f], "float");
    E_p[f] = mrc_fld_get_as(_E_p[f], "float");
    E_m[f] = mrc_fld_get_as(_E_m[f], "float");
  }

  // Reconstruction    UijkE  = u_ijk + (dxu_)ijk    UijkW = u_ijk - (dxu_)ijk
  mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
    for (int m = 0; m <= _UU1; m++) {
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

      // reconstruction for B compoents. e.g. i E-W location, Bx is already in correct locatio n
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

  // find currents at cell faces
  float *bdx3 = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD3);
  float *bdy3 = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD3);
  float *bdz3 = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD3);

  for (int i = 0; i < 3; i++) {    
    mrc_fld_foreach(u, ix,iy,iz, 1, 1) {	
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
      MRC_F3(E_p[i], 0, ix,iy,iz) = 
	-((MRC_F3(u_p[i], _RV1Y, ix,iy,iz) * MRC_F3(u_p[i], _B1Z, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz))- 	
	  (MRC_F3(u_p[i], _B1Y, ix,iy,iz) * MRC_F3(u_p[i], _RV1Z, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz)))+
	  eta *  MRC_F3(u_p[i],_JX, ix,iy,iz)
	+ d_i * ((MRC_F3(u_p[i], _JY, ix,iy,iz) * MRC_F3(u_p[i], _B1Z, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz))- 	
		 (MRC_F3(u_p[i], _B1Y, ix,iy,iz) * MRC_F3(u_p[i], _JZ, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz)));  
  
      MRC_F3(E_p[i], 1, ix,iy,iz) = 
	((MRC_F3(u_p[i], _RV1X, ix,iy,iz) * MRC_F3(u_p[i], _B1Z, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz))- 	
	 (MRC_F3(u_p[i], _B1X, ix,iy,iz) * MRC_F3(u_p[i], _RV1Z, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz)))+
	 eta * MRC_F3(u_p[i],_JY, ix,iy,iz)
	- d_i * ((MRC_F3(u_p[i], _JX, ix,iy,iz) * MRC_F3(u_p[i], _B1Z, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz))- 	
		 (MRC_F3(u_p[i], _B1X, ix,iy,iz) * MRC_F3(u_p[i], _JZ, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz)));
      
      MRC_F3(E_p[i], 2, ix,iy,iz) = 
	-((MRC_F3(u_p[i], _RV1X, ix,iy,iz) * MRC_F3(u_p[i], _B1Y, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz))- 	
	  (MRC_F3(u_p[i], _B1X, ix,iy,iz) * MRC_F3(u_p[i], _RV1Y, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz)))+
	  eta * MRC_F3(u_p[i],_JZ, ix,iy,iz)
        + d_i * ((MRC_F3(u_p[i], _JX, ix,iy,iz) * MRC_F3(u_p[i], _B1Y, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz))- 	
		     (MRC_F3(u_p[i], _B1X, ix,iy,iz) * MRC_F3(u_p[i], _JY, ix,iy,iz) / MRC_F3(u_p[i], _RR1, ix,iy,iz)));

	
      // _m 
      MRC_F3(E_m[i], 0, ix,iy,iz) =  
	-((MRC_F3(u_m[i], _RV1Y, ix,iy,iz) * MRC_F3(u_m[i], _B1Z, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz))- 	
	  (MRC_F3(u_m[i], _B1Y, ix,iy,iz) * MRC_F3(u_m[i], _RV1Z, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz)))+
   	  eta * MRC_F3(u_m[i],_JX, ix,iy,iz)
	+ d_i * ((MRC_F3(u_m[i], _JY, ix,iy,iz) * MRC_F3(u_m[i], _B1Z, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz))- 	
		 (MRC_F3(u_m[i], _B1Y, ix,iy,iz) * MRC_F3(u_m[i], _JZ, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz)));
      
      MRC_F3(E_m[i], 1, ix,iy,iz) = 
	((MRC_F3(u_m[i], _RV1X, ix,iy,iz) * MRC_F3(u_m[i], _B1Z, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz))- 	
	 (MRC_F3(u_m[i], _B1X, ix,iy,iz) * MRC_F3(u_m[i], _RV1Z, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz)))+
	 eta * MRC_F3(u_m[i],_JY, ix,iy,iz)
	- d_i * ((MRC_F3(u_m[i], _JX, ix,iy,iz) * MRC_F3(u_m[i], _B1Z, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz))- 	
		 (MRC_F3(u_m[i], _B1X, ix,iy,iz) * MRC_F3(u_m[i], _JZ, ix,iy,iz) / MRC_F3(u_m[i], _RR1, ix,iy,iz)));
      
      MRC_F3(E_m[i], 2, ix,iy,iz) =  
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
    mrc_fld_put_as(E_p[f], _E_p[f]);
    mrc_fld_put_as(E_m[f], _E_m[f]);
  }
}

// ----------------------------------------------------------------------
// calc_fluxes

static void
calc_fluxes(struct ggcm_mhd *mhd, struct mrc_fld *_flux[3],
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
    float ppp = (gamma - 1.f) *
      (MRC_F3(u_p[0], _UU1, ix-1,iy,iz) - .5f * rhoip * (sqr(MRC_F3(u_p[0], _RV1X, ix-1,iy,iz)) +
						     sqr(MRC_F3(u_p[0], _RV1Y, ix-1,iy,iz)) +
						     sqr(MRC_F3(u_p[0], _RV1Z, ix-1,iy,iz)))- 
       .5f * (sqr(MRC_F3(u_p[0], _B1X, ix-1,iy,iz)) +
	      sqr(MRC_F3(u_p[0], _B1Y, ix-1,iy,iz)) +
	      sqr(MRC_F3(u_p[0], _B1Z, ix-1,iy,iz))));
       
    float ppm = (gamma - 1.f) *
      (MRC_F3(u_m[0], _UU1, ix,iy,iz) - .5f * rhoim * (sqr(MRC_F3(u_m[0], _RV1X, ix,iy,iz)) +
						       sqr(MRC_F3(u_m[0], _RV1Y, ix,iy,iz)) +
						       sqr(MRC_F3(u_m[0], _RV1Z, ix,iy,iz)))-       
       .5f * (sqr(MRC_F3(u_m[0], _B1X, ix,iy,iz)) +
	      sqr(MRC_F3(u_m[0], _B1Y, ix,iy,iz)) +
	      sqr(MRC_F3(u_m[0], _B1Z, ix,iy,iz))));


    if (ppm < 0) { 
      ppm = 0; 
    }
    if (ppp < 0) { 
      ppp = 0; 
    } 
    
    float csp = sqrtf((gamma * ppp) / (MRC_F3(u_p[0], _RR1, ix-1,iy,iz)));
    float csm = sqrtf((gamma * ppm) / (MRC_F3(u_m[0], _RR1, ix,iy,iz)));   
       
    float cAp = sqrtf( (sqr(MRC_F3(u_p[0], _B1X, ix-1,iy,iz))+
			sqr(MRC_F3(u_p[0], _B1Y, ix-1,iy,iz))+
			sqr(MRC_F3(u_p[0], _B1Z, ix-1,iy,iz)))/ MRC_F3( u_p[0], _RR1, ix-1,iy,iz) );

    float cAm = sqrtf( (sqr(MRC_F3(u_m[0], _B1X, ix,iy,iz))+
			sqr(MRC_F3(u_m[0], _B1Y, ix,iy,iz))+
			sqr(MRC_F3(u_m[0], _B1Z, ix,iy,iz)))/ MRC_F3(u_m[0], _RR1, ix,iy,iz) );
 						      
    float tmpp = sqr(csp) + sqr(cAp);
    float cfp = sqrtf( 0.5 * ( tmpp + sqrtf( sqr( sqr(cAp) + sqr(csp) )
					     - (4. * mpermi * sqr(csp * MRC_F3(u, _B1X, ix-1,iy,iz)) /  
						MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) )) );      
    float tmpm = sqr(csm) + sqr(cAm);
    float cfm = sqrtf( 0.5 * ( tmpm + sqrtf( sqr( sqr(cAm) + sqr(csm) )
					     - (4.* mpermi * sqr(csm * MRC_F3(u, _B1X, ix,iy,iz)) /  
						MRC_F3(u_m[0], _RR1, ix,iy,iz))  )) );
    
    /*
    ap = fmaxf(fmaxf((MRC_F3(u_p[0], _RV1X, ix-1,iy,iz) / MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) + cfp,
		     (MRC_F3(u_m[0], _RV1X, ix,iy,iz) / MRC_F3(u_m[0], _RR1, ix,iy,iz)) + cfm), 0.f);

    am = fminf(fminf( (MRC_F3(u_p[0], _RV1X, ix-1,iy,iz) / MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) - cfp,
		     (MRC_F3(u_m[0], _RV1X, ix,iy,iz) / MRC_F3(u_m[0], _RR1,ix,iy,iz)) - cfm), 0.f);
    */

    // cw =  ( |B|*pi )  /  ( e * rho * delta x ) 
    float cwp = d_i * cAp * sqrtf(1./MRC_F3(u_p[0], _RR1, ix-1,iy,iz)) * bdx1[ix+2]  ;
    float cwm = d_i * cAm * sqrtf(1./MRC_F3(u_m[0], _RR1, ix,iy,iz)) * bdx1[ix+2]  ; 

#if incws == 0
    cwp = 0;
    cwm = 0;
#endif

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


    if (ppm < 0) { 
      ppm = 0; 
    }
    if (ppp < 0) { 
      ppp = 0; 
    } 

    
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
    
    float bp = fmaxf(fmaxf( (MRC_F3(u_p[1], _RV1Y, ix,iy-1,iz) / MRC_F3(u_p[1], _RR1, ix,iy-1,iz)) + cfp + cwp,
			    (MRC_F3(u_m[1], _RV1Y, ix,iy,iz) / MRC_F3(u_m[1], _RR1, ix,iy,iz)) + cfm + cwm), 0.f);
    float bm = fminf(fminf( (MRC_F3(u_p[1], _RV1Y, ix,iy-1,iz) / MRC_F3(u_p[1], _RR1, ix,iy-1,iz)) - cfp - cwp,
			    (MRC_F3(u_m[1], _RV1Y, ix,iy,iz) / MRC_F3(u_m[1], _RR1, ix,iy,iz)) - cfm - cwm), 0.f);
    
#if KT == 1
    bp = fmaxf(bp,-bm);
    bm=-bp;
#endif    

   // Coeffiecents cp and cm 

    rhoip = 1.f / MRC_F3(u_p[2], _RR1, ix,iy,iz);
    rhoim = 1.f / MRC_F3(u_m[2], _RR1, ix,iy,iz-1); 


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
       
    float tmpbe =  .5f * (sqr(MRC_F3(u_p[2], _B1X, ix,iy,iz)) +
	      sqr(MRC_F3(u_p[2], _B1Y, ix,iy,iz)) +
			  sqr(MRC_F3(u_p[2], _B1Z, ix,iy,iz)));

    float tmpke = rhoip ;// 0.5f * rhoip * (sqr(MRC_F3(u_p[2], _RV1X, ix,iy,iz)) +
    //sqr(MRC_F3(u_p[2], _RV1Y, ix,iy,iz)) +
    //				   sqr(MRC_F3(u_p[2], _RV1Z, ix,iy,iz)));  
    float tmpee = MRC_F3(u_p[2], _UU1, ix,iy,iz); 

    csp = sqrtf((gamma * ppp) * rhoip);
    csm = sqrtf((gamma * ppm) * rhoim);
   			
    cAp = sqrtf( (sqr(MRC_F3(u_p[2], _B1X, ix,iy,iz-1))+
		  sqr(MRC_F3(u_p[2], _B1Y, ix,iy,iz-1))+
		  sqr(MRC_F3(u_p[2], _B1Z, ix,iy,iz-1))) * rhoip);// MRC_F3(u_p[2], _RR1, ix,iy,iz-1) );
    
    cAm = sqrtf( (sqr(MRC_F3(u_m[2], _B1X, ix,iy,iz))+
		  sqr(MRC_F3(u_m[2], _B1Y, ix,iy,iz))+
		  sqr(MRC_F3(u_m[2], _B1Z, ix,iy,iz))) * rhoim); // MRC_F3(u_m[2], _RR1, ix,iy,iz) );

    if (ppm < 0) { 
      ppm = 0; 
    }
    if (ppp < 0) { 
      ppp = 0; 
    } 

    tmpp = sqr(csp) + sqr(cAp);
    cfp = sqrtf( 0.5f * (tmpp + sqrtf( sqr( sqr(cAp) + sqr(csp) ) - 
				       (4. * mpermi * sqr(csp * MRC_F3(u_p[2], _B1Z, ix,iy,iz-1)) *  
					rhoip )) ));      
    if (!isfinite(cfp)) {
      mprintf("ix %d %d %d cfp = csp %g cAp %g rr %g ppp %g ppm %g tmpbe %g tmpke %g tmpee %g \n", ix,iy,iz, csp, cAp, MRC_F3(u_p[2], _RR1, ix,iy,iz-1), ppp, ppm, tmpbe, tmpke, tmpee);
    }

    tmpm = sqr(csm) + sqr(cAm);
    cfm = sqrtf( 0.5f * (tmpm + sqrtf( sqr( sqr(cAm) + sqr(csm) ) -
				       (4.* mpermi * sqr(csm * MRC_F3(u_m[2], _B1Z, ix,iy,iz)) *  
					rhoim)) ));

    /*
    cp = fmaxf(fmaxf( (MRC_F3(u_p[2], _RV1Z, ix,iy,iz-1) / MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) + cfp,
 		      (MRC_F3(u_m[2], _RV1Z, ix,iy,iz) / MRC_F3(u_m[2], _RR1, ix,iy,iz)) + cfm), 0.f);
    cm = fminf(fminf( (MRC_F3(u_p[2], _RV1Z, ix,iy,iz-1) / MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) - cfp,
		      (MRC_F3(u_m[2], _RV1Z, ix,iy,iz) / MRC_F3(u_m[2], _RR1, ix,iy,iz)) - cfm), 0.f);
    */

    
    //cwp = d_i * cAp * M_PI * sqrtf(MRC_F3(u_p[2], _RR1, ix,iy,iz-1)) * bdz1[iz+2]  ;
    //cwm = d_i * cAm * M_PI * sqrtf(MRC_F3(u_m[2], _RR1, ix,iy,iz)) * bdz1[iz+2]  ; 

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
    //    bp = 1e-3; bm = -1e-3;
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
    //    cp = 1e-3; cm = -1e-3;
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
#if 1
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
	assert(isfinite(FLUX(flux, 0, m, ix,iy,iz)));
	assert(isfinite(FLUX(flux, 1, m, ix,iy,iz)));
	if (!isfinite(FLUX(flux, 2, m, ix,iy,iz))) {
	  mprintf("ix %d %d %d cp %g cm %g\n", ix,iy,iz, cp, cm);
	}
	assert(isfinite(FLUX(flux, 2, m, ix,iy,iz)));
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

// ----------------------------------------------------------------------
// calc_cweno_fluxes
//
// calculates CWENO fluxes on faces in flux_E, from the original state
// vector u (which is cell centered / on the Yee grid)
// flux[0-4] are the fluid vars
// flux[5-7] are E-field "fluxes" (not B-field!)

static void
calc_cweno_fluxes(struct ggcm_mhd *mhd, struct mrc_fld *_flux[3],
		  struct mrc_fld *_u)
{
  //initialize cell surface center variables			    
  struct mrc_fld *_u_delta[3], *_u_p[3], *_u_m[3], *_E_p[3], *_E_m[3];
  struct mrc_fld *_flux_p[3], *_flux_m[3];
  for (int f = 0; f < 3; f++) {
    _u_delta[f] = ggcm_mhd_get_fields(mhd, "u_delta", _B1Z + 1);
    _u_p[f] = ggcm_mhd_get_fields(mhd, "u_p", _JZ + 1);
    _u_m[f] = ggcm_mhd_get_fields(mhd, "u_m", _JZ + 1);
    _E_p[f] = ggcm_mhd_get_fields(mhd, "E_p", 3);
    _E_m[f] = ggcm_mhd_get_fields(mhd, "E_m", 3);
    _flux_p[f] = ggcm_mhd_get_fields(mhd, "flux_p", 8);
    _flux_m[f] = ggcm_mhd_get_fields(mhd, "flux_m", 8);
  }
  
  ggcm_mhd_fill_ghosts(mhd, _u, 0, mhd->time);

  calc_u_delta(_u_delta, _u);

  /* ggcm_mhd_fill_ghosts(mhd, _u_delta[0], 0, mhd->time); */
  /* ggcm_mhd_fill_ghosts(mhd, _u_delta[1], 0, mhd->time); */
  /* ggcm_mhd_fill_ghosts(mhd, _u_delta[2], 0, mhd->time); */

  calc_u_pm(mhd, _u_p, _u_m, _E_p, _E_m, _u, _u_delta);
  
#ifdef DEBUG
  {
    static struct ggcm_mhd_diag *diag;
    static int cnt;
    if (!diag) {
      diag = ggcm_mhd_diag_create(ggcm_mhd_comm(mhd));
      ggcm_mhd_diag_set_param_obj(diag, "mhd", mhd);
      ggcm_mhd_diag_set_param_string(diag, "run", "u_p0");
      ggcm_mhd_diag_set_param_string(diag, "fields", "rr1:rv1:uu1:b1");
      ggcm_mhd_diag_setup(diag);
      ggcm_mhd_diag_view(diag);
    }
    ggcm_mhd_diag_run_now(diag, _u_p[0], DIAG_TYPE_3D, cnt++);
  }
  {
    static struct ggcm_mhd_diag *diag;
    static int cnt;
    if (!diag) {
      diag = ggcm_mhd_diag_create(ggcm_mhd_comm(mhd));
      ggcm_mhd_diag_set_param_obj(diag, "mhd", mhd);
      ggcm_mhd_diag_set_param_string(diag, "run", "u_m0");
      ggcm_mhd_diag_set_param_string(diag, "fields", "rr1:rv1:uu1:b1");
      ggcm_mhd_diag_setup(diag);
      ggcm_mhd_diag_view(diag);
    }
    ggcm_mhd_diag_run_now(diag, _u_m[0], DIAG_TYPE_3D, cnt++);
  }
#endif

  // calculate fluxes per face (small f's) using reconstructed 
  // variables U^N(SWETB) and B^N(SWETB) = (Bx,By,Bz)^N(SEWTB)
  for (int f = 0; f < 3; f++) {
    calc_fluxes_per_face(_flux_p, mhd, _u_p[f], f);
    calc_fluxes_per_face(_flux_m, mhd, _u_m[f], f);
  }
  
  calc_fluxes(mhd, _flux, _flux_p, _flux_m, _u, _u_p, _u_m, _E_p, _E_m);

  for (int f = 0; f < 3; f++) {
    mrc_fld_destroy(_flux_p[f]);
    mrc_fld_destroy(_flux_m[f]);
    mrc_fld_destroy(_u_delta[f]);
    mrc_fld_destroy(_u_p[f]);
    mrc_fld_destroy(_u_m[f]);
    mrc_fld_destroy(_E_p[f]);
    mrc_fld_destroy(_E_m[f]);
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
calc_fct_rhs(struct ggcm_mhd *mhd, struct mrc_fld *_rhs, struct mrc_fld *_flux[3])
{  
  // compute edge centered electric fields by interpolation (iv)  here tmp_fld are edge centered 
  // so that e.g.  MRC_F3(tmp_fld, 0, 0, 0, 0) is E_x 0,-1/2,-1/2  
  //         i.e.  MRC_F3(tmp_fld, 0, ix,iy,iz) is E_x ix,iy-1/2,iz-1/2   
  //           and MRC_F3(tmp_fld, 1, 0, 0, 0) is E_y -1/2,0,-1/2  etc etc 
                       

  struct mrc_ddc *ddc = mrc_domain_get_ddc(mhd->domain);
  mrc_ddc_set_param_int(ddc, "max_n_fields", 3);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct mrc_fld *_E_ec = ggcm_mhd_get_fields(mhd, "E_ec", 3);

  struct mrc_fld *E_ec = mrc_fld_get_as(_E_ec, "float");
  struct mrc_fld *flux[3];
  for (int f = 0; f < 3; f++) {
    flux[f] = mrc_fld_get_as(_flux[f], "float");
  }

  //initialize cell edge center Electric field structure      
  mrc_fld_foreach(E_ec, ix,iy,iz, 1, 1) { 
    MRC_F3(E_ec, 0, ix,iy,iz) = .25f*(- FLUX(flux, 1, _EZ, ix  ,iy  ,iz  )
				      - FLUX(flux, 1, _EZ, ix  ,iy  ,iz-1) 
				      + FLUX(flux, 2, _EY, ix  ,iy  ,iz  )
				      + FLUX(flux, 2, _EY, ix  ,iy-1,iz  ));    
    MRC_F3(E_ec, 1, ix,iy,iz) = .25f*(- FLUX(flux, 2, _EX, ix  ,iy  ,iz  )
				      - FLUX(flux, 2, _EX, ix-1,iy  ,iz  )
				      + FLUX(flux, 0, _EZ, ix  ,iy  ,iz  )
				      + FLUX(flux, 0, _EZ, ix  ,iy  ,iz-1));
    MRC_F3(E_ec, 2, ix,iy,iz) = .25f*(- FLUX(flux, 0, _EY, ix  ,iy  ,iz  )
				      - FLUX(flux, 0, _EY, ix  ,iy-1,iz  )
				      + FLUX(flux, 1, _EX, ix  ,iy  ,iz  )
				      + FLUX(flux, 1, _EX, ix-1,iy  ,iz  ));    
  } mrc_fld_foreach_end;
  
  for (int f = 0; f < 3; f++) {
    mrc_fld_put_as(flux[f], _flux[f]);
  }
  mrc_fld_put_as(E_ec, _E_ec);

  //  fill_ghost_fld(mhd, E_ec);

  struct mrc_fld *rhs = mrc_fld_get_as(_rhs, "mhd_fc_float");
  E_ec = mrc_fld_get_as(_E_ec, "float");

  mrc_fld_foreach(rhs, ix, iy,  iz, 1, 1) {
    B1X(rhs, ix, iy, iz) =  
      (-((MRC_F3(E_ec, 2, ix, iy+1, iz) - MRC_F3(E_ec, 2, ix, iy, iz)) /
	 (.5f*( MRC_CRDY(crds, iy+1) - MRC_CRDY(crds, iy-1)))) 
       +((MRC_F3(E_ec, 1, ix, iy, iz+1) - MRC_F3(E_ec, 1, ix, iy, iz)) /
	 (.5f*( MRC_CRDZ(crds, iz+1) - MRC_CRDZ(crds, iz-1))))); 
    
    B1Y(rhs, ix, iy, iz) = 
      (-((MRC_F3(E_ec, 0, ix, iy, iz+1) - MRC_F3(E_ec, 0, ix, iy, iz)) /
	 (.5f*( MRC_CRD(crds, 2, iz+1) - MRC_CRD(crds, 2, iz-1))))
       +((MRC_F3(E_ec, 2, ix+1, iy, iz) - MRC_F3(E_ec, 2, ix, iy, iz)) /
	 (.5f*(MRC_CRD(crds, 0, ix+1) - MRC_CRD(crds, 0, ix-1))))); 
    
    B1Z(rhs, ix, iy, iz) = 
      (-((MRC_F3( E_ec, 1, ix+1, iy, iz) - MRC_F3(E_ec, 1, ix, iy, iz)) /
	 (.5f*( MRC_CRD(crds, 0, ix+1) - MRC_CRD(crds, 0, ix-1))))  
       +((MRC_F3( E_ec, 0, ix, iy+1, iz) - MRC_F3(E_ec, 0, ix, iy, iz)) /
	 (.5f*( MRC_CRD(crds, 1, iy+1) - MRC_CRD(crds, 1, iy-1))))); 
  } mrc_fld_foreach_end;

  mrc_fld_put_as(E_ec, _E_ec);
  mrc_fld_put_as(rhs, _rhs);

  mrc_fld_destroy(_E_ec);
}

static void
ggcm_mhd_step_cweno_calc_rhs(struct ggcm_mhd_step *step, struct mrc_fld *rhs,
			     struct mrc_fld *fld)
{
  struct ggcm_mhd *mhd = step->mhd;  
  
  //fill_ghost_fld(mhd, fld);

#ifdef DEBUG
  {
    static struct ggcm_mhd_diag *diag;
    static int cnt;
    if (!diag) {
      diag = ggcm_mhd_diag_create(ggcm_mhd_comm(mhd));
      ggcm_mhd_diag_set_param_obj(diag, "mhd", mhd);
      ggcm_mhd_diag_set_param_string(diag, "run", "fld");
      ggcm_mhd_diag_set_param_string(diag, "fields", "rr1:rv1:uu1:b1:divb:pp_full:b");
      ggcm_mhd_diag_setup(diag);
      ggcm_mhd_diag_view(diag);
    }
    ggcm_mhd_diag_run_now(diag, fld, DIAG_TYPE_3D, cnt++);
  }
#endif

  struct mrc_fld *flux[3];
  for (int f = 0; f < 3; f++) {
    flux[f] = ggcm_mhd_get_fields(mhd, "flux", 8);
  }

  calc_cweno_fluxes(mhd, flux, fld);
  calc_neg_divg(mhd, rhs, flux);
  calc_fct_rhs(mhd, rhs, flux);

  struct mrc_fld *r = mrc_fld_get_as(rhs, "mhd_fc_float");
  assert(mhd->fld->_data_type == MRC_NT_FLOAT);
  struct mrc_fld *f = mrc_fld_get_as(mhd->fld, mrc_fld_type(mhd->fld));

  mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
    for (int m = 0; m < 8; m++) {
      MRC_F3(r, m, ix,iy,iz) *= MRC_F3(f, _YMASK, ix,iy,iz);
    }
  } mrc_fld_foreach_end;

  mrc_fld_put_as(r, rhs);
  mrc_fld_put_as(f, mhd->fld);


#ifdef DEBUG
  {
    static struct ggcm_mhd_diag *diag;
    static int cnt;
    if (!diag) {
      diag = ggcm_mhd_diag_create(ggcm_mhd_comm(mhd));
      ggcm_mhd_diag_set_param_obj(diag, "mhd", mhd);
      ggcm_mhd_diag_set_param_string(diag, "run", "rhs");
      ggcm_mhd_diag_set_param_string(diag, "fields", "rr1:rv1:uu1:b1:divb");
      ggcm_mhd_diag_setup(diag);
      ggcm_mhd_diag_view(diag);
    }
    ggcm_mhd_diag_run_now(diag, rhs, DIAG_TYPE_3D, cnt++);
  }
#endif

  for (int f = 0; f < 3; f++) {
    mrc_fld_destroy(flux[f]);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "cweno"

struct ggcm_mhd_step_ops ggcm_mhd_step_cweno_ops = {
  .name        = "cweno",
  .calc_rhs    = ggcm_mhd_step_cweno_calc_rhs,
};

