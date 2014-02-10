#include "ggcm_mhd_step_cweno_private.h" 

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>

// ----------------------------------------------------------------------
// limit_0
//

void 
limit_0(struct mrc_fld *u_delta[3], struct mrc_fld *u)
{
  //None 
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
      if (!isfinite(MRC_F3(u_delta[0], _B1X + m, ix,iy,iz))) {
	MRC_F3(u_delta[0], _B1X + m, ix,iy,iz) = 0.f;
      }
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
    mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
	for (int m = 0; m <= _UU1; m++) {    
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
	}
	for (int m = 0; m <= 2; m++) {    
	  float rl =  ( B1XYZ(u, m, ix+dind[0], iy+dind[1], iz+dind[2])  - B1XYZ(u, m, ix, iy, iz));
	  float rr =  ( B1XYZ(u, m, ix, iy, iz)  - B1XYZ(u, m, ix-dind[0], iy-dind[1], iz-dind[2]));
	  if ( rl*rr > 0 ) {
	    if ( fabsf(rr)<fabsf(rl) ) {
	      MRC_F3(u_delta[i], _B1X+m, ix,iy,iz) = 0.5*rr ;            
	    } else { 
	      MRC_F3(u_delta[i], _B1X+m, ix,iy,iz) = 0.5*rl ;    
	    }
	  } else { 
	    MRC_F3(u_delta[i], _B1X+m, ix,iy,iz) = 0.0; 
	  }
	}
      } mrc_fld_foreach_end;      
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

    mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
      for (int m = 0; m <= _UU1; m++) {    
	float rl =  2.0*( MRC_F3(u, m, ix+dind[0], iy+dind[1], iz+dind[2])  - MRC_F3(u, m, ix, iy, iz));
	float rr =  2.0*( MRC_F3(u, m, ix, iy, iz)  - MRC_F3(u, m, ix-dind[0], iy-dind[1], iz-dind[2]));
	float cen = 0.5*( MRC_F3(u, m, ix+dind[0], iy+dind[1], iz+dind[2])  - 
			  MRC_F3(u, m, ix-dind[0], iy-dind[1], iz-dind[2]) ) ; 	
	if ( (cen*rl*rr) > 0) {
	  MRC_F3(u_delta[i], m, ix,iy,iz) = 0.5*sign(cen)*fmin(fabsf(cen),fmin(fabsf(rl),fabsf(rr))) ; 
	} else {
	  MRC_F3(u_delta[i], m, ix,iy,iz) = 0.0 ; 	  
	  //MRC_F3(u_delta[i], m, ix,iy,iz) = 0.0 ; 
	  if (!isfinite(MRC_F3(u_delta[i], m, ix,iy,iz))) MRC_F3(u_delta[i], m, ix,iy,iz) = 0.f;
	}
      }
    } mrc_fld_foreach_end;      

    mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
      for (int m = 0; m <= 3; m++) {    
	float rl =  2.0*( B1XYZ(u, _B1X+m, ix+dind[0], iy+dind[1], iz+dind[2])  - B1XYZ(u, _B1X+m, ix, iy, iz));
	float rr =  2.0*( B1XYZ(u, _B1X+m, ix, iy, iz)  - B1XYZ(u, _B1X+m, ix-dind[0], iy-dind[1], iz-dind[2]));
	float cen = 0.5*((B1XYZ(u, _B1X+m, ix+dind[0], iy+dind[1], iz+dind[2])  - 
			  B1XYZ(u, _B1X+m, ix-dind[0], iy-dind[1], iz-dind[2]))) ; 	
	if ( (cen*rl*rr) > 0) {
	  MRC_F3(u_delta[i], m, ix,iy,iz) = 0.5*sign(cen)*fmin(fabsf(cen),fmin(fabsf(rl),fabsf(rr))) ; 
	} else {
	  MRC_F3(u_delta[i], m, ix,iy,iz) = 0.0 ; 	  
	  if (!isfinite(MRC_F3(u_delta[i], m, ix,iy,iz))) MRC_F3(u_delta[i], m, ix,iy,iz) = 0.f;
	}
      }	
    } mrc_fld_foreach_end;          
    dind[i] = 0;
  }

}

// ----------------------------------------------------------------------
// limit_4
//

void 
limit_4(struct mrc_fld *u_delta[3], struct mrc_fld *u)
{
   // generalised minmod limiter with parameter(Van Leer 1979)
  int dind[3] = {0, 0, 0};
  float theta = 2.0f; 
  for (int i = 0; i < 3; i++) {    
    dind[i] = 1; 
    mrc_fld_foreach(u, ix,iy,iz, 1, 1) {
      for (int m = 0; m <= _UU1; m++) {    
	float rl =  theta *( MRC_F3(u, m, ix+dind[0], iy+dind[1], iz+dind[2])  - MRC_F3(u, m, ix, iy, iz));
	float rr =  theta *( MRC_F3(u, m, ix, iy, iz)  - MRC_F3(u, m, ix-dind[0], iy-dind[1], iz-dind[2]));
	float cen = .5f *( MRC_F3(u, m, ix+dind[0], iy+dind[1], iz+dind[2])  - 
			   MRC_F3(u, m, ix-dind[0], iy-dind[1], iz-dind[2]) ) ; 		  
	if ( (cen*rl*rr) > 0) {
	  MRC_F3(u_delta[i], m, ix,iy,iz) = 0.5*sign(cen)*fmin(fabsf(cen),fmin(fabsf(rl),fabsf(rr))) ; 
	} else {
	  MRC_F3(u_delta[i], m, ix,iy,iz) = 0.0 ; 
	}	  
	if (!isfinite(MRC_F3(u_delta[i], m, ix,iy,iz))) MRC_F3(u_delta[i], m, ix,iy,iz) = 0.f;	 
      }      
      for (int m = 0; m <3; m++) {    
	float rl =  theta *( B1XYZ(u, _B1X +m, ix+dind[0], iy+dind[1], iz+dind[2])  - B1XYZ(u, _B1X +m, ix, iy, iz));
	float rr =  theta *( B1XYZ(u, _B1X +m, ix, iy, iz)  - B1XYZ(u, _B1X +m, ix-dind[0], iy-dind[1], iz-dind[2]));
	float cen = .5f *( B1XYZ(u, _B1X +m, ix+dind[0], iy+dind[1], iz+dind[2])  - 
			   B1XYZ(u, _B1X +m, ix-dind[0], iy-dind[1], iz-dind[2]) ) ; 		  
	if ( (cen*rl*rr) > 0) {	  
	  MRC_F3(u_delta[i], _B1X +m, ix,iy,iz) = 0.5*sign(cen)*fmin(fabsf(cen),fmin(fabsf(rl),fabsf(rr))) ; 
	} else {
	  MRC_F3(u_delta[i], _B1X +m, ix,iy,iz) = 0.0 ; 
	}	  
	if (!isfinite(MRC_F3(u_delta[i], _B1X +m, ix,iy,iz))) MRC_F3(u_delta[i], _B1X +m, ix,iy,iz) = 0.f;	 
      }
    } mrc_fld_foreach_end;        
    dind[i] = 0;
  }
}


// ----------------------------------------------------------------------
// calc_u_delta

void
calc_u_delta(struct mrc_fld *u_delta[3], struct mrc_fld *u, struct ggcm_mhd *mhd)
{
  // --- find u_delta
#if LMTR == 0 
  limit_0(u_delta,u); 
#elif LMTR == 1
  limit_1(u_delta, u);
#elif LMTR == 2
  limit_2(u_delta, u);
#elif LMTR == 3
  limit_3(u_delta, u);
#elif LMTR == 4
  limit_4(u_delta, u);
#endif
  
}
