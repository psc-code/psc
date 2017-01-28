#include "ggcm_mhd_step_cweno_private.h" 

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>

// ----------------------------------------------------------------------
// calc_u_pm  
//
// calculate point values of variables at cell surface centers using 
// linear reconstrcution with TVD slope limiters 
// (Ziegler 2004 section 3.2) 

void
calc_u_pm(struct ggcm_mhd *mhd, struct mrc_fld *u_p[3], struct mrc_fld *u_m[3],
	  struct mrc_fld *E_p[3], struct mrc_fld *E_m[3],
	  struct mrc_fld *u, struct mrc_fld *u_delta[3])
{
  float d_i = mhd->par.d_i;
  float eta = mhd->par.diffco;
 
  // Reconstruction    UijkE  = u_ijk + (dxu_)ijk    UijkW = u_ijk - (dxu_)ijk
  mrc_fld_foreach(u, ix,iy,iz, 2, 1) {
    for (int m = 0; m <= UU; m++) {
      // defined like this, both u_p and u_m are coplaner when indices are the same    
      MRC_F3(u_p[0], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) + MRC_F3(u_delta[0], m, ix,iy,iz);
      MRC_F3(u_p[1], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) + MRC_F3(u_delta[1], m, ix,iy,iz);
      MRC_F3(u_p[2], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) + MRC_F3(u_delta[2], m, ix,iy,iz);
      MRC_F3(u_m[0], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) - MRC_F3(u_delta[0], m, ix,iy,iz);
      MRC_F3(u_m[1], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) - MRC_F3(u_delta[1], m, ix,iy,iz);
      MRC_F3(u_m[2], m, ix,iy,iz) = MRC_F3(u, m, ix,iy,iz) - MRC_F3(u_delta[2], m, ix,iy,iz);
    }
  } mrc_fld_foreach_end;

  
     mrc_fld_foreach(u, ix,iy,iz, 2, 1) {
    if  (MRC_F3(u_p[0], RR, ix,iy,iz) <= 0.f) {
      MRC_F3(u_p[0], RR, ix,iy,iz) = RMIN;
    }
    if  (MRC_F3(u_p[1], RR, ix,iy,iz) <= 0.f) {
      MRC_F3(u_p[1], RR, ix,iy,iz) = RMIN;
    }
    if  (MRC_F3(u_p[2], RR, ix,iy,iz) <= 0.f) {
      MRC_F3(u_p[2], RR, ix,iy,iz) = RMIN;
    }    
    if (MRC_F3(u_m[0], RR, ix,iy,iz) <= 0.f) { 
      MRC_F3(u_m[0], RR, ix,iy,iz) = RMIN;
    }
    if  (MRC_F3(u_m[1], RR, ix,iy,iz) <= 0.f) {
      MRC_F3(u_m[1], RR, ix,iy,iz) = RMIN;
    }    
    if (MRC_F3(u_m[2], RR, ix,iy,iz) <= 0.f) { 
      MRC_F3(u_m[2], RR, ix,iy,iz) = RMIN;
    }
   
  } mrc_fld_foreach_end;
  

  // calculation of cell surface center averages for B 
  // note that B is staggered here.. so that index 1 for B is 1-1/2 i.e. cell surface for V
  mrc_fld_foreach(u, ix,iy,iz, 2, 1) {
    for (int i = 0; i < 3; i++) {

      int cycl[5]={0,1,2,0,1};
      int dind[3]={0,0,0};
      int ip1= cycl[i+1];
      int ip2= cycl[i+2];

      // reconstruction for B compoents. e.g. i E-W location, Bx is already in correct locatio n
      // but transverse components are first reconstructed to cell edges then averaged to cell surfaces

      // _p  i-> 0:E 1:   
      dind[i]=1; 
      MRC_F3(u_p[i], BX+i, ix,iy,iz) = MRC_F3(u, BX+i, ix+dind[0],iy+dind[1],iz+dind[2]);
      // +1 --> Bxi+1/2y,z because staggered grid 
      dind[i]=0;

      dind[ip1]=1;
      MRC_F3(u_p[i], BX+ip1, ix,iy,iz) =
      	(0.5*(MRC_F3(u, BX+ip1, ix+dind[0],iy+dind[1],iz+dind[2]) +
      	      MRC_F3(u_delta[i], BX+ip1, ix+dind[0],iy+dind[1],iz+dind[2]) +
      	      MRC_F3(u, BX+ip1, ix,iy,iz) + MRC_F3(u_delta[i], BX+ip1, ix,iy,iz)));
      dind[ip1]=0;
       
      dind[ip2]=1;
      MRC_F3(u_p[i], BX+ip2, ix,iy,iz) =
      	(0.5*(MRC_F3(u, BX+ip2, ix+dind[0],iy+dind[1],iz+dind[2] ) +
      	      MRC_F3(u_delta[i], BX+ip2, ix+dind[0],iy+dind[1],iz+dind[2]) +
      	      MRC_F3(u, BX+ip2, ix,iy,iz) + MRC_F3(u_delta[i], BX+ip2, ix,iy,iz)));
      dind[ip2]=0;
      
      // _m
      dind[i]=0; 
      MRC_F3(u_m[i], BX+i, ix,iy,iz) = MRC_F3(u, BX+i, ix,iy,iz);
      //  +0 --> Bxi-1/2y,z because staggered grid
      
      dind[ip1]=1;
      MRC_F3(u_m[i], BX+ip1, ix,iy,iz) = 
	(0.5*(MRC_F3(u, BX+ip1, ix+dind[0],iy+dind[1],iz+dind[2]) -
	      MRC_F3(u_delta[i], BX+ip1, ix+dind[0],iy+dind[1],iz+dind[2]) +
	      MRC_F3(u, BX+ip1, ix,iy,iz) - MRC_F3(u_delta[i], BX+ip1, ix,iy,iz)));
      dind[ip1]=0;

      dind[ip2]=1;
      MRC_F3(u_m[i], BX+ip2, ix,iy,iz) =
	(0.5*(MRC_F3(u, BX+ip2, ix+dind[0],iy+dind[1],iz+dind[2] ) -
	      MRC_F3(u_delta[i], BX+ip2, ix+dind[0],iy+dind[1],iz+dind[2]) +
	      MRC_F3(u, BX+ip2, ix,iy,iz) - MRC_F3(u_delta[i], BX+ip2, ix,iy,iz)));
      dind[ip2]=0;
    }
  } mrc_fld_foreach_end;

  // find currents at cell faces
  float *bdx3 = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD3);
  float *bdy3 = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD3);
  float *bdz3 = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD3);

  for (int i = 0; i < 3; i++) {    
    mrc_fld_foreach(u, ix,iy,iz, 2, 2) {	
      // _p
      MRC_F3(u_p[i],_JX, ix,iy,iz) = 
	0.5*((MRC_F3(u,BZ, ix,iy+1,iz) - MRC_F3(u,BZ, ix,iy-1,iz)) * bdy3[iy] - 
	     (MRC_F3(u,BY, ix,iy,iz+1) - MRC_F3(u,BY, ix,iy,iz-1)) * bdz3[iz]);     

      MRC_F3(u_p[i],_JY, ix,iy,iz) =
	0.5*((MRC_F3(u,BX, ix,iy,iz+1) - MRC_F3(u,BX, ix,iy,iz-1)) * bdz3[iz] -
	     (MRC_F3(u,BZ, ix+1,iy,iz) - MRC_F3(u,BZ, ix-1,iy,iz)) * bdx3[ix]);       

      MRC_F3(u_p[i],_JZ, ix,iy,iz) = 
	0.5*((MRC_F3(u,BY, ix+1,iy,iz) - MRC_F3(u,BY, ix-1,iy,iz)) * bdx3[ix] - 
	     (MRC_F3(u,BX, ix,iy+1,iz) - MRC_F3(u,BX, ix,iy-1,iz)) * bdy3[iy]); 

     
      // _m 
      MRC_F3(u_m[i], _JX, ix,iy,iz) = MRC_F3(u_p[i], _JX, ix-1,iy,iz) ; 
      MRC_F3(u_m[i], _JY, ix,iy,iz) = MRC_F3(u_p[i], _JY, ix,iy-1,iz) ;
      MRC_F3(u_m[i], _JZ, ix,iy,iz) = MRC_F3(u_p[i], _JZ, ix,iy,iz-1) ;      
    } mrc_fld_foreach_end;
  }
  
  // Calculation of cell surface center  averages for E using v and B just calculated 
  //  E^N(SWETB) = -V^N(SWETB) X B^N(SWETB)      + eta*J  [  + di J x B ] 
  mrc_fld_foreach(u, ix, iy, iz, 2, 2) {  
    for (int i = 0; i < 3; i++) {     
      // _p 
      MRC_F3(E_p[i], 0, ix,iy,iz) = 
	-((MRC_F3(u_p[i], RVY, ix,iy,iz) * MRC_F3(u_p[i], BZ, ix,iy,iz) / MRC_F3(u_p[i], RR, ix,iy,iz))- 	
	  (MRC_F3(u_p[i], BY, ix,iy,iz) * MRC_F3(u_p[i], RVZ, ix,iy,iz) / MRC_F3(u_p[i], RR, ix,iy,iz)))+
	  eta *  MRC_F3(u_p[i],_JX, ix,iy,iz)
	+ d_i * ((MRC_F3(u_p[i], _JY, ix,iy,iz) * MRC_F3(u_p[i], BZ, ix,iy,iz) / MRC_F3(u_p[i], RR, ix,iy,iz))- 	
		 (MRC_F3(u_p[i], BY, ix,iy,iz) * MRC_F3(u_p[i], _JZ, ix,iy,iz) / MRC_F3(u_p[i], RR, ix,iy,iz)));  
  
      MRC_F3(E_p[i], 1, ix,iy,iz) = 
	((MRC_F3(u_p[i], RVX, ix,iy,iz) * MRC_F3(u_p[i], BZ, ix,iy,iz) / MRC_F3(u_p[i], RR, ix,iy,iz))- 	
	 (MRC_F3(u_p[i], BX, ix,iy,iz) * MRC_F3(u_p[i], RVZ, ix,iy,iz) / MRC_F3(u_p[i], RR, ix,iy,iz)))+
	 eta * MRC_F3(u_p[i],_JY, ix,iy,iz)
	- d_i * ((MRC_F3(u_p[i], _JX, ix,iy,iz) * MRC_F3(u_p[i], BZ, ix,iy,iz) / MRC_F3(u_p[i], RR, ix,iy,iz))- 	
		 (MRC_F3(u_p[i], BX, ix,iy,iz) * MRC_F3(u_p[i], _JZ, ix,iy,iz) / MRC_F3(u_p[i], RR, ix,iy,iz)));
      
      MRC_F3(E_p[i], 2, ix,iy,iz) = 
	-((MRC_F3(u_p[i], RVX, ix,iy,iz) * MRC_F3(u_p[i], BY, ix,iy,iz) / MRC_F3(u_p[i], RR, ix,iy,iz))- 	
	  (MRC_F3(u_p[i], BX, ix,iy,iz) * MRC_F3(u_p[i], RVY, ix,iy,iz) / MRC_F3(u_p[i], RR, ix,iy,iz)))+
	  eta * MRC_F3(u_p[i],_JZ, ix,iy,iz)
        + d_i * ((MRC_F3(u_p[i], _JX, ix,iy,iz) * MRC_F3(u_p[i], BY, ix,iy,iz) / MRC_F3(u_p[i], RR, ix,iy,iz))- 	
		     (MRC_F3(u_p[i], BX, ix,iy,iz) * MRC_F3(u_p[i], _JY, ix,iy,iz) / MRC_F3(u_p[i], RR, ix,iy,iz)));

	
      // _m 
      MRC_F3(E_m[i], 0, ix,iy,iz) =  
	-((MRC_F3(u_m[i], RVY, ix,iy,iz) * MRC_F3(u_m[i], BZ, ix,iy,iz) / MRC_F3(u_m[i], RR, ix,iy,iz))- 	
	  (MRC_F3(u_m[i], BY, ix,iy,iz) * MRC_F3(u_m[i], RVZ, ix,iy,iz) / MRC_F3(u_m[i], RR, ix,iy,iz)))+
   	  eta * MRC_F3(u_m[i],_JX, ix,iy,iz)
	+ d_i * ((MRC_F3(u_m[i], _JY, ix,iy,iz) * MRC_F3(u_m[i], BZ, ix,iy,iz) / MRC_F3(u_m[i], RR, ix,iy,iz))- 	
		 (MRC_F3(u_m[i], BY, ix,iy,iz) * MRC_F3(u_m[i], _JZ, ix,iy,iz) / MRC_F3(u_m[i], RR, ix,iy,iz)));
      
      MRC_F3(E_m[i], 1, ix,iy,iz) = 
	((MRC_F3(u_m[i], RVX, ix,iy,iz) * MRC_F3(u_m[i], BZ, ix,iy,iz) / MRC_F3(u_m[i], RR, ix,iy,iz))- 	
	 (MRC_F3(u_m[i], BX, ix,iy,iz) * MRC_F3(u_m[i], RVZ, ix,iy,iz) / MRC_F3(u_m[i], RR, ix,iy,iz)))+
	 eta * MRC_F3(u_m[i],_JY, ix,iy,iz)
	- d_i * ((MRC_F3(u_m[i], _JX, ix,iy,iz) * MRC_F3(u_m[i], BZ, ix,iy,iz) / MRC_F3(u_m[i], RR, ix,iy,iz))- 	
		 (MRC_F3(u_m[i], BX, ix,iy,iz) * MRC_F3(u_m[i], _JZ, ix,iy,iz) / MRC_F3(u_m[i], RR, ix,iy,iz)));
      
      MRC_F3(E_m[i], 2, ix,iy,iz) =  
	-((MRC_F3(u_m[i], RVX, ix,iy,iz) * MRC_F3(u_m[i], BY, ix,iy,iz) / MRC_F3(u_m[i], RR, ix,iy,iz))- 	
	  (MRC_F3(u_m[i], BX, ix,iy,iz) * MRC_F3(u_m[i], RVY, ix,iy,iz) / MRC_F3(u_m[i], RR, ix,iy,iz)))+
	  eta * MRC_F3(u_m[i],_JZ, ix,iy,iz) 
        + d_i * ((MRC_F3(u_m[i], _JX, ix,iy,iz) * MRC_F3(u_m[i], BY, ix,iy,iz) / MRC_F3(u_m[i], RR, ix,iy,iz))- 	
		 (MRC_F3(u_m[i], BX, ix,iy,iz) * MRC_F3(u_m[i], _JY, ix,iy,iz) / MRC_F3(u_m[i], RR, ix,iy,iz)));   
    }    
  } mrc_fld_foreach_end;

}
