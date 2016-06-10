#include "ggcm_mhd_step_cweno_private.h" 

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>

// ----------------------------------------------------------------------
// calc_u_cweno
//
// calculate point values of variables at cell surface centers using 
// thrid order CWENO reconstrcution 
// (Kurganov & Levy 2000 section 2.1)  (Ziegler 2004) 

void
calc_u_cweno(struct ggcm_mhd *mhd, struct mrc_fld *u_p[3], struct mrc_fld *u_m[3],
	  struct mrc_fld *E_p[3], struct mrc_fld *E_m[3],
	  struct mrc_fld *u, struct mrc_fld *u_delta[3])
{
  float d_i = mhd->par.d_i;
  float eta = mhd->par.diffco;
  
  float *bdx1 = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD1);
  float *bdy1 = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD1); 
  float *bdz1 = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD1); 
  
  float *bdx3 = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD3);
  float *bdy3 = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD3); 
  float *bdz3 = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD3); 
  
  float cl = 0.25, cr = 0.25, cc = 0.5; 
  float eps = 1e-18, pp = 2.0;   
  
  float Ir0, Il0, Ic0, ar0, al0, ac0, at0; 
  float Ir1, Il1, Ic1, ar1, al1, ac1, at1; 
  float Ir2, Il2, Ic2, ar2, al2, ac2, at2; 
  
  float duR0, duL0, dm0, dp0; 
  float duR1, duL1, dm1, dp1; 
  float duR2, duL2, dm2, dp2; 
  
  
  mrc_fld_foreach(u, ix,iy,iz, 2, 1) {    
    for (int m = 0; m <= UU; m++){ 

      //reuse u_j+1-u_j/dx  for adjacent Pm and Pp to half the cost of computation. 	
      //Pj(x) = wl Pl + wr Pr + wc  Pc 
      
      // pre-compute some quantities for re-use
      duR0 = MRC_F3(u, m, ix+1,iy,iz) - MRC_F3(u, m,ix,iy,iz);
      duL0 = MRC_F3(u, m, ix,iy,iz) - MRC_F3(u, m,ix-1,iy,iz);
      dm0 = (duR0-duL0); 
      dp0 = (duR0+duL0); 
      
      duR1 = MRC_F3(u, m, ix,iy+1,iz) - MRC_F3(u, m,ix,iy,iz);
      duL1 = MRC_F3(u, m,ix,iy,iz) - MRC_F3(u, m,ix,iy-1,iz);
      dm1 = (duR1-duL1); 
      dp1 = (duR1+duL1);
      
      duR2 = MRC_F3(u, m, ix,iy,iz+1) - MRC_F3(u, m,ix,iy,iz);
      duL2 = MRC_F3(u, m,ix,iy,iz) - MRC_F3(u, m,ix,iy,iz-1);
      dm2 = (duR2-duL2); 
      dp2 = (duR2+duL2); 
      
      // compute smoothness indicators and weightes using only density. 
      // computing smoothness indicators for all quantities increasaes numerical viscosity (c.f. 
      if ( m == 0 ) { 
	Ir0 = duR0*duR0;
	Il0 = duL0*duL0;
	Ic0 = (13./3.)*dm0*dm0 + (0.25)*dp0*dp0;
	ar0 = cr / pow(eps+Ir0,pp); 
	al0 = cl / pow(eps+Il0,pp);
	ac0 = cc / pow(eps+Ic0,pp);
	at0 = ar0+al0+ac0 ; 
	
	Ir1 = duR1*duR1;
	Il1 = duL1*duL1;
	Ic1 = (13./3.)*dm1*dm1 + (0.25)*dp1*dp1 ;
	ar1 = cr / pow(eps+Ir1,pp); 
	al1 = cl / pow(eps+Il1,pp);
	ac1 = cc / pow(eps+Ic1,pp);
	at1 = ar1+al1+ac1 ; 	
	
	Ir2 = duR2*duR2;
	Il2 = duL2*duL2;
	Ic2 = (13./3.)*dm2*dm2 + (0.25)*dp2*dp2 ;
	ar2 = cr / pow(eps+Ir2,pp); 
	al2 = cl / pow(eps+Il2,pp);
	ac2 = cc / pow(eps+Ic2,pp);
	at2 = ar2+al2+ac2 ; 
      }
    
      // x 
      MRC_F3(u_p[0], m, ix,iy,iz) =    
	(ar0 * ( MRC_F3(u, m, ix,iy,iz) + 0.5*duR0 * bdx3[ix] / (bdx1[ix]) ) + 
	 al0 * ( MRC_F3(u, m, ix,iy,iz) + 0.5*duL0 * bdx3[ix] / (bdx1[ix]) ) + 
	 ac0 * ( MRC_F3(u, m, ix,iy,iz) - dm0 / 12. - dm1 / 12. -  dm2 / 12. + 
		dp0 * 0.25 * bdx3[ix] / (bdx1[ix])  +  dm0 * 0.25 * ( pow( bdx3[ix] / bdx1[ix],2.)) ))/at0 ;
      MRC_F3(u_m[0], m, ix,iy,iz) =  
	(ar0 * ( MRC_F3(u, m, ix,iy,iz) + 0.5*duR0 * bdx3[ix] / (-bdx1[ix]) ) + 
	 al0 * ( MRC_F3(u, m, ix,iy,iz) + 0.5*duL0 * bdx3[ix] / (-bdx1[ix]) ) + 
	 ac0 * ( MRC_F3(u, m, ix,iy,iz) - dm0 / 12. - dm1 / 12. - dm2 / 12. + 
		dp0 * 0.25 * bdx3[ix] / (-bdx1[ix])  +  dm0 * 0.25 * ( pow( bdx3[ix] / bdx1[ix],2.)) )) /at0 ; 

       // y 		 
        MRC_F3(u_p[1], m, ix,iy,iz) =  
	  (ar1 * ( MRC_F3(u, m, ix,iy,iz) +  0.5*duR1 * bdy3[iy] / bdy1[iy] ) + 
	   al1 * ( MRC_F3(u, m, ix,iy,iz) +  0.5*duL1 * bdy3[iy] / bdy1[iy] ) + 
	   ac1 * ( MRC_F3(u, m, ix,iy,iz) - dm1 / 12. - dm0 / 12. -  dm2 / 12.  + 
		  dp1 * 0.25 * bdy3[iy] / bdy1[iy]  +  dm1 * 0.25  * ( pow( bdy3[iy] / bdy1[iy],2.)) ))/at1 ; 
        MRC_F3(u_m[1], m, ix,iy,iz) =    
	   (ar1 * ( MRC_F3(u, m, ix,iy,iz) +  0.5*duR1 * bdy3[iy] / (-bdy1[iy]) ) + 
	    al1 * ( MRC_F3(u, m, ix,iy,iz) +  0.5*duL1 * bdy3[iy] / (-bdy1[iy]) ) + 
	    ac1 * ( MRC_F3(u, m, ix,iy,iz) - dm2 / 12. -  dm0 / 12. - dm1 / 12.  + 
		   dp1 * 0.25 * bdy3[iy] / (-bdy1[iy])  +  dm1 * 0.25 * ( pow( bdy3[iy] / bdy1[iy],2.)) ))/at1 ;	

       // z		
        MRC_F3(u_p[2], m, ix,iy,iz) =  
	  (ar2 * ( MRC_F3(u, m, ix,iy,iz) + 0.5* duR2 * bdz3[iz] / bdz1[iz] ) + 
	   al2 * ( MRC_F3(u, m, ix,iy,iz) + 0.5*duL2 * bdz3[iz] / bdz1[iz] ) + 
	   ac2 * ( MRC_F3(u, m, ix,iy,iz) - dm2 / 12. - dm0 / 12. - dm1 / 12. + 
		  dp2 * 0.25 * bdz3[iz] / bdz1[iz]  +  dm2 * 0.25 * ( pow( bdz3[iz] / bdz1[iz],2.)) ))/at2 ;
        MRC_F3(u_m[2], m, ix,iy,iz) =    
	  (ar2 * ( MRC_F3(u, m, ix,iy,iz) + 0.5* duR2 * bdz3[iz] / (-bdz1[iz]) ) + 
	   al2 * ( MRC_F3(u, m, ix,iy,iz) + 0.5* duL2 * bdz3[iz] / (-bdz1[iz]) ) + 
	   ac2 * ( MRC_F3(u, m, ix,iy,iz) - dm2 / 12. - dm0 / 12. - dm1 / 12.  + 
		  dp2 * 0.25 * bdz3[iz] / (-bdz1[iz])  +  dm2 * 0.25 * ( pow( bdz3[iz] / bdz1[iz],2.)) ))/at2 ;	
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
  /*
  float *bdx3 = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD3);
  float *bdy3 = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD3);
  float *bdz3 = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD3);
  */
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
