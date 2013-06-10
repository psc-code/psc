
#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_commu.h"
#include "ggcm_mhd_crds.h"

#include <mrc_domain.h>
#include <assert.h>

// ======================================================================
// ggcm_mhd_bnd subclass "conducting"

// FIXME: now hardcoded for conducting walls at both bounds in y 
// should work ok for Fadeev ic as currently written but would want more flexibility. 

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_conducting_fill_ghosts

enum {
  _EX = _B1Z + 1,
  _EY,
  _EZ,
  _JX,
  _JY,
  _JZ,
  __NR_FLDS,
};

// compute one sided finite diff approx. 
#define OSDy2l(fld, i, ix,iy,iz,s)				     \
  ((-3.0*MRC_F3(fld, i, ix,iy,iz) + 4.0* MRC_F3(fld, i, ix,iy+1,iz)  \
   - MRC_F3(fld, i, ix, iy+2, iz) ) / s ) 
#define OSDy2h(fld, i, ix,iy,iz,s)				     \
  (( 3.0*MRC_F3(fld, i, ix,iy,iz) - 4.0* MRC_F3(fld, i, ix,iy-1,iz)  \
   + MRC_F3(fld, i, ix, iy-2, iz) ) / s ) 

static void
ggcm_mhd_bnd_conducting_fill_ghosts(struct ggcm_mhd_bnd *bnd, int m,
				float bntim)
{
  struct ggcm_mhd *mhd = bnd->mhd;

  ggcm_mhd_commu_run(mhd->commu, m, m+8);

  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, "float");
  const int *dims = mrc_fld_dims(fld);
  int nx = dims[0], ny = dims[1], nz = dims[2];
  int sw = SW_2; 
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, 0, &info);
  int bc[3];
  mrc_domain_get_param_int(mhd->domain, "bcx", &bc[0]); // FIXME in libmrc
  mrc_domain_get_param_int(mhd->domain, "bcy", &bc[1]);
  mrc_domain_get_param_int(mhd->domain, "bcz", &bc[2]);


  //struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  //    bdy1[i] = 1.f / (fyy1[i+1] - fyy1[i]);
  //  float *bd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD1);
  float *bd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD1);
  //  float *bd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD1);

  //   bdx2[i] = .5f * (fxx1[i+1] - fxx1[i-1]);
  float *bd2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD2);
  float *bd2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD2);
  float *bd2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD2);

  // fill boundary values with closes non-zero value 
  // FIXME: should this be necessary? 
  
  bd2x[-2] = bd2x[0]; 
  bd2x[-1] = bd2x[0]; 
  bd2x[nx] = bd2x[nx-1];
  bd2x[nx+1] = bd2x[nx-1]; 

  bd2y[-2] = bd2y[0]; 
  bd2y[-1] = bd2y[0]; 
  bd2y[ny] = bd2y[ny-1];
  bd2y[ny+1] = bd2y[ny-1]; 

  bd2z[-2] = bd2z[0]; 
  bd2z[-1] = bd2z[0]; 
  bd2z[nz] = bd2z[nz-1];
  bd2z[nz+1] = bd2z[nz-1]; 
  
//-----------------------------------------------------------------------------------//
// lower bound 
//-----------------------------------------------------------------------------------//
  if (bc[1] != BC_PERIODIC && info.off[1] == 0) { // x lo    
    // transverse magnetic extrapolated
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {
	// bd1y and bd2y indices offset by one 
	MRC_F3(fld,_B1X, ix,-1,iz) = MRC_F3(fld, _B1X, ix,0,iz) -
	  (1./bd1y[-1])*OSDy2l(fld, _B1X, ix,0,iz,2.*bd2y[0] );
        MRC_F3(fld,_B1Z, ix,-1,iz) =  MRC_F3(fld, _B1Z, ix,0,iz) -
	  (1./bd1y[-1])*OSDy2l(fld, _B1Z, ix,0,iz,2.*bd2y[0] );
      }
    }
    // set normal magnetic field component for divB=0
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {
	MRC_F3(fld,_B1Y, ix,-1,iz) = MRC_F3(fld,_B1Y, ix,0,iz) + bd2y[0] *  
	  ((MRC_F3(fld,_B1X, ix,0,iz) - MRC_F3(fld,_B1X, ix-1,0,iz) ) / bd2x[ix] + 
	   (MRC_F3(fld,_B1Z, ix,0,iz) - MRC_F3(fld,_B1Z, ix,0,iz-1) ) / bd2z[iz]);
      }
    }	
    // transverse magnetic field extrapolated 
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {       
	// bd1y and bd2y indices offset by one 
	MRC_F3(fld,_B1X, ix,-2,iz) = MRC_F3(fld,_B1X, ix,-1,iz) -
	  (1./bd1y[-2])*OSDy2l(fld, _B1X, ix,-1,iz,2.*bd2y[-1] );
	MRC_F3(fld,_B1Z, ix,-2,iz) = MRC_F3(fld,_B1Z, ix,-1,iz) -
	  (1./bd1y[-2])*OSDy2l(fld, _B1Z, ix,-1,iz,2.*bd2y[-1] );
      }
    }
    // set normal magnetic field component for divB=0
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {
	MRC_F3(fld,_B1Y, ix,-2,iz) = MRC_F3(fld,_B1Y, ix,-1,iz) + bd2y[-1] * 
	  ((MRC_F3(fld,_B1X, ix,-1,iz) - MRC_F3(fld,_B1X, ix-1,-1,iz) ) / bd2x[ix] + 
           (MRC_F3(fld,_B1Z, ix,-1,iz) - MRC_F3(fld,_B1Z, ix,-1,iz-1) ) / bd2z[iz]);      
      }
    }	
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {	
	// impenetrable wall 
	MRC_F3(fld,_RV1Y, ix,-2,iz) = -MRC_F3(fld,_RV1Y, ix,1,iz);	
	MRC_F3(fld,_RV1Y, ix,-1,iz) = -MRC_F3(fld,_RV1Y, ix,0,iz);
	MRC_F3(fld,_RR1, ix,-1,iz) = MRC_F3(fld,_RR1, ix,0,iz);
	MRC_F3(fld,_RR1, ix,-2,iz) = MRC_F3(fld,_RR1, ix,1,iz);
	
	// the rest are extrapolations 
	MRC_F3(fld,_RV1X, ix,-1,iz) = 2.*MRC_F3(fld,_RV1X, ix,0,iz)-MRC_F3(fld,_RV1X, ix,1,iz);	
	MRC_F3(fld,_RV1X, ix,-2,iz) = 2.*MRC_F3(fld,_RV1X, ix,-1,iz)-MRC_F3(fld,_RV1X, ix,0,iz);
	
	MRC_F3(fld,_RV1Z, ix,-1,iz) = 2.*MRC_F3(fld,_RV1Z, ix,0,iz)-MRC_F3(fld,_RV1Z, ix,1,iz);	
	MRC_F3(fld,_RV1Z, ix,-2,iz) = 2.*MRC_F3(fld,_RV1Z, ix,-1,iz)-MRC_F3(fld,_RV1Z, ix,0,iz);
	
	MRC_F3(fld,_UU1, ix,-1,iz) = 2.*MRC_F3(fld,_UU1, ix,0,iz)-MRC_F3(fld,_UU1, ix,1,iz);	
	MRC_F3(fld,_UU1, ix,-2,iz) = 2.*MRC_F3(fld,_UU1, ix,-1,iz)-MRC_F3(fld,_UU1, ix,0,iz);	
      }
    }
  }
  //-----------------------------------------------------------------------------------//
  // upper bound 
  //-----------------------------------------------------------------------------------//
  if (bc[1] != BC_PERIODIC && info.off[1] + info.ldims[1] == gdims[1]) { // x hi
    //  transverse magnetic field extrapolated
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {
	MRC_F3(fld,_B1X, ix,ny,iz) = MRC_F3(fld, _B1X, ix,ny-1,iz) + 
	  (1./bd1y[ny])*OSDy2h(fld, _B1X, ix,ny-1,iz,2.*bd2y[ny+1]); 
	MRC_F3(fld,_B1Z, ix,ny,iz) = MRC_F3(fld, _B1Z, ix,ny-1,iz) +
	  (1./bd1y[ny])*OSDy2h(fld, _B1Z, ix,ny-1,iz,2.*bd2y[ny+1]); 
      }
    }
    // set normal magnetic field component for divB=0
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {
	MRC_F3(fld,_B1Y, ix,ny-1,iz) = MRC_F3(fld,_B1Y, ix,ny-2,iz) - bd2y[ny-1] *  
	  ((MRC_F3(fld,_B1X, ix,ny-1,iz) - MRC_F3(fld,_B1X, ix-1,ny-1,iz) ) / bd2x[ix] +  										
	   (MRC_F3(fld,_B1Z, ix,ny-1,iz) - MRC_F3(fld,_B1Z, ix,ny-1,iz-1) ) / bd2z[iz]);
      }
    }	
    //  transverse magnetic field extrapolated
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {
	MRC_F3(fld,_B1X, ix,ny+1,iz) = MRC_F3(fld, _B1X, ix,ny,iz) + 
	  (1./bd1y[ny])* OSDy2h(fld, _B1X, ix,ny+1,iz,2.*bd2y[ny+1] );
        MRC_F3(fld,_B1Z, ix,ny+1,iz) = MRC_F3(fld, _B1Z, ix,ny,iz) + 
	  (1./bd1y[ny])*OSDy2h(fld, _B1Z, ix,ny+1,iz,2.*bd2y[ny+1]);
      }
    }
    // set normal magnetic field component for divB=0
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {
	MRC_F3(fld,_B1Y, ix,ny,iz) = MRC_F3(fld,_B1Y, ix,ny-1,iz) - bd2y[ny] *  
	  ((MRC_F3(fld,_B1X, ix,ny,iz) - MRC_F3(fld,_B1X, ix-1,ny,iz) ) / bd2x[ix] + 
	   (MRC_F3(fld,_B1Z, ix,ny,iz) - MRC_F3(fld,_B1Z, ix,ny,iz-1) ) / bd2z[iz]);
      }
    }	
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {

	// impenetrable wall 	
	MRC_F3(fld,_RV1Y, ix,ny+1,iz) = -MRC_F3(fld,_RV1Y, ix,ny-2,iz);	
	MRC_F3(fld,_RV1Y, ix,ny,iz) = -MRC_F3(fld,_RV1Y, ix,ny-1,iz);
	MRC_F3(fld,_RR1, ix,ny+1,iz) = MRC_F3(fld,_RR1, ix,ny-2,iz);	
	MRC_F3(fld,_RR1, ix,ny,iz) = MRC_F3(fld,_RR1, ix,ny-1,iz);

	// the rest are extrapolations 
	MRC_F3(fld,_RV1X, ix,ny,iz) = MRC_F3(fld,_RV1X, ix,ny-1,iz) + 
	  (1./bd1y[ny-1]) * OSDy2h(fld, _RV1X, ix,ny-1,iz,2.*bd2y[ny-1]);  	
	MRC_F3(fld,_RV1X, ix,ny+1,iz) = MRC_F3(fld,_RV1X, ix,ny,iz) + 
	  (1./bd1y[ny]) * OSDy2h(fld, _RV1X, ix,ny,iz,2.*bd2y[ny]);  	
	
	MRC_F3(fld,_RV1Z, ix,ny,iz) = MRC_F3(fld,_RV1Z, ix,ny-1,iz) + 
	  (1./bd1y[ny-1]) * OSDy2h(fld, _RV1Z, ix,ny-1,iz,2.*bd2y[ny-1]);  	
	MRC_F3(fld,_RV1Z, ix,ny+1,iz) = MRC_F3(fld,_RV1Z, ix,ny,iz) +
	  (1./bd1y[ny]) * OSDy2h(fld, _RV1Z, ix,ny,iz,2.*bd2y[ny]);  	

	MRC_F3(fld,_UU1, ix,ny,iz) = MRC_F3(fld,_UU1, ix,ny-1,iz) + 
	  (1./bd1y[ny-1]) * OSDy2h(fld, _UU1, ix,ny-1,iz,2.*bd2y[ny-1]);  	
	MRC_F3(fld,_UU1, ix,ny+1,iz) = MRC_F3(fld,_UU1, ix,ny,iz) +
	  (1./bd1y[ny]) * OSDy2h(fld, _UU1, ix,ny,iz,2.*bd2y[ny]);  	
      }
    }
  }

  mrc_fld_put_as(fld, mhd->fld);
}


// ----------------------------------------------------------------------
// ggcm_mhd_bnd_conducting_ops

struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_conducting_ops = {
  .name        = "conducting",
  .fill_ghosts = ggcm_mhd_bnd_conducting_fill_ghosts,
};

