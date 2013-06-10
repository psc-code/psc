
#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_commu.h"
#include "ggcm_mhd_crds.h"

#include "ggcm_mhd_crds.h"

#include <mrc_domain.h>
#include <assert.h>

// ======================================================================
// ggcm_mhd_bnd subclass "conducting_x"

// FIXME: now hardcoded for conducting_x walls at both bounds in y 
// should work ok for Fadeev ic as currently written but would want more flexibility. 

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_conducting_x_fill_ghosts

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
#define OSDx2l(fld, i, ix,iy,iz,s)				     \
  ((-3.0*MRC_F3(fld, i, ix,iy,iz) + 4.0* MRC_F3(fld, i, ix+1,iy,iz)  \
   - MRC_F3(fld, i, ix+2, iy, iz) ) / s ) 
#define OSDx2h(fld, i, ix,iy,iz,s)				     \
  (( 3.0*MRC_F3(fld, i, ix,iy,iz) - 4.0* MRC_F3(fld, i, ix-1,iy,iz)  \
   + MRC_F3(fld, i, ix-2, iy, iz) ) / s ) 

static void
ggcm_mhd_bnd_conducting_x_fill_ghosts(struct ggcm_mhd_bnd *bnd, int m,
				float bntim)
{
  struct ggcm_mhd *mhd = bnd->mhd;

  ggcm_mhd_commu_run(mhd->commu, m, m+8);

  struct mrc_fld *f3 = mrc_fld_get_as(mhd->fld, "float");
  const int *dims = mrc_fld_dims(f3);
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
  float *bd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD1);
  //  float *bd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD1);
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
  if (bc[0] != BC_PERIODIC && info.off[0] == 0) { // x lo    
    // transverse magnetic extrapolated
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int iy = -sw; iy < ny + sw; iy++) {
	// bd1y and bd2y indices offset by one 
	MRC_F3(f3,_B1Y, -1,iy,iz) = MRC_F3(f3, _B1Y, 0,iy,iz) -
	  (1./bd1x[-1])*OSDx2l(f3, _B1Y, 0,iy,iz,2.*bd2x[0] );
        MRC_F3(f3,_B1Z, -1,iy,iz) =  MRC_F3(f3, _B1Z, 0,iy,iz) -
	  (1./bd1x[-1])*OSDx2l(f3, _B1Z, 0,iy,iz,2.*bd2x[0] );
      }
    }

    // set normal magnetic field component for divB=0
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int iy = -sw; iy < ny + sw; iy++) {
	MRC_F3(f3,_B1X, -1,iy,iz) = MRC_F3(f3,_B1X, 0,iy,iz) + bd2x[0] *  
	  ((MRC_F3(f3,_B1Y, 0,iy,iz) - MRC_F3(f3,_B1Y, 0,iy-1,iz) ) / bd2y[iy] + 
	   (MRC_F3(f3,_B1Z, 0,iy,iz) - MRC_F3(f3,_B1Z, 0,iy,iz-1) ) / bd2z[iz]);
      }
    }	
    // transverse magnetic field extrapolated 
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int iy = -sw; iy < ny + sw; iy++) {       
	// bd1x and bd2y indices offset by one 
	MRC_F3(f3,_B1Y, -2,iy,iz) = MRC_F3(f3,_B1Y, -1,iy,iz) -
	  (1./bd1x[-2])*OSDx2l(f3, _B1Y, -1,iy,iz,2.*bd2x[-1] );
	MRC_F3(f3,_B1Z, -2,iy,iz) = MRC_F3(f3,_B1Z, -1,iy,iz) -
	  (1./bd1x[-2])*OSDx2l(f3, _B1Z, -1,iy,iz,2.*bd2x[-1] );
      }
    }
    // set normal magnetic field component for divB=0
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int iy = -sw; iy < ny + sw; iy++) {
	MRC_F3(f3,_B1X, -2,iy,iz) = MRC_F3(f3,_B1X, -1,iy,iz) + bd2x[-1] * 
	  ((MRC_F3(f3,_B1Y, -1,iy,iz) - MRC_F3(f3,_B1Y, -1,iy-1,iz) ) / bd2y[iy] + 
           (MRC_F3(f3,_B1Z, -1,iy,iz) - MRC_F3(f3,_B1Z, -1,iy,iz-1) ) / bd2z[iz]);      
      }
    }	
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int iy = -sw; iy < ny + sw; iy++) {	
	// impenetrable wall 
	MRC_F3(f3,_RV1X, -2,iy,iz) = -MRC_F3(f3,_RV1X, 1,iy,iz);	
	MRC_F3(f3,_RV1X, -1,iy,iz) = -MRC_F3(f3,_RV1X, 0,iy,iz);
	MRC_F3(f3,_RR1, -1,iy,iz) = MRC_F3(f3,_RR1, 0,iy,iz);
	MRC_F3(f3,_RR1, -2,iy,iz) = MRC_F3(f3,_RR1, 1,iy,iz);
	
	// the rest are extrapolations 
	MRC_F3(f3,_RV1Y, -1,iy,iz) = 2.*MRC_F3(f3,_RV1Y, 0,iy,iz)-MRC_F3(f3,_RV1Y, 1,iy,iz);	
	MRC_F3(f3,_RV1Y, -2,iy,iz) = 2.*MRC_F3(f3,_RV1Y, -1,iy,iz)-MRC_F3(f3,_RV1Y, 0,iy,iz);
	
	MRC_F3(f3,_RV1Z, -1,iy,iz) = 2.*MRC_F3(f3,_RV1Z, 0,iy,iz)-MRC_F3(f3,_RV1Z, 1,iy,iz);	
	MRC_F3(f3,_RV1Z, -2,iy,iz) = 2.*MRC_F3(f3,_RV1Z, -1,iy,iz)-MRC_F3(f3,_RV1Z, 0,iy,iz);
	
	MRC_F3(f3,_UU1, -1,iy,iz) = 2.*MRC_F3(f3,_UU1, 0,iy,iz)-MRC_F3(f3,_UU1, 1,iy,iz);	
	MRC_F3(f3,_UU1, -2,iy,iz) = 2.*MRC_F3(f3,_UU1, -1,iy,iz)-MRC_F3(f3,_UU1, 0,iy,iz);	
      }
    }
  }
  //-----------------------------------------------------------------------------------//
  // upper bound 
  //-----------------------------------------------------------------------------------//
  if (bc[0] != BC_PERIODIC && info.off[0] + info.ldims[0] == gdims[0]) { // x hi
    //  transverse magnetic field extrapolated
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int iy = -sw; iy < ny + sw; iy++) {
	MRC_F3(f3,_B1Y, nx,iy,iz) = MRC_F3(f3, _B1Y, nx-1,iy,iz) + 
	  (1./bd1x[nx])*OSDx2h(f3, _B1Y, nx-1,iy,iz,2.*bd2x[nx+1]); 
	MRC_F3(f3,_B1Z, nx,iy,iz) = MRC_F3(f3, _B1Z, nx-1,iy,iz) +
	  (1./bd1x[nx])*OSDx2h(f3, _B1Z, nx-1,iy,iz,2.*bd2x[nx+1]); 
      }
    }
    // set normal magnetic field component for divB=0
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int iy = -sw; iy < ny + sw; iy++) {
	MRC_F3(f3,_B1X, nx-1,iy,iz) = MRC_F3(f3,_B1X, nx-2,iy,iz) - bd2x[nx-1] *
	  ((MRC_F3(f3,_B1Y, nx-1,iy,iz) - MRC_F3(f3,_B1Y, nx-1,iy-1,iz) ) / bd2y[iy] +
	   (MRC_F3(f3,_B1Z, nx-1,iy,iz) - MRC_F3(f3,_B1Z, nx-1,iy,iz-1) ) / bd2z[iz]);
      }
    }	
    //  transverse magnetic field extrapolated
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int iy = -sw; iy < ny + sw; iy++) {
	MRC_F3(f3,_B1Y, nx+1,iy,iz) = MRC_F3(f3, _B1Y, nx,iy,iz) + 
	  (1./bd1x[nx])* OSDx2h(f3, _B1Y, nx+1,iy,iz,2.*bd2x[nx+1] );
        MRC_F3(f3,_B1Z, nx+1,iy,iz) = MRC_F3(f3, _B1Z, nx,iy,iz) + 
	  (1./bd1x[nx])*OSDx2h(f3, _B1Z, nx+1,iy,iz,2.*bd2x[nx+1]);
      }
    }
    // set normal magnetic field component for divB=0
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int iy = -sw; iy < ny + sw; iy++) {
	MRC_F3(f3,_B1X, nx,iy,iz) = MRC_F3(f3,_B1X, nx-1,iy,iz) - bd2x[nx] *  
	  ((MRC_F3(f3,_B1Y, nx,iy,iz) - MRC_F3(f3,_B1Y, nx,iy-1,iz) ) / bd2y[iy] + 
	   (MRC_F3(f3,_B1Z, nx,iy,iz) - MRC_F3(f3,_B1Z, nx,iy,iz-1) ) / bd2z[iz]);
      }
    }	
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int iy = -sw; iy < ny + sw; iy++) {

	// impenetrable wall 	
	MRC_F3(f3,_RV1X, nx+1,iy,iz) = -MRC_F3(f3,_RV1X, nx-2,iy,iz);	
	MRC_F3(f3,_RV1X, nx,iy,iz) = -MRC_F3(f3,_RV1X, nx-1,iy,iz);
	MRC_F3(f3,_RR1, nx+1,iy,iz) = MRC_F3(f3,_RR1, nx-2,iy,iz);	
	MRC_F3(f3,_RR1, nx,iy,iz) = MRC_F3(f3,_RR1, nx-1,iy,iz);

	// the rest are extrapolations 
	MRC_F3(f3,_RV1Y, nx,iy,iz) = MRC_F3(f3,_RV1Y, nx-1,iy,iz) + 
	  (1./bd1x[nx-1]) * OSDx2h(f3, _RV1Y, nx-1,iy,iz,2.*bd2x[nx-1]);  	
	MRC_F3(f3,_RV1Y, nx+1,iy,iz) = MRC_F3(f3,_RV1Y, nx,iy,iz) + 
	  (1./bd1x[nx]) * OSDx2h(f3, _RV1Y, nx,iy,iz,2.*bd2x[nx]);  	
	
	MRC_F3(f3,_RV1Z, nx,iy,iz) = MRC_F3(f3,_RV1Z, nx-1,iy,iz) + 
	  (1./bd1x[nx-1]) * OSDx2h(f3, _RV1Z, nx-1,iy,iz,2.*bd2x[nx-1]);  	
	MRC_F3(f3,_RV1Z, nx+1,iy,iz) = MRC_F3(f3,_RV1Z, nx,iy,iz) +
	  (1./bd1x[nx]) * OSDx2h(f3, _RV1Z, nx,iy,iz,2.*bd2x[nx]);  	

	MRC_F3(f3,_UU1, nx,iy,iz) = MRC_F3(f3,_UU1, nx-1,iy,iz) + 
	  (1./bd1x[nx-1]) * OSDx2h(f3, _UU1, nx-1,iy,iz,2.*bd2x[nx-1]);  	
	MRC_F3(f3,_UU1, nx+1,iy,iz) = MRC_F3(f3,_UU1, nx,iy,iz) +
	  (1./bd1x[nx]) * OSDx2h(f3, _UU1, nx,iy,iz,2.*bd2x[nx]);  	
      }
    }
  }

  mrc_fld_put_as(f3, mhd->fld);
}


// ----------------------------------------------------------------------
// ggcm_mhd_bnd_conducting_x_ops

struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_conducting_x_ops = {
  .name        = "conducting_x",
  .fill_ghosts = ggcm_mhd_bnd_conducting_x_fill_ghosts,
};

