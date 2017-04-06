
#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"

#include <mrc_domain.h>
#include <mrc_fld_as_float.h>

#include <assert.h>

// ======================================================================
// ggcm_mhd_bnd subclass "open_x"

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_open_x_fill_ghosts

enum {
  _EX = BZ + 1,
  _EY,
  _EZ,
  _JX,
  _JY,
  _JZ,
  __NR_FLDS,
};

static void
ggcm_mhd_bnd_open_x_fill_ghosts_scons(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld,
				      float bntim)
{
  struct ggcm_mhd *mhd = bnd->mhd;

  struct mrc_fld *f3 = mrc_fld_get_as(fld, FLD_TYPE);

  assert(mrc_fld_nr_patches(f3) == 1);
  int p = 0;

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int dy = (gdims[1] > 1), dz = (gdims[2] > 1);
  int sw_y = gdims[1] > 1 ? SW_2 : 0; 
  int sw_z = gdims[2] > 1 ? SW_2 : 0; 

  const int *dims = mrc_fld_dims(f3);
  int nx = dims[0], ny = dims[1], nz = dims[2];
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, 0, &info);

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
  if (info.off[0] == 0) { // x lo    
    // transverse magnetic extrapolated
    for (int iz = -sw_z; iz < nz + sw_z; iz++) {
      for (int iy = -sw_y; iy < ny + sw_y; iy++) {
	// bd1y and bd2y indices offset by one 
	M3(f3,BY, -1,iy,iz, p) = M3(f3,BY, 0,iy,iz, p); 
        M3(f3,BZ, -1,iy,iz, p) = M3(f3,BZ, 0,iy,iz, p);
      }
    }

    // set normal magnetic field component for divB=0
    for (int iz = -sw_z; iz < nz + sw_z; iz++) {
      for (int iy = -sw_y; iy < ny + sw_y; iy++) {
	M3(f3,BX, -1,iy,iz, p) = M3(f3,BX, 0,iy,iz, p) + bd2x[0] *  
	  ((M3(f3,BY, 0,iy,iz, p) - M3(f3,BY, 0,iy-dy,iz, p) ) / bd2y[iy] + 
	   (M3(f3,BZ, 0,iy,iz, p) - M3(f3,BZ, 0,iy,iz-dz, p) ) / bd2z[iz]);
      }
    }	
    // transverse magnetic field extrapolated 
    for (int iz = -sw_z; iz < nz + sw_z; iz++) {
      for (int iy = -sw_y; iy < ny + sw_y; iy++) {       
	// bd1x and bd2y indices offset by one 
	M3(f3,BY, -2,iy,iz, p) = M3(f3,BY, 1,iy,iz, p);
	M3(f3,BZ, -2,iy,iz, p) = M3(f3,BZ, 1,iy,iz, p);
      }
    }
    // set normal magnetic field component for divB=0
    for (int iz = -sw_z; iz < nz + sw_z; iz++) {
      for (int iy = -sw_y; iy < ny + sw_y; iy++) {
	M3(f3,BX, -2,iy,iz, p) = M3(f3,BX, -1,iy,iz, p) + bd2x[-1] * 
	  ((M3(f3,BY, -1,iy,iz, p) - M3(f3,BY, -1,iy-dy,iz, p) ) / bd2y[iy] + 
           (M3(f3,BZ, -1,iy,iz, p) - M3(f3,BZ, -1,iy,iz-dz, p) ) / bd2z[iz]);      
      }
    }	
    for (int iz = -sw_z; iz < nz + sw_z; iz++) {
      for (int iy = -sw_y; iy < ny + sw_y; iy++) {	

	M3(f3,RVX, -1,iy,iz, p) = M3(f3,RVX, 0,iy,iz, p);
	M3(f3,RVX, -2,iy,iz, p) = M3(f3,RVX, 1,iy,iz, p);	

	M3(f3,RR, -1,iy,iz, p) = M3(f3,RR, 0,iy,iz, p);
	M3(f3,RR, -2,iy,iz, p) = M3(f3,RR, 1,iy,iz, p);

	M3(f3,RVY, -1,iy,iz, p) = M3(f3,RVY, 0,iy,iz, p); 
	M3(f3,RVY, -2,iy,iz, p) = M3(f3,RVY, 1,iy,iz, p);
	
	M3(f3,RVZ, -1,iy,iz, p) = M3(f3,RVZ, 0,iy,iz, p);
	M3(f3,RVZ, -2,iy,iz, p) = M3(f3,RVZ, 1,iy,iz, p);
	
	M3(f3,UU, -1,iy,iz, p) = M3(f3,UU, 0,iy,iz, p);
	M3(f3,UU, -2,iy,iz, p) = M3(f3,UU, 1,iy,iz, p);

      }
    }
  }
  //-----------------------------------------------------------------------------------//
  // upper bound 
  //-----------------------------------------------------------------------------------//
  if (info.off[0] + info.ldims[0] == gdims[0]) { // x hi
    //  transverse magnetic field extrapolated
    for (int iz = -sw_z; iz < nz + sw_z; iz++) {
      for (int iy = -sw_y; iy < ny + sw_y; iy++) {
	M3(f3,BY, nx,iy,iz, p) = M3(f3, BY, nx-1,iy,iz, p);
	M3(f3,BZ, nx,iy,iz, p) = M3(f3, BZ, nx-1,iy,iz, p);
      }
    }
    // set normal magnetic field component for divB=0
    for (int iz = -sw_z; iz < nz + sw_z; iz++) {
      for (int iy = -sw_y; iy < ny + sw_y; iy++) {
	M3(f3,BX, nx-1,iy,iz, p) = M3(f3,BX, nx-2,iy,iz, p) - bd2x[nx-1] *
	  ((M3(f3,BY, nx-1,iy,iz, p) - M3(f3,BY, nx-1,iy-dy,iz, p) ) / bd2y[iy] +
	   (M3(f3,BZ, nx-1,iy,iz, p) - M3(f3,BZ, nx-1,iy,iz-dz, p) ) / bd2z[iz]);
      }
    }	
    //  transverse magnetic field extrapolated
    for (int iz = -sw_z; iz < nz + sw_z; iz++) {
      for (int iy = -sw_y; iy < ny + sw_y; iy++) {
	M3(f3,BY, nx+1,iy,iz, p) = M3(f3, BY, nx-2,iy,iz, p);
        M3(f3,BZ, nx+1,iy,iz, p) = M3(f3, BZ, nx-2,iy,iz, p); 
      }
    }
    // set normal magnetic field component for divB=0
    for (int iz = -sw_z; iz < nz + sw_z; iz++) {
      for (int iy = -sw_y; iy < ny + sw_y; iy++) {
	M3(f3,BX, nx,iy,iz, p) = M3(f3,BX, nx-1,iy,iz, p) - bd2x[nx] *  
	  ((M3(f3,BY, nx,iy,iz, p) - M3(f3,BY, nx,iy-dy,iz, p) ) / bd2y[iy] + 
	   (M3(f3,BZ, nx,iy,iz, p) - M3(f3,BZ, nx,iy,iz-dz, p) ) / bd2z[iz]);
      }
    }	
    for (int iz = -sw_z; iz < nz + sw_z; iz++) {
      for (int iy = -sw_y; iy < ny + sw_y; iy++) {

	M3(f3,RVX, nx+1,iy,iz, p) = M3(f3,RVX, nx-2,iy,iz, p);	
	M3(f3,RVX, nx,iy,iz, p) = M3(f3,RVX, nx-1,iy,iz, p);

	M3(f3,RR, nx+1,iy,iz, p) = M3(f3,RR, nx-2,iy,iz, p);	
	M3(f3,RR, nx,iy,iz, p) = M3(f3,RR, nx-1,iy,iz, p);

	M3(f3,RVY, nx+1,iy,iz, p) = M3(f3,RVY, nx-2,iy,iz, p); 
	M3(f3,RVY, nx,iy,iz, p) = M3(f3,RVY, nx-1,iy,iz, p);

	M3(f3,RVZ, nx+1,iy,iz, p) = M3(f3,RVZ, nx-2,iy,iz, p);
	M3(f3,RVZ, nx,iy,iz, p) = M3(f3,RVZ, nx-1,iy,iz, p);

	M3(f3,UU, nx+1,iy,iz, p) = M3(f3,UU, nx-2,iy,iz, p);
	M3(f3,UU, nx,iy,iz, p) = M3(f3,UU, nx-1,iy,iz, p); 

      }
    }
  }

  mrc_fld_put_as(f3, fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_open_x_fill_ghosts_fcons_cc

static void
ggcm_mhd_bnd_open_x_fill_ghosts_fcons_cc(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld,
					 float bntim)
{
  struct ggcm_mhd *mhd = bnd->mhd;

  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int sw_y = gdims[1] > 1 ? SW_2 : 0; 
  int sw_z = gdims[2] > 1 ? SW_2 : 0; 

  assert(mrc_fld_nr_patches(f) == 1);
  int p = 0;
  const int *ldims = mrc_fld_spatial_dims(f);
  int mx = ldims[0];
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, 0, &info);
  int nr_comps = mrc_fld_nr_comps(f);

    if (info.off[0] == 0) { // x lo FIXME, there is (?) a function for this
    for (int k = -sw_z; k < ldims[2] + sw_z; k++) {
      for (int j = -sw_y; j < ldims[1] + sw_y; j++) {	
	for (int m = 0; m < nr_comps; m++) {
	  M3(f, m, -1,j,k, p) = M3(f, m, 0,j,k, p);
	  M3(f, m, -2,j,k, p) = M3(f, m, 1,j,k, p);
	}
      }
    }
  }

  if (info.off[0] + info.ldims[0] == gdims[0]) { // x hi
    for (int k = -sw_z; k < ldims[2] + sw_z; k++) {
      for (int j = -sw_y; j < ldims[1] + sw_y; j++) {	
	for (int m = 0; m < nr_comps; m++) {
	  M3(f, m, mx  ,j,k, p) = M3(f, m, mx-1,j,k, p);
	  M3(f, m, mx+1,j,k, p) = M3(f, m, mx-2,j,k, p);
	}
      }
    }
  }

  mrc_fld_put_as(f, fld);
}

static void
ggcm_mhd_bnd_open_x_fill_ghosts(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld,
				float bntim)
{
  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);

  int bcx;
  mrc_domain_get_param_int(bnd->mhd->domain, "bcx", &bcx);
  assert(bcx == BC_NONE);

  if (mhd_type == MT_SCONS_FC) {
    ggcm_mhd_bnd_open_x_fill_ghosts_scons(bnd, fld, bntim);
  } else if (mhd_type == MT_FCONS_CC) {
    ggcm_mhd_bnd_open_x_fill_ghosts_fcons_cc(bnd, fld, bntim);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_open_x_ops

struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_open_x_ops = {
  .name        = "open_x",
  .fill_ghosts = ggcm_mhd_bnd_open_x_fill_ghosts,
};

