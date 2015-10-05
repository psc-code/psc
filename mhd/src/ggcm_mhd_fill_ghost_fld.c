#include "ggcm_mhd_step_cweno_private.h" 

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_fld_as_float.h>
 

// ----------------------------------------------------------------------
// fill_ghost_fld 
// 
// This fills ghost cells for fld objects that have been duplicated in the 
// time-stepping routine (c.f. mrc_ts_rk2.c). Without this, zero-values at 
// boundaries will give inf and nan values for non-periodic boudnaries.

void __unused 
ggcm_mhd_fill_ghost_fld(struct ggcm_mhd *mhd, struct mrc_fld *_fld)
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
	      F3(fld,m, ix,iy,iz) = F3(f3,m, ix,iy,iz); 
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
	      F3(fld,m, ix,iy,iz) = F3(f3,m, ix,iy,iz);
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
	      F3(fld,m, ix,iy,iz) = F3(f3,m, ix,iy,iz); 
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
	      F3(fld,m, ix,iy,iz) = F3(f3,m, ix,iy,iz);
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
	      F3(fld,m, ix,iy,iz) = F3(f3,m, ix,iy,iz); 
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
	      F3(fld,m, ix,iy,iz) = F3(f3,m, ix,iy,iz);
		}		
	  }
	}
       }
    }    
  }  

  mrc_fld_put_as(f3, mhd->fld);
  mrc_fld_put_as(fld, _fld);
}

