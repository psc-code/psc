
#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_commu.h"
#include "ggcm_mhd_flds.h"

#include <mrc_domain.h>

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

static void
ggcm_mhd_bnd_conducting_fill_ghosts(struct ggcm_mhd_bnd *bnd, int m,
				float bntim)
{
  struct ggcm_mhd *mhd = bnd->mhd;

  ggcm_mhd_commu_run(mhd->commu, m, m+8);

  struct ggcm_mhd_flds *flds = ggcm_mhd_flds_get_as(mhd->flds_base, "c");
  struct mrc_f3 *f3 = ggcm_mhd_flds_get_mrc_f3(flds);
  const int *dims = mrc_f3_dims(f3);
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

  if (bc[1] != BC_PERIODIC && info.off[1] == 0) { // x lo
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {

	MRC_F3(f3,_JX, ix,-2,iz) = 0.0 ;
	MRC_F3(f3,_JX, ix,-1,iz) = 0.0 ;
	MRC_F3(f3,_JX, ix, 0,iz) = 0.0 ;

	MRC_F3(f3,_JZ, ix,-2,iz) = 0.0 ;
	MRC_F3(f3,_JZ, ix,-1,iz) = 0.0 ;
	MRC_F3(f3,_JZ, ix, 0,iz) = 0.0 ;

	MRC_F3(f3,_EX, ix,-2,iz) = 0.0 ;
	MRC_F3(f3,_EX, ix,-1,iz) = 0.0 ;
	
	MRC_F3(f3,_EZ, ix,-2,iz) = 0.0 ;
	MRC_F3(f3,_EZ, ix,-1,iz) = 0.0 ;
	
	B1Y(f3, ix,-1,iz) = 0.0;
	B1Y(f3, ix, 0,iz) = 0.0;

	MRC_F3(f3,_RV1Y, ix,-2,iz) = 0.0;
	MRC_F3(f3,_RV1Y, ix,-1,iz) = 0.0;
	MRC_F3(f3,_RV1Y, ix,0,iz) = 0.0;
	MRC_F3(f3,_RV1Y, ix,1,iz) = 0.0;


      }
    }
  }
    
 if (bc[1] != BC_PERIODIC && info.off[1] + info.ldims[1] == gdims[1]) { // x hi
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int ix = -sw; ix < nx + sw; ix++) {

	MRC_F3(f3,_JX, ix,ny+1,iz) = 0.0 ;
	MRC_F3(f3,_JX, ix,ny,iz) = 0.0 ;
	MRC_F3(f3,_JX, ix,ny-1,iz) = 0.0 ;

	MRC_F3(f3,_JZ, ix,ny+1,iz) = 0.0 ;
	MRC_F3(f3,_JZ, ix,ny,iz) = 0.0 ;
	MRC_F3(f3,_JZ, ix,ny-1,iz) = 0.0 ;

	MRC_F3(f3,_EX, ix,ny+1,iz) = 0.0 ;
	MRC_F3(f3,_EX, ix,ny,iz) = 0.0 ;	
	MRC_F3(f3,_EZ, ix,ny+1,iz) = 0.0 ;
	MRC_F3(f3,_EZ, ix,ny,iz) = 0.0 ;
       
	MRC_F3(f3,_JX, ix,ny+1,iz) = 0.0 ;
	MRC_F3(f3,_JX, ix,ny,iz) = 0.0 ;	
	MRC_F3(f3,_JZ, ix,ny+1,iz) = 0.0 ;
	MRC_F3(f3,_JZ, ix,ny,iz) = 0.0 ;

	B1Y(f3, ix,ny+1,iz) = 0.0; 
	B1Y(f3, ix,ny,iz) = 0.0; 

	MRC_F3(f3,_RV1Y, ix,ny+1,iz) = 0.0;
	MRC_F3(f3,_RV1Y, ix,ny,iz) = 0.0; 
	MRC_F3(f3,_RV1Y, ix,ny-1,iz) = 0.0;
      }
    }
  }

  ggcm_mhd_flds_put_as(flds, mhd->flds_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_conducting_ops

struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_conducting_ops = {
  .name        = "conducting",
  .fill_ghosts = ggcm_mhd_bnd_conducting_fill_ghosts,
};

