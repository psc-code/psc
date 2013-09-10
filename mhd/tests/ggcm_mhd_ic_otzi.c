
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_domain.h>
#include <mrc_crds.h>
#include <math.h>

// ======================================================================
// ggcm_mhd_ic subclass "otzi"

// ----------------------------------------------------------------------
// ggcm_mhd_ic_otzi_run

static void
ggcm_mhd_ic_otzi_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_fld *fld = mhd->fld;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);  

  float gamma = mhd->par.gamm;

  mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {
    float r[3];
    r[0] = MRC_CRD(crds, 0, ix);
    r[1] = MRC_CRD(crds, 1, iy);
    r[2] = MRC_CRD(crds, 2, iz);
    MRC_F3(fld, _RR1, ix, iy, iz) = sqr(gamma) ;
    MRC_F3(fld, _RV1X , ix, iy, iz) = sqr(gamma)*sin( r[1] ) ;// - sqr(gamma)*sin( r[1] ) ;
    MRC_F3(fld, _RV1Y , ix, iy, iz) = -sqr(gamma)*sin(r[0]);
    MRC_F3(fld, _RV1Z , ix, iy, iz) = 0.0;
    MRC_F3(fld, _B1X , ix, iy, iz) = sin( 2. *r[1] );  
    MRC_F3(fld, _B1Y , ix, iy, iz) = -sin( r[0] );
    MRC_F3(fld, _B1Z , ix, iy, iz) = 0.0;
    MRC_F3(fld, _UU1 , ix, iy, iz) = (gamma) / (gamma - 1.f) +
      (.5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	      sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	      sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / (sqr(gamma)))+
      (0.5f) * (sqr(MRC_F3(fld, _B1X, ix,iy,iz)) +
		sqr(MRC_F3(fld, _B1Y, ix,iy,iz)) +
		sqr(MRC_F3(fld, _B1Z, ix,iy,iz)));
    
  } mrc_fld_foreach_end;


#if 0
  mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
    float r[3];
    r[0] = MRC_CRD(crds, 0, ix);
    r[1] = MRC_CRD(crds, 1, iy);
    r[2] = MRC_CRD(crds, 2, iz);
    MRC_F3(fld, _RR1 , ix, iy, iz) = 25. / (36.*M_PI );
    MRC_F3(fld, _RV1X, ix, iy, iz) = - sin(2. * M_PI * r[1] ) * MRC_F3(fld, _RR1 , ix, iy, iz);
    MRC_F3(fld, _RV1Y, ix, iy, iz) = sin(2. * M_PI * r[0] )* MRC_F3(fld, _RR1 , ix, iy, iz);
    MRC_F3(fld, _RV1Z, ix, iy, iz) = 0.0;
    MRC_F3(fld, _B1X , ix, iy, iz) = - sqrt(1./(4.*M_PI)) * sin( 2. * M_PI * r[1] ); 
    MRC_F3(fld, _B1Y , ix, iy, iz) =   sqrt(1./(4.*M_PI)) * sin( 4. * M_PI * r[0] );
    MRC_F3(fld, _B1Z , ix, iy, iz) = 0.0;
    MRC_F3(fld, _UU1 , ix, iy, iz) = ((5./(12. * (M_PI))) / (gamma - 1.f)) +
      (.5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	      sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	      sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / MRC_F3(fld, _RR1, ix, iy, iz))+
      (0.5f) * (sqr(MRC_F3(fld, _B1X, ix,iy,iz)) +
		sqr(MRC_F3(fld, _B1Y, ix,iy,iz)) +
		sqr(MRC_F3(fld, _B1Z, ix,iy,iz)));
  } mrc_fld_foreach_end;
#endif
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_otzi_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_otzi_ops = {
  .name        = "otzi",
  .run         = ggcm_mhd_ic_otzi_run,
};
