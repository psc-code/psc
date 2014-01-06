
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_domain.h>
#include <mrc_fld.h>
#include <mrc_fld_as_float.h>

#include <math.h>

// ======================================================================
// ggcm_mhd_ic subclass "sound"

// ----------------------------------------------------------------------
// ggcm_mhd_ic_sound_run

static void
ggcm_mhd_ic_sound_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);  

  float gamma = mhd->par.gamm;
  float cs = sqrt(gamma);
  float pert = 1e-3;
  float xl[3], xh[3];
  mrc_crds_get_param_float3(crds, "l", xl);
  mrc_crds_get_param_float3(crds, "h", xh);
  float k = 2. * M_PI / (xh[0] - xl[0]);

  mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {
    float xx = MRC_CRDX(crds, ix);

    float rr = 1 + pert / cs * sin(k * xx);
    float vx = pert * sin(k * xx);
    float pp = 1 + pert * gamma / cs * sin(k * xx);
    RR1(fld, ix,iy,iz) = rr;
    PP1(fld, ix,iy,iz) = pp;
    V1X(fld, ix,iy,iz) = vx;
  } mrc_fld_foreach_end;

  mrc_fld_put_as(fld, mhd->fld);

  ggcm_mhd_convert_fc_from_primitive(mhd, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_sound_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_sound_ops = {
  .name        = "wave_sound",
  .run         = ggcm_mhd_ic_sound_run,
};


// ======================================================================
// ggcm_mhd class "sound"

// ----------------------------------------------------------------------
// ggcm_mhd_sound_create

static void
ggcm_mhd_sound_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);
}

static struct ggcm_mhd_ops ggcm_mhd_sound_ops = {
  .name             = "sound",
  .create           = ggcm_mhd_sound_create,
};

// ======================================================================
// main

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_sound_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_sound_ops);  
 
  return ggcm_mhd_main(&argc, &argv);
}
