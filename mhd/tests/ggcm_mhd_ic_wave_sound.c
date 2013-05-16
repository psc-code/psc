
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_flds.h"

#include <mrc_domain.h>
#include <mrc_fld.h>
#include <math.h>

// ======================================================================
// ggcm_mhd_ic subclass "wave_sound"

// ----------------------------------------------------------------------
// ggcm_mhd_ic_wave_sound_run

static void
ggcm_mhd_ic_wave_sound_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd *gmhd = ic->mhd;
  struct mrc_f3 *f3 = ggcm_mhd_flds_get_mrc_f3(gmhd->flds_base);
  struct mrc_crds *crds = mrc_domain_get_crds(gmhd->domain);  

  float gamma = gmhd->par.gamm;
  float cs = sqrt(gamma);
  float pert = 1e-3;
  float xl[3], xh[3];
  mrc_crds_get_xl_xh(crds, xl, xh);
  float k = 2. * M_PI / (xh[0] - xl[0]);

  mrc_f3_foreach(f3, ix, iy, iz, 1, 1) {
    float xx = MRC_CRDX(crds, ix);

    float rr = 1 + pert / cs * sin(k * xx);
    float vx = pert * sin(k * xx);
    float pp = 1 + pert * gamma / cs * sin(k * xx);
    MRC_F3(f3, _RR1, ix, iy, iz) = rr;
    MRC_F3(f3, _RV1X , ix, iy, iz) = rr * vx;
    MRC_F3(f3, _UU1 , ix, iy, iz) = pp / (gamma - 1);
  } mrc_f3_foreach_end;
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_wave_sound_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_wave_sound_ops = {
  .name        = "wave_sound",
  .run         = ggcm_mhd_ic_wave_sound_run,
};

