
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_crds.h"

#include <mrc_profile.h>
#include <mrc_io.h>

#include <math.h>
#include <assert.h>

#include <mrc_fld_as_float.h>

#include "mhd_sc.c"

// ----------------------------------------------------------------------
// newstep_c2
//
// like newstep_c, but doesn't need any of the prep work
// (no primvar, primbb, ymask, zmask)

void
newstep_c2(struct ggcm_mhd *mhd, float *dtn)
{
  static int PR;
  if (!PR) {
    PR = prof_register("newstep_c2", 1., 0, 0);
  }
  prof_start(PR);

  struct mrc_fld *f = mrc_fld_get_as(mhd->fld, "float");
  float *fd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FD1);
  float *fd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FD1);
  float *fd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FD1);
  float *fx2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FX2);
  float *fx2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FX2);
  float *fx2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FX2);

  float splim2 = sqr(mhd->par.speedlimit / mhd->par.vvnorm);
  float isphere2 = sqr(mhd->par.isphere);
  float gamm   = mhd->par.gamm;
  float d_i    = mhd->par.d_i;
  float thx    = mhd->par.thx;
  float eps    = 1e-9f;
  float dt     = 1e10f;
  float va02i  = 1.f / sqr(mhd->par.speedlimit / mhd->par.vvnorm);
  float epsz   = 1e-15f;
  float s      = gamm - 1.f;

  float two_pi_d_i = 2. * M_PI * d_i;
  bool have_hall = d_i > 0.f;
  mrc_fld_foreach(f, ix, iy, iz, 0, 0) {
    float hh = fmaxf(fmaxf(fd1x[ix], fd1y[iy]), fd1z[iz]);
    float rri = 1.f / fabsf(MRC_F3(f,_RR1, ix,iy,iz)); // FIXME abs necessary?
    float bb = 
      sqr(.5f*(MRC_F3(f,_B1X, ix,iy,iz)+MRC_F3(f,_B1X, ix-1,iy,iz))) + 
      sqr(.5f*(MRC_F3(f,_B1Y, ix,iy,iz)+MRC_F3(f,_B1Y, ix,iy-1,iz))) +
      sqr(.5f*(MRC_F3(f,_B1Z, ix,iy,iz)+MRC_F3(f,_B1Z, ix,iy,iz-1)));
    if (have_hall) {
      bb *= 1 + sqr(two_pi_d_i * hh);
    }
    float vv1 = fminf(bb * rri, splim2);
    
    float rv2 = 
      sqr(MRC_F3(f,_RV1X, ix,iy,iz)) +
      sqr(MRC_F3(f,_RV1Y, ix,iy,iz)) +
      sqr(MRC_F3(f,_RV1Z, ix,iy,iz));
    float rvv = rri * rv2;
    float pp = s * (MRC_F3(f,_UU1, ix,iy,iz) - .5f * rvv);
    float vv2 = gamm * fmaxf(0.f, pp) * rri;
    float vv3 = rri * sqrtf(rv2);
    float vv = sqrtf(vv1 + vv2) + vv3;
    vv = fmaxf(eps, vv);
    
    float ymask = 1.f;
    if (fx2x[ix] + fx2y[iy] + fx2z[iz] < isphere2)
      ymask = 0.f;
    
    float rrm = fmaxf(epsz, bb * va02i);
    float zmask = ymask * fminf(1.f, MRC_F3(f, _RR1, ix,iy,iz) / rrm);
    
    float tt = thx / fmaxf(eps, hh*vv*zmask);
    dt = fminf(dt, tt);
  } mrc_fld_foreach_end;

  mrc_fld_put_as(f, mhd->fld);

  MPI_Allreduce(&dt, dtn, 1, MPI_FLOAT, MPI_MIN, mhd->obj.comm);

  prof_stop(PR);
}

void
newstep(struct ggcm_mhd *mhd, float *dtn)
{
  ggcm_mhd_fill_ghosts(mhd, mhd->fld, _RR1, mhd->time);

  int mhd_type;
  mrc_fld_get_param_int(mhd->fld, "mhd_type", &mhd_type);
  if (mhd_type == MT_SEMI_CONSERVATIVE_GGCM) {
    newstep_c2(mhd, dtn);
  } else if (mhd_type == MT_SEMI_CONSERVATIVE) {
    zmaskn(mhd, mhd->fld);
    *dtn = newstep_sc(mhd, mhd->fld);
  } else {
    assert(0);
  }
}
