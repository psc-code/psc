
#include "ggcm_mhd_step_private.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_crds.h"

#include <mrc_domain.h>

#include <math.h>

// TODO:
// - handle various resistivity models
// - handle limit2, limit3
// - handle lowmask

#define REPS (1.e-10f)

enum {
  LIMIT_NONE,
  LIMIT_1,
};

static void
rmaskn_c(struct ggcm_mhd *mhd)
{
  struct mrc_fld *f = mhd->fld;

  float diffco = mhd->par.diffco;
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, 0, &info);

  float *fx1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FX1);

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    MRC_F3(f,_RMASK, ix,iy,iz) = 0.f;
    float xxx = fx1x[ix];
    if (xxx < -15.f)
      continue;
    if (iy + info.off[1] < 2)
      continue;
    if (iz + info.off[2] < 2)
      continue;
    if (ix + info.off[0] >= gdims[0] - 2)
      continue;
    if (iy + info.off[1] >= gdims[1] - 2)
      continue;
    if (iz + info.off[2] >= gdims[2] - 2)
      continue;

    MRC_F3(f, _RMASK, ix,iy,iz) = diffco * MRC_F3(f, _ZMASK, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static void
vgflrr_c(struct ggcm_mhd *mhd)
{
  struct mrc_fld *f = mhd->fld;

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    float a = MRC_F3(f,_RR, ix,iy,iz);
    MRC_F3(f,_TMP1, ix,iy,iz) = a * MRC_F3(f,_VX, ix,iy,iz);
    MRC_F3(f,_TMP2, ix,iy,iz) = a * MRC_F3(f,_VY, ix,iy,iz);
    MRC_F3(f,_TMP3, ix,iy,iz) = a * MRC_F3(f,_VZ, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static void
vgflrvx_c(struct ggcm_mhd *mhd)
{
  struct mrc_fld *f = mhd->fld;

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    float a = MRC_F3(f,_RR, ix,iy,iz) * MRC_F3(f,_VX, ix,iy,iz);
    MRC_F3(f,_TMP1, ix,iy,iz) = a * MRC_F3(f,_VX, ix,iy,iz);
    MRC_F3(f,_TMP2, ix,iy,iz) = a * MRC_F3(f,_VY, ix,iy,iz);
    MRC_F3(f,_TMP3, ix,iy,iz) = a * MRC_F3(f,_VZ, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static void
vgflrvy_c(struct ggcm_mhd *mhd)
{
  struct mrc_fld *f = mhd->fld;

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    float a = MRC_F3(f,_RR, ix,iy,iz) * MRC_F3(f,_VY, ix,iy,iz);
    MRC_F3(f,_TMP1, ix,iy,iz) = a * MRC_F3(f,_VX, ix,iy,iz);
    MRC_F3(f,_TMP2, ix,iy,iz) = a * MRC_F3(f,_VY, ix,iy,iz);
    MRC_F3(f,_TMP3, ix,iy,iz) = a * MRC_F3(f,_VZ, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static void
vgflrvz_c(struct ggcm_mhd *mhd)
{
  struct mrc_fld *f = mhd->fld;

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    float a = MRC_F3(f,_RR, ix,iy,iz) * MRC_F3(f,_VZ, ix,iy,iz);
    MRC_F3(f,_TMP1, ix,iy,iz) = a * MRC_F3(f,_VX, ix,iy,iz);
    MRC_F3(f,_TMP2, ix,iy,iz) = a * MRC_F3(f,_VY, ix,iy,iz);
    MRC_F3(f,_TMP3, ix,iy,iz) = a * MRC_F3(f,_VZ, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static void
vgfluu_c(struct ggcm_mhd *mhd)
{
  struct mrc_fld *f = mhd->fld;

  float gamma = mhd->par.gamm;
  float s = gamma / (gamma - 1.f);
  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    float ep = s * MRC_F3(f,_PP, ix,iy,iz) +
      .5f * MRC_F3(f,_RR, ix,iy,iz) * (sqr(MRC_F3(f,_VX, ix,iy,iz)) + 
				       sqr(MRC_F3(f,_VY, ix,iy,iz)) + 
				       sqr(MRC_F3(f,_VZ, ix,iy,iz)));
    MRC_F3(f,_TMP1, ix,iy,iz) = ep * MRC_F3(f,_VX, ix,iy,iz);
    MRC_F3(f,_TMP2, ix,iy,iz) = ep * MRC_F3(f,_VY, ix,iy,iz);
    MRC_F3(f,_TMP3, ix,iy,iz) = ep * MRC_F3(f,_VZ, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static void
fluxl_c(struct ggcm_mhd *mhd, int m)
{
  struct mrc_fld *f = mhd->fld;

  mrc_fld_foreach(f, ix,iy,iz, 1, 0) {
    float aa = MRC_F3(f,m, ix,iy,iz);
    float cmsv = MRC_F3(f,_CMSV, ix,iy,iz);
    MRC_F3(f,_FLX, ix,iy,iz) =
      .5f * ((MRC_F3(f,_TMP1, ix  ,iy,iz) + MRC_F3(f,_TMP1, ix+1,iy,iz)) -
	     .5f * (MRC_F3(f,_CMSV, ix+1,iy,iz) + cmsv) * (MRC_F3(f,m, ix+1,iy,iz) - aa));
    MRC_F3(f,_FLY, ix,iy,iz) =
      .5f * ((MRC_F3(f,_TMP2, ix,iy  ,iz) + MRC_F3(f,_TMP2, ix,iy+1,iz)) -
	     .5f * (MRC_F3(f,_CMSV, ix,iy+1,iz) + cmsv) * (MRC_F3(f,m, ix,iy+1,iz) - aa));
    MRC_F3(f,_FLZ, ix,iy,iz) =
      .5f * ((MRC_F3(f,_TMP3, ix,iy,iz  ) + MRC_F3(f,_TMP3, ix,iy,iz+1)) -
	     .5f * (MRC_F3(f,_CMSV, ix,iy,iz+1) + cmsv) * (MRC_F3(f,m, ix,iy,iz+1) - aa));
  } mrc_fld_foreach_end;
}

static void
fluxb_c(struct ggcm_mhd *mhd, int m)
{
  struct mrc_fld *f = mhd->fld;

  float s1 = 1.f/12.f;
  float s7 = 7.f * s1;

  mrc_fld_foreach(f, ix,iy,iz, 1, 0) {
    float fhx = (s7 * (MRC_F3(f, _TMP1, ix  ,iy,iz) + MRC_F3(f, _TMP1, ix+1,iy,iz)) -
		 s1 * (MRC_F3(f, _TMP1, ix-1,iy,iz) + MRC_F3(f, _TMP1, ix+2,iy,iz)));
    float fhy = (s7 * (MRC_F3(f, _TMP2, ix,iy  ,iz) + MRC_F3(f, _TMP2, ix,iy+1,iz)) -
		 s1 * (MRC_F3(f, _TMP2, ix,iy-1,iz) + MRC_F3(f, _TMP2, ix,iy+2,iz)));
    float fhz = (s7 * (MRC_F3(f, _TMP3, ix,iy,iz  ) + MRC_F3(f, _TMP3, ix,iy,iz+1)) -
		 s1 * (MRC_F3(f, _TMP3, ix,iy,iz-1) + MRC_F3(f, _TMP3, ix,iy,iz+2)));

    float aa = MRC_F3(f,m, ix,iy,iz);
    float cmsv = MRC_F3(f,_CMSV, ix,iy,iz);
    float flx =
      .5f * ((MRC_F3(f,_TMP1, ix  ,iy,iz) + MRC_F3(f,_TMP1, ix+1,iy,iz)) -
	     .5f * (MRC_F3(f,_CMSV, ix+1,iy,iz) + cmsv) * (MRC_F3(f,m, ix+1,iy,iz) - aa));
    float fly =
      .5f * ((MRC_F3(f,_TMP2, ix,iy  ,iz) + MRC_F3(f,_TMP2, ix,iy+1,iz)) -
	     .5f * (MRC_F3(f,_CMSV, ix,iy+1,iz) + cmsv) * (MRC_F3(f,m, ix,iy+1,iz) - aa));
    float flz = 
      .5f * ((MRC_F3(f,_TMP3, ix,iy,iz  ) + MRC_F3(f,_TMP3, ix,iy,iz+1)) -
	     .5f * (MRC_F3(f,_CMSV, ix,iy,iz+1) + cmsv) * (MRC_F3(f,m, ix,iy,iz+1) - aa));

    float cx = MRC_F3(f, _CX, ix,iy,iz);
    MRC_F3(f, _FLX, ix,iy,iz) = cx * flx + (1.f - cx) * fhx;
    float cy = MRC_F3(f, _CY, ix,iy,iz);
    MRC_F3(f, _FLY, ix,iy,iz) = cy * fly + (1.f - cy) * fhy;
    float cz = MRC_F3(f, _CZ, ix,iy,iz);
    MRC_F3(f, _FLZ, ix,iy,iz) = cz * flz + (1.f - cz) * fhz;
  } mrc_fld_foreach_end;
}

static void
pushn_c(struct ggcm_mhd *mhd, int ma, int mc, float dt)
{
  struct mrc_fld *f = mhd->fld;
  float *fd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FD1);
  float *fd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FD1);
  float *fd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FD1);

  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    float s = dt * MRC_F3(f,_YMASK, ix,iy,iz);
    MRC_F3(f,mc, ix,iy,iz) = MRC_F3(f,ma, ix,iy,iz)
      - s * (fd1x[ix] * (MRC_F3(f,_FLX, ix,iy,iz) - MRC_F3(f,_FLX, ix-1,iy,iz)) +
	     fd1y[iy] * (MRC_F3(f,_FLY, ix,iy,iz) - MRC_F3(f,_FLY, ix,iy-1,iz)) +
	     fd1z[iz] * (MRC_F3(f,_FLZ, ix,iy,iz) - MRC_F3(f,_FLZ, ix,iy,iz-1)));
  } mrc_fld_foreach_end;
}

static void
pushpp_c(struct ggcm_mhd *mhd, float dt, int m)
{
  struct mrc_fld *f = mhd->fld;
  float *fd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FD1);
  float *fd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FD1);
  float *fd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FD1);

  float dth = -.5f * dt;
  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    float fpx = fd1x[ix] * (MRC_F3(f, _PP, ix+1,iy,iz) - MRC_F3(f, _PP, ix-1,iy,iz));
    float fpy = fd1y[iy] * (MRC_F3(f, _PP, ix,iy+1,iz) - MRC_F3(f, _PP, ix,iy-1,iz));
    float fpz = fd1z[iz] * (MRC_F3(f, _PP, ix,iy,iz+1) - MRC_F3(f, _PP, ix,iy,iz-1));
    float z = dth * MRC_F3(f,_ZMASK, ix,iy,iz);
    MRC_F3(f, m + _RV1X, ix,iy,iz) += z * fpx;
    MRC_F3(f, m + _RV1Y, ix,iy,iz) += z * fpy;
    MRC_F3(f, m + _RV1Z, ix,iy,iz) += z * fpz;
  } mrc_fld_foreach_end;
}

static void
vgrs(struct mrc_fld *f, int m, float s)
{
  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    MRC_F3(f, m, ix,iy,iz) = s;
  } mrc_fld_foreach_end;
}

static void
vgrv(struct mrc_fld *f, int m_to, int m_from)
{
  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    MRC_F3(f, m_to, ix,iy,iz) = MRC_F3(f, m_from, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static inline void
limit1a(struct mrc_fld *f, int m, int ix, int iy, int iz, int IX, int IY, int IZ, int C)
{
  const float reps = 0.003;
  const float seps = -0.001;
  const float teps = 1.e-25;

  // Harten/Zwas type switch
  float aa = MRC_F3(f, m, ix,iy,iz);
  float a1 = MRC_F3(f, m, ix+IX,iy+IY,iz+IZ);
  float a2 = MRC_F3(f, m, ix-IX,iy-IY,iz-IZ);
  float d1 = aa - a2;
  float d2 = a1 - aa;
  float s1 = fabsf(d1);
  float s2 = fabsf(d2);
  float f1 = fabsf(a1) + fabsf(a2) + fabsf(aa);
  float s5 = s1 + s2 + reps*f1 + teps;
  float r3 = fabsf(s1 - s2) / s5; // edge condition
  float f2 = seps * f1 * f1;
  if (d1 * d2 < f2) {
    r3 = 1.f;
  }
  r3 = r3 * r3;
  r3 = r3 * r3;
  r3 = fminf(2.f * r3, 1.);
  MRC_F3(f, C, ix   ,iy   ,iz   ) = fmaxf(MRC_F3(f, C, ix   ,iy   ,iz   ), r3);
  MRC_F3(f, C, ix-IX,iy-IY,iz-IZ) = fmaxf(MRC_F3(f, C, ix-IX,iy-IY,iz-IZ), r3);
}

static void
limit1_c(struct mrc_fld *f, int m, float time, float timelo, int C)
{
  if (time < timelo) {
    vgrs(f, C + 0, 1.f);
    vgrs(f, C + 1, 1.f);
    vgrs(f, C + 2, 1.f);
    return;
  }

  mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
/* .if (limit_aspect_low) then */
/* .call lowmask(0,0,0,tl1) */
    limit1a(f, m, ix,iy,iz, 1,0,0, C + 0);
    limit1a(f, m, ix,iy,iz, 0,1,0, C + 1);
    limit1a(f, m, ix,iy,iz, 0,0,1, C + 2);
  } mrc_fld_foreach_end;
}

static void
vgfl_c(struct ggcm_mhd *mhd, int m)
{
  switch (m) {
  case _RR1:  return vgflrr_c(mhd);
  case _RV1X: return vgflrvx_c(mhd);
  case _RV1Y: return vgflrvy_c(mhd);
  case _RV1Z: return vgflrvz_c(mhd);
  case _UU1:  return vgfluu_c(mhd);
  default: assert(0);
  }
}

static void
pushfv_c(struct ggcm_mhd *mhd, int m, float dt, int m_prev, int m_curr, int m_next,
	 int limit)
{
  struct mrc_fld *f = mhd->fld;

  vgfl_c(mhd, m);
  if (limit == LIMIT_NONE) {
    fluxl_c(mhd, m_curr + m);
  } else {
    vgrv(f, _CX, _BX); vgrv(f, _CY, _BY); vgrv(f, _CY, _BY);
    limit1_c(f, m_curr + m, mhd->time, mhd->par.timelo, _CX);
    fluxb_c(mhd, m_curr + m);
  }

  pushn_c(mhd, m_prev + m, m_next + m, dt);
}

static void
currbb_c(struct ggcm_mhd *mhd, int m_curr)
{
  enum { _XJX = _BX, _XJY = _BY, _XJZ = _BZ };

  struct mrc_fld *f = mhd->fld;
  float *bd4x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD4);
  float *bd4y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD4);
  float *bd4z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD4);

  mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
    MRC_F3(f, _XJX, ix,iy,iz) =
      (MRC_F3(f, m_curr + _B1Z, ix,iy+1,iz) - MRC_F3(f, m_curr + _B1Z, ix,iy,iz)) * bd4y[iy] -
      (MRC_F3(f, m_curr + _B1Y, ix,iy,iz+1) - MRC_F3(f, m_curr + _B1Y, ix,iy,iz)) * bd4z[iz];
    MRC_F3(f, _XJY, ix,iy,iz) =
      (MRC_F3(f, m_curr + _B1X, ix,iy,iz+1) - MRC_F3(f, m_curr + _B1X, ix,iy,iz)) * bd4z[iz] -
      (MRC_F3(f, m_curr + _B1Z, ix+1,iy,iz) - MRC_F3(f, m_curr + _B1Z, ix,iy,iz)) * bd4x[ix];
    MRC_F3(f, _XJZ, ix,iy,iz) =
      (MRC_F3(f, m_curr + _B1Y, ix+1,iy,iz) - MRC_F3(f, m_curr + _B1Y, ix,iy,iz)) * bd4x[ix] -
      (MRC_F3(f, m_curr + _B1X, ix,iy+1,iz) - MRC_F3(f, m_curr + _B1X, ix,iy,iz)) * bd4y[iy];
    MRC_F3(f,_TMP1, ix,iy,iz) = .5f * (MRC_F3(f, m_curr + _B1X, ix,iy,iz) +
				       MRC_F3(f, m_curr + _B1X, ix-1,iy,iz));
    MRC_F3(f,_TMP2, ix,iy,iz) = .5f * (MRC_F3(f, m_curr + _B1Y, ix,iy,iz) +
				       MRC_F3(f, m_curr + _B1Y, ix,iy-1,iz));
    MRC_F3(f,_TMP3, ix,iy,iz) = .5f * (MRC_F3(f, m_curr + _B1Z, ix,iy,iz) +
				       MRC_F3(f, m_curr + _B1Z, ix,iy,iz-1));
  } mrc_fld_foreach_end;
}

static void
push_ej_c(struct ggcm_mhd *mhd, float dt, int m_curr, int m_next)
{
  enum { _XJX = _BX, _XJY = _BY, _XJZ = _BZ };

  currbb_c(mhd, m_curr);
	
  struct mrc_fld *f = mhd->fld;

  float s1 = .25f * dt;
  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    float z = MRC_F3(f,_ZMASK, ix,iy,iz);
    float s2 = s1 * z;
    float cx = (MRC_F3(f,_XJX, ix  ,iy  ,iz  ) +
		MRC_F3(f,_XJX, ix  ,iy-1,iz  ) +
		MRC_F3(f,_XJX, ix  ,iy  ,iz-1) +
		MRC_F3(f,_XJX, ix  ,iy-1,iz-1));
    float cy = (MRC_F3(f,_XJY, ix  ,iy  ,iz  ) +
		MRC_F3(f,_XJY, ix-1,iy  ,iz  ) +
		MRC_F3(f,_XJY, ix  ,iy  ,iz-1) +
		MRC_F3(f,_XJY, ix-1,iy  ,iz-1));
    float cz = (MRC_F3(f,_XJZ, ix  ,iy  ,iz  ) +
		MRC_F3(f,_XJZ, ix-1,iy  ,iz  ) +
		MRC_F3(f,_XJZ, ix  ,iy-1,iz  ) +
		MRC_F3(f,_XJZ, ix-1,iy-1,iz  ));
    float ffx = s2 * (cy * MRC_F3(f, _TMP3, ix,iy,iz) -
		      cz * MRC_F3(f, _TMP2, ix,iy,iz));
    float ffy = s2 * (cz * MRC_F3(f, _TMP1, ix,iy,iz) -
		      cx * MRC_F3(f, _TMP3, ix,iy,iz));
    float ffz = s2 * (cx * MRC_F3(f, _TMP2, ix,iy,iz) -
		      cy * MRC_F3(f, _TMP1, ix,iy,iz));
    float duu = (ffx * MRC_F3(f, _VX, ix,iy,iz) +
		 ffy * MRC_F3(f, _VY, ix,iy,iz) +
		 ffz * MRC_F3(f, _VZ, ix,iy,iz));

    MRC_F3(f, m_next + _RV1X, ix,iy,iz) += ffx;
    MRC_F3(f, m_next + _RV1Y, ix,iy,iz) += ffy;
    MRC_F3(f, m_next + _RV1Z, ix,iy,iz) += ffz;
    MRC_F3(f, m_next + _UU1 , ix,iy,iz) += duu;
  } mrc_fld_foreach_end;
}

static inline float
bcthy3f(float s1, float s2)
{
  if (s1 > 0.f && fabsf(s2) > REPS) {
/* .if(calce_aspect_low) then */
/* .call lowmask(IX, 0, 0,tl1) */
/* .call lowmask( 0,IY, 0,tl2) */
/* .call lowmask( 0, 0,IZ,tl3) */
/* .call lowmask(IX,IY,IZ,tl4) */
/*       tt=tt*(1.0-max(tl1,tl2,tl3,tl4)) */
    return s1 / s2;
  }
  return 0.f;
}


static void
bcthy3z_NL1(struct ggcm_mhd *mhd, int XX, int YY, int ZZ, int IX, int IY, int IZ,
	    int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	    int FF, float dt, int m_curr)
{
  struct mrc_fld *f = mhd->fld;

  float *bd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD1);
  float *bd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD1);
  float *bd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD1);

  float *bd2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD2);
  float *bd2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD2);
  float *bd2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD2);

  // d_z B_y, d_y B_z
  mrc_fld_foreach(f, ix,iy,iz, 2, 1) {
    float bd1[3] = { bd1x[ix], bd1y[iy], bd1z[iz] };

    MRC_F3(f, _TMP1, ix,iy,iz) = bd1[ZZ] * 
      (MRC_F3(f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2) - MRC_F3(f, m_curr + _B1X + YY, ix,iy,iz));
    MRC_F3(f, _TMP2, ix,iy,iz) = bd1[YY] * 
      (MRC_F3(f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1) - MRC_F3(f, m_curr + _B1X + ZZ, ix,iy,iz));
  } mrc_fld_foreach_end;

  // harmonic average if same sign
  mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
    float s1, s2;
    s1 = MRC_F3(f, _TMP1, ix,iy,iz) * MRC_F3(f, _TMP1, ix-JX2,iy-JY2,iz-JZ2);
    s2 = MRC_F3(f, _TMP1, ix,iy,iz) + MRC_F3(f, _TMP1, ix-JX2,iy-JY2,iz-JZ2);
    MRC_F3(f, _TMP3, ix,iy,iz) = bcthy3f(s1, s2);
    s1 = MRC_F3(f, _TMP2, ix,iy,iz) * MRC_F3(f, _TMP2, ix-JX1,iy-JY1,iz-JZ1);
    s2 = MRC_F3(f, _TMP2, ix,iy,iz) + MRC_F3(f, _TMP2, ix-JX1,iy-JY1,iz-JZ1);
    MRC_F3(f, _TMP4, ix,iy,iz) = bcthy3f(s1, s2);
  } mrc_fld_foreach_end;

  // edge centered velocity
  mrc_fld_foreach(f, ix,iy,iz, 1, 0) {
    float vvYY = .25f * (MRC_F3(f, _VX + YY, ix   ,iy   ,iz   ) + 
			 MRC_F3(f, _VX + YY, ix   ,iy+IY,iz+IZ) +
			 MRC_F3(f, _VX + YY, ix+IX,iy   ,iz+IZ) +
			 MRC_F3(f, _VX + YY, ix+IX,iy+IY,iz   ));
    MRC_F3(f, _TMP1, ix,iy,iz) = vvYY; /* - d_i * vcurrYY */
    float vvZZ = .25f * (MRC_F3(f, _VX + ZZ, ix   ,iy   ,iz   ) + 
			 MRC_F3(f, _VX + ZZ, ix   ,iy+IY,iz+IZ) +
			 MRC_F3(f, _VX + ZZ, ix+IX,iy   ,iz+IZ) +
			 MRC_F3(f, _VX + ZZ, ix+IX,iy+IY,iz   ));
    MRC_F3(f, _TMP2, ix,iy,iz) = vvZZ; /* - d_i * vcurrZZ */
  } mrc_fld_foreach_end;

  float diffmul=1.0;
  if (mhd->time < 600.f) { // no anomalous res at startup
    diffmul = 0.f;
  }

  // edge centered E = - v x B (+ dissipation)
  mrc_fld_foreach(f, ix,iy,iz, 1, 0) {
    float bd2[3] = { bd2x[ix], bd2y[iy], bd2z[iz] };
    float bd2p[3] = { bd2x[ix+1], bd2y[iy+1], bd2z[iz+1] };
    float e1, vv;
    vv = MRC_F3(f, _TMP1, ix,iy,iz);
    if (vv > 0.f) {
      e1 = MRC_F3(f, m_curr + _B1X + ZZ, ix,iy,iz) +
	MRC_F3(f, _TMP4, ix,iy,iz) * (bd2[YY] - dt*vv);
    } else {
      e1 = MRC_F3(f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1) -
	MRC_F3(f, _TMP4, ix+JX1,iy+JY1,iz+JZ1) * (bd2p[YY] + dt*vv);
    }
    float ttmp1 = e1 * vv;

    vv = MRC_F3(f, _TMP2, ix,iy,iz);
    if (vv > 0.f) {
      e1 = MRC_F3(f, m_curr + _B1X + YY, ix,iy,iz) +
	MRC_F3(f, _TMP3, ix,iy,iz) * (bd2[ZZ] - dt*vv);
    } else {
      e1 = MRC_F3(f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2) -
	MRC_F3(f, _TMP3, ix+JX2,iy+JY2,iz+JZ2) * (bd2p[ZZ] + dt*vv);
    }
    float ttmp2 = e1 * vv;

    float t1m = MRC_F3(f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1) - MRC_F3(f, m_curr + _B1X + ZZ, ix,iy,iz);
    float t1p = fabsf(MRC_F3(f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1)) + fabsf(MRC_F3(f, m_curr + _B1X + ZZ, ix,iy,iz));
    float t2m = MRC_F3(f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2) - MRC_F3(f, m_curr + _B1X + YY, ix,iy,iz);
    float t2p = fabsf(MRC_F3(f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2)) + fabsf(MRC_F3(f, m_curr + _B1X + YY, ix,iy,iz));
    float tp = t1p + t2p + REPS;
    float tpi = diffmul / tp;
    float d1 = sqr(t1m * tpi);
    float d2 = sqr(t2m * tpi);
    if (d1 < mhd->par.diffth) d1 = 0.;
    if (d2 < mhd->par.diffth) d2 = 0.;
    ttmp1 -= d1 * t1m * MRC_F3(f, _RMASK, ix,iy,iz);
    ttmp2 -= d2 * t2m * MRC_F3(f, _RMASK, ix,iy,iz);
    MRC_F3(f, _RESIS, ix,iy,iz) += fabsf(d1+d2) * MRC_F3(f, _ZMASK, ix,iy,iz);
    MRC_F3(f, FF, ix,iy,iz) = ttmp1 - ttmp2;
  } mrc_fld_foreach_end;
}

static void
calce_nl1_c(struct ggcm_mhd *mhd, float dt, int m_curr)
{
  bcthy3z_NL1(mhd, 0,1,2, 0,1,1, 0,1,0, 0,0,1, _FLX, dt, m_curr);
  bcthy3z_NL1(mhd, 1,2,0, 1,0,1, 0,0,1, 1,0,0, _FLY, dt, m_curr);
  bcthy3z_NL1(mhd, 2,0,1, 1,1,0, 1,0,0, 0,1,0, _FLZ, dt, m_curr);
}

static void
calce_c(struct ggcm_mhd *mhd, float dt, int m_curr)
{
  switch (mhd->par.magdiffu) {
  case MAGDIFFU_NL1:
    return calce_nl1_c(mhd, dt, m_curr);
  default:
    assert(0);
  }
}

static void
bpush_c(struct ggcm_mhd *mhd, float dt, int m_prev, int m_next)
{
  struct mrc_fld *f = mhd->fld;
  float *bd3x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD3);
  float *bd3y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD3);
  float *bd3z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD3);

  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    MRC_F3(f, m_next + _B1X, ix,iy,iz) = MRC_F3(f, m_prev + _B1X, ix,iy,iz) +
      dt * (bd3y[iy] * (MRC_F3(f,_FLZ, ix,iy,iz) - MRC_F3(f,_FLZ, ix,iy-1,iz)) -
	    bd3z[iz] * (MRC_F3(f,_FLY, ix,iy,iz) - MRC_F3(f,_FLY, ix,iy,iz-1)));
    MRC_F3(f, m_next + _B1Y, ix,iy,iz) = MRC_F3(f, m_prev + _B1Y, ix,iy,iz) +
      dt * (bd3z[iz] * (MRC_F3(f,_FLX, ix,iy,iz) - MRC_F3(f,_FLX, ix,iy,iz-1)) -
	    bd3x[ix] * (MRC_F3(f,_FLZ, ix,iy,iz) - MRC_F3(f,_FLZ, ix-1,iy,iz)));
    MRC_F3(f, m_next + _B1Z, ix,iy,iz) = MRC_F3(f, m_prev + _B1Z, ix,iy,iz) +
      dt * (bd3x[ix] * (MRC_F3(f,_FLY, ix,iy,iz) - MRC_F3(f,_FLY, ix-1,iy,iz)) -
	    bd3y[iy] * (MRC_F3(f,_FLX, ix,iy,iz) - MRC_F3(f,_FLX, ix,iy-1,iz)));
  } mrc_fld_foreach_end;
}

static void
pushstage_c(struct ggcm_mhd *mhd, float dt, int m_prev, int m_curr, int m_next,
	    int limit)
{
  rmaskn_c(mhd);

  if (limit != LIMIT_NONE) {
    struct mrc_fld *f = mhd->fld;

    vgrs(f, _BX, 0.f); vgrs(f, _BY, 0.f); vgrs(f, _BZ, 0.f);
    limit1_c(f, _PP, mhd->time, mhd->par.timelo, _BX);
    // limit2, 3
  }

  pushfv_c(mhd, _RR1 , dt, m_prev, m_curr, m_next, limit);
  pushfv_c(mhd, _RV1X, dt, m_prev, m_curr, m_next, limit);
  pushfv_c(mhd, _RV1Y, dt, m_prev, m_curr, m_next, limit);
  pushfv_c(mhd, _RV1Z, dt, m_prev, m_curr, m_next, limit);
  pushfv_c(mhd, _UU1 , dt, m_prev, m_curr, m_next, limit);

  pushpp_c(mhd, dt, m_next);

#if 0
  if (magdiffu.eq.magdiffu_nl1) call calc_resis_nl1(bxB,byB,bzB);
  if (magdiffu.eq.magdiffu_res1) call calc_resis_res1(bxB,byB,bzB);
  if (magdiffu.eq.magdiffu_const) call calc_resis_const(bxB,byB,bzB);
#endif

  push_ej_c(mhd, dt, m_curr, m_next);
  calce_c(mhd, dt, m_curr);
  bpush_c(mhd, dt, m_prev, m_next);
}

// ======================================================================
// ggcm_mhd_step subclass "c"
//
// this class will do full predictor / corrector steps,
// ie., including primvar() etc.

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_pred

static void
ggcm_mhd_step_c_pred(struct ggcm_mhd_step *step)
{
  primvar_c(step->mhd, _RR1);
  primbb_c(step->mhd, _RR1);
  zmaskn_c(step->mhd);

  float dth = .5f * step->mhd->dt;
  pushstage_c(step->mhd, dth, _RR1, _RR1, _RR2, LIMIT_NONE);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_corr

static void
ggcm_mhd_step_c_corr(struct ggcm_mhd_step *step)
{
  primvar_c(step->mhd, _RR2);
  primbb_c(step->mhd, _RR2);
  zmaskn_c(step->mhd);
  pushstage_c(step->mhd, step->mhd->dt, _RR1, _RR2, _RR1, LIMIT_1);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "c"

struct ggcm_mhd_step_ops ggcm_mhd_step_c_ops = {
  .name        = "c",
  .pred        = ggcm_mhd_step_c_pred,
  .corr        = ggcm_mhd_step_c_corr,
};
