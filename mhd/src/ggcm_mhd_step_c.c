
#include "ggcm_mhd_step_private.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_crds.h"

#include <mrc_domain.h>
#include <mrc_profile.h>
#include <mrc_fld_as_float.h>

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
  float diff_swbnd = mhd->par.diff_swbnd;
  int diff_obnd = mhd->par.diff_obnd;
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, 0, &info);

  float *fx1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FX1);

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    F3(f,_RMASK, ix,iy,iz) = 0.f;
    float xxx = fx1x[ix];
    if (xxx < diff_swbnd)
      continue;
    if (iy + info.off[1] < diff_obnd)
      continue;
    if (iz + info.off[2] < diff_obnd)
      continue;
    if (ix + info.off[0] >= gdims[0] - diff_obnd)
      continue;
    if (iy + info.off[1] >= gdims[1] - diff_obnd)
      continue;
    if (iz + info.off[2] >= gdims[2] - diff_obnd)
      continue;
    F3(f, _RMASK, ix,iy,iz) = diffco * F3(f, _ZMASK, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static void
vgflrr_c(struct ggcm_mhd *mhd)
{
  struct mrc_fld *f = mhd->fld;

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    float a = F3(f,_RR, ix,iy,iz);
    F3(f,_TMP1, ix,iy,iz) = a * F3(f,_VX, ix,iy,iz);
    F3(f,_TMP2, ix,iy,iz) = a * F3(f,_VY, ix,iy,iz);
    F3(f,_TMP3, ix,iy,iz) = a * F3(f,_VZ, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static void
vgflrvx_c(struct ggcm_mhd *mhd)
{
  struct mrc_fld *f = mhd->fld;

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    float a = F3(f,_RR, ix,iy,iz) * F3(f,_VX, ix,iy,iz);
    F3(f,_TMP1, ix,iy,iz) = a * F3(f,_VX, ix,iy,iz);
    F3(f,_TMP2, ix,iy,iz) = a * F3(f,_VY, ix,iy,iz);
    F3(f,_TMP3, ix,iy,iz) = a * F3(f,_VZ, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static void
vgflrvy_c(struct ggcm_mhd *mhd)
{
  struct mrc_fld *f = mhd->fld;

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    float a = F3(f,_RR, ix,iy,iz) * F3(f,_VY, ix,iy,iz);
    F3(f,_TMP1, ix,iy,iz) = a * F3(f,_VX, ix,iy,iz);
    F3(f,_TMP2, ix,iy,iz) = a * F3(f,_VY, ix,iy,iz);
    F3(f,_TMP3, ix,iy,iz) = a * F3(f,_VZ, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static void
vgflrvz_c(struct ggcm_mhd *mhd)
{
  struct mrc_fld *f = mhd->fld;

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    float a = F3(f,_RR, ix,iy,iz) * F3(f,_VZ, ix,iy,iz);
    F3(f,_TMP1, ix,iy,iz) = a * F3(f,_VX, ix,iy,iz);
    F3(f,_TMP2, ix,iy,iz) = a * F3(f,_VY, ix,iy,iz);
    F3(f,_TMP3, ix,iy,iz) = a * F3(f,_VZ, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static void
vgfluu_c(struct ggcm_mhd *mhd)
{
  struct mrc_fld *f = mhd->fld;

  float gamma = mhd->par.gamm;
  float s = gamma / (gamma - 1.f);
  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    float ep = s * F3(f,_PP, ix,iy,iz) +
      .5f * F3(f,_RR, ix,iy,iz) * (sqr(F3(f,_VX, ix,iy,iz)) + 
				   sqr(F3(f,_VY, ix,iy,iz)) + 
				   sqr(F3(f,_VZ, ix,iy,iz)));
    F3(f,_TMP1, ix,iy,iz) = ep * F3(f,_VX, ix,iy,iz);
    F3(f,_TMP2, ix,iy,iz) = ep * F3(f,_VY, ix,iy,iz);
    F3(f,_TMP3, ix,iy,iz) = ep * F3(f,_VZ, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static void
fluxl_c(struct ggcm_mhd *mhd, int m)
{
  struct mrc_fld *f = mhd->fld;

  mrc_fld_foreach(f, ix,iy,iz, 1, 0) {
    float aa = F3(f,m, ix,iy,iz);
    float cmsv = F3(f,_CMSV, ix,iy,iz);
    F3(f,_FLX, ix,iy,iz) =
      .5f * ((F3(f,_TMP1, ix  ,iy,iz) + F3(f,_TMP1, ix+1,iy,iz)) -
	     .5f * (F3(f,_CMSV, ix+1,iy,iz) + cmsv) * (F3(f,m, ix+1,iy,iz) - aa));
    F3(f,_FLY, ix,iy,iz) =
      .5f * ((F3(f,_TMP2, ix,iy  ,iz) + F3(f,_TMP2, ix,iy+1,iz)) -
	     .5f * (F3(f,_CMSV, ix,iy+1,iz) + cmsv) * (F3(f,m, ix,iy+1,iz) - aa));
    F3(f,_FLZ, ix,iy,iz) =
      .5f * ((F3(f,_TMP3, ix,iy,iz  ) + F3(f,_TMP3, ix,iy,iz+1)) -
	     .5f * (F3(f,_CMSV, ix,iy,iz+1) + cmsv) * (F3(f,m, ix,iy,iz+1) - aa));
  } mrc_fld_foreach_end;
}

static void
fluxb_c(struct ggcm_mhd *mhd, int m)
{
  struct mrc_fld *f = mhd->fld;

  float s1 = 1.f/12.f;
  float s7 = 7.f * s1;

  mrc_fld_foreach(f, ix,iy,iz, 1, 0) {
    float fhx = (s7 * (F3(f, _TMP1, ix  ,iy,iz) + F3(f, _TMP1, ix+1,iy,iz)) -
		 s1 * (F3(f, _TMP1, ix-1,iy,iz) + F3(f, _TMP1, ix+2,iy,iz)));
    float fhy = (s7 * (F3(f, _TMP2, ix,iy  ,iz) + F3(f, _TMP2, ix,iy+1,iz)) -
		 s1 * (F3(f, _TMP2, ix,iy-1,iz) + F3(f, _TMP2, ix,iy+2,iz)));
    float fhz = (s7 * (F3(f, _TMP3, ix,iy,iz  ) + F3(f, _TMP3, ix,iy,iz+1)) -
		 s1 * (F3(f, _TMP3, ix,iy,iz-1) + F3(f, _TMP3, ix,iy,iz+2)));

    float aa = F3(f,m, ix,iy,iz);
    float cmsv = F3(f,_CMSV, ix,iy,iz);
    float flx =
      .5f * ((F3(f,_TMP1, ix  ,iy,iz) + F3(f,_TMP1, ix+1,iy,iz)) -
	     .5f * (F3(f,_CMSV, ix+1,iy,iz) + cmsv) * (F3(f,m, ix+1,iy,iz) - aa));
    float fly =
      .5f * ((F3(f,_TMP2, ix,iy  ,iz) + F3(f,_TMP2, ix,iy+1,iz)) -
	     .5f * (F3(f,_CMSV, ix,iy+1,iz) + cmsv) * (F3(f,m, ix,iy+1,iz) - aa));
    float flz = 
      .5f * ((F3(f,_TMP3, ix,iy,iz  ) + F3(f,_TMP3, ix,iy,iz+1)) -
	     .5f * (F3(f,_CMSV, ix,iy,iz+1) + cmsv) * (F3(f,m, ix,iy,iz+1) - aa));

    float cx = F3(f, _CX, ix,iy,iz);
    F3(f, _FLX, ix,iy,iz) = cx * flx + (1.f - cx) * fhx;
    float cy = F3(f, _CY, ix,iy,iz);
    F3(f, _FLY, ix,iy,iz) = cy * fly + (1.f - cy) * fhy;
    float cz = F3(f, _CZ, ix,iy,iz);
    F3(f, _FLZ, ix,iy,iz) = cz * flz + (1.f - cz) * fhz;
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
    float s = dt * F3(f,_YMASK, ix,iy,iz);
    F3(f,mc, ix,iy,iz) = F3(f,ma, ix,iy,iz)
      - s * (fd1x[ix] * (F3(f,_FLX, ix,iy,iz) - F3(f,_FLX, ix-1,iy,iz)) +
	     fd1y[iy] * (F3(f,_FLY, ix,iy,iz) - F3(f,_FLY, ix,iy-1,iz)) +
	     fd1z[iz] * (F3(f,_FLZ, ix,iy,iz) - F3(f,_FLZ, ix,iy,iz-1)));
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
    float fpx = fd1x[ix] * (F3(f, _PP, ix+1,iy,iz) - F3(f, _PP, ix-1,iy,iz));
    float fpy = fd1y[iy] * (F3(f, _PP, ix,iy+1,iz) - F3(f, _PP, ix,iy-1,iz));
    float fpz = fd1z[iz] * (F3(f, _PP, ix,iy,iz+1) - F3(f, _PP, ix,iy,iz-1));
    float z = dth * F3(f,_ZMASK, ix,iy,iz);
    F3(f, m + _RV1X, ix,iy,iz) += z * fpx;
    F3(f, m + _RV1Y, ix,iy,iz) += z * fpy;
    F3(f, m + _RV1Z, ix,iy,iz) += z * fpz;
  } mrc_fld_foreach_end;
}

static void
vgrs(struct mrc_fld *f, int m, float s)
{
  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    F3(f, m, ix,iy,iz) = s;
  } mrc_fld_foreach_end;
}

static void
vgrv(struct mrc_fld *f, int m_to, int m_from)
{
  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    F3(f, m_to, ix,iy,iz) = F3(f, m_from, ix,iy,iz);
  } mrc_fld_foreach_end;
}

static inline void
limit1a(struct mrc_fld *f, int m, int ix, int iy, int iz, int IX, int IY, int IZ, int C)
{
  const float reps = 0.003;
  const float seps = -0.001;
  const float teps = 1.e-25;

  // Harten/Zwas type switch
  float aa = F3(f, m, ix,iy,iz);
  float a1 = F3(f, m, ix+IX,iy+IY,iz+IZ);
  float a2 = F3(f, m, ix-IX,iy-IY,iz-IZ);
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
  F3(f, C, ix   ,iy   ,iz   ) = fmaxf(F3(f, C, ix   ,iy   ,iz   ), r3);
  F3(f, C, ix-IX,iy-IY,iz-IZ) = fmaxf(F3(f, C, ix-IX,iy-IY,iz-IZ), r3);
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
    F3(f, _XJX, ix,iy,iz) =
      (F3(f, m_curr + _B1Z, ix,iy+1,iz) - F3(f, m_curr + _B1Z, ix,iy,iz)) * bd4y[iy] -
      (F3(f, m_curr + _B1Y, ix,iy,iz+1) - F3(f, m_curr + _B1Y, ix,iy,iz)) * bd4z[iz];
    F3(f, _XJY, ix,iy,iz) =
      (F3(f, m_curr + _B1X, ix,iy,iz+1) - F3(f, m_curr + _B1X, ix,iy,iz)) * bd4z[iz] -
      (F3(f, m_curr + _B1Z, ix+1,iy,iz) - F3(f, m_curr + _B1Z, ix,iy,iz)) * bd4x[ix];
    F3(f, _XJZ, ix,iy,iz) =
      (F3(f, m_curr + _B1Y, ix+1,iy,iz) - F3(f, m_curr + _B1Y, ix,iy,iz)) * bd4x[ix] -
      (F3(f, m_curr + _B1X, ix,iy+1,iz) - F3(f, m_curr + _B1X, ix,iy,iz)) * bd4y[iy];
    F3(f,_TMP1, ix,iy,iz) = .5f * (F3(f, m_curr + _B1X, ix,iy,iz) +
				   F3(f, m_curr + _B1X, ix-1,iy,iz));
    F3(f,_TMP2, ix,iy,iz) = .5f * (F3(f, m_curr + _B1Y, ix,iy,iz) +
				   F3(f, m_curr + _B1Y, ix,iy-1,iz));
    F3(f,_TMP3, ix,iy,iz) = .5f * (F3(f, m_curr + _B1Z, ix,iy,iz) +
				   F3(f, m_curr + _B1Z, ix,iy,iz-1));
  } mrc_fld_foreach_end;
}

  // bx12,by12,bz12,currx,curry,currz, tmp1,tmp2,tmp3,flx,fly,flz,zmask;
static void
curbc_c(struct ggcm_mhd *mhd, int m_curr)
{ 
  enum { _TX = _TMP1, _TY = _TMP2, _TZ = _TMP3 };
  enum { _BBX = _FLX, _BBY = _FLY, _BBZ = _FLZ };

  struct mrc_fld *f = mhd->fld;
  float *bd4x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD4);
  float *bd4y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD4);
  float *bd4z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD4);

  mrc_fld_foreach(f, ix,iy,iz, 2, 1) {
    F3(f, _TX, ix,iy,iz) =
      (F3(f, m_curr + _B1Z, ix,iy+1,iz) - F3(f, m_curr + _B1Z, ix,iy,iz)) * bd4y[iy] -
      (F3(f, m_curr + _B1Y, ix,iy,iz+1) - F3(f, m_curr + _B1Y, ix,iy,iz)) * bd4z[iz];
    F3(f, _TY, ix,iy,iz) =
      (F3(f, m_curr + _B1X, ix,iy,iz+1) - F3(f, m_curr + _B1X, ix,iy,iz)) * bd4z[iz] -
      (F3(f, m_curr + _B1Z, ix+1,iy,iz) - F3(f, m_curr + _B1Z, ix,iy,iz)) * bd4x[ix];
    F3(f, _TZ, ix,iy,iz) =
      (F3(f, m_curr + _B1Y, ix+1,iy,iz) - F3(f, m_curr + _B1Y, ix,iy,iz)) * bd4x[ix] -
      (F3(f, m_curr + _B1X, ix,iy+1,iz) - F3(f, m_curr + _B1X, ix,iy,iz)) * bd4y[iy];

    // FIXME, not needed (?), and staggering is suspicious
    F3(f, _BBX, ix,iy,iz) = .5f * (F3(f, m_curr + _B1X, ix,iy,iz) +
				   F3(f, m_curr + _B1X, ix+1,iy,iz));
    F3(f, _BBY, ix,iy,iz) = .5f * (F3(f, m_curr + _B1Y, ix,iy,iz) +
				   F3(f, m_curr + _B1Y, ix,iy+1,iz));
    F3(f, _BBZ, ix,iy,iz) = .5f * (F3(f, m_curr + _B1Z, ix,iy,iz) +
				   F3(f, m_curr + _B1Z, ix,iy,iz+1));
  } mrc_fld_foreach_end;

  mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
    float s = .25f * F3(f, _ZMASK, ix, iy, iz);
    F3(f, _CURRX, ix,iy,iz) = s * (F3(f, _TX, ix,iy  ,iz  ) + F3(f, _TX, ix,iy-1,iz  ) +
				   F3(f, _TX, ix,iy  ,iz-1) + F3(f, _TX, ix,iy-1,iz-1));
    F3(f, _CURRY, ix,iy,iz) = s * (F3(f, _TY, ix  ,iy,iz  ) + F3(f, _TY, ix-1,iy,iz  ) +
				   F3(f, _TY, ix  ,iy,iz-1) + F3(f, _TY, ix-1,iy,iz-1));
    F3(f, _CURRZ, ix,iy,iz) = s * (F3(f, _TZ, ix  ,iy  ,iz) + F3(f, _TZ, ix-1,iy  ,iz) +
				   F3(f, _TZ, ix  ,iy-1,iz) + F3(f, _TZ, ix-1,iy-1,iz));
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
    float z = F3(f,_ZMASK, ix,iy,iz);
    float s2 = s1 * z;
    float cx = (F3(f,_XJX, ix  ,iy  ,iz  ) +
		F3(f,_XJX, ix  ,iy-1,iz  ) +
		F3(f,_XJX, ix  ,iy  ,iz-1) +
		F3(f,_XJX, ix  ,iy-1,iz-1));
    float cy = (F3(f,_XJY, ix  ,iy  ,iz  ) +
		F3(f,_XJY, ix-1,iy  ,iz  ) +
		F3(f,_XJY, ix  ,iy  ,iz-1) +
		F3(f,_XJY, ix-1,iy  ,iz-1));
    float cz = (F3(f,_XJZ, ix  ,iy  ,iz  ) +
		F3(f,_XJZ, ix-1,iy  ,iz  ) +
		F3(f,_XJZ, ix  ,iy-1,iz  ) +
		F3(f,_XJZ, ix-1,iy-1,iz  ));
    float ffx = s2 * (cy * F3(f, _TMP3, ix,iy,iz) -
		      cz * F3(f, _TMP2, ix,iy,iz));
    float ffy = s2 * (cz * F3(f, _TMP1, ix,iy,iz) -
		      cx * F3(f, _TMP3, ix,iy,iz));
    float ffz = s2 * (cx * F3(f, _TMP2, ix,iy,iz) -
		      cy * F3(f, _TMP1, ix,iy,iz));
    float duu = (ffx * F3(f, _VX, ix,iy,iz) +
		 ffy * F3(f, _VY, ix,iy,iz) +
		 ffz * F3(f, _VZ, ix,iy,iz));

    F3(f, m_next + _RV1X, ix,iy,iz) += ffx;
    F3(f, m_next + _RV1Y, ix,iy,iz) += ffy;
    F3(f, m_next + _RV1Z, ix,iy,iz) += ffz;
    F3(f, m_next + _UU1 , ix,iy,iz) += duu;
  } mrc_fld_foreach_end;
}

static void
res1_const_c(struct ggcm_mhd *mhd)
{
  // resistivity comes in ohm*m
  int res1border = 4;
  float eta0i = 1.0/53.5848e6;
  float diffsphere2 = sqr(mhd->par.diffsphere);
  float diff = mhd->par.diffco * eta0i;

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, 0, &info);

  struct mrc_fld *f = mhd->fld;
  float *fx2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FX2);
  float *fx2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FX2);
  float *fx2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FX2);

  mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
    F3(f, _RESIS, ix,iy,iz) = 0.f;
    float r2 = fx2x[ix] + fx2y[iy] + fx2z[iz];
    if (r2 < diffsphere2)
      continue;
    if (iy + info.off[1] < res1border)
      continue;
    if (iz + info.off[2] < res1border)
      continue;
    if (ix + info.off[0] >= gdims[0] - res1border)
      continue;
    if (iy + info.off[1] >= gdims[1] - res1border)
      continue;
    if (iz + info.off[2] >= gdims[2] - res1border)
      continue;

    F3(f, _RESIS, ix,iy,iz) = diff;
  } mrc_fld_foreach_end;
}

static void
calc_resis_const_c(struct ggcm_mhd *mhd, int m_curr)
{
  curbc_c(mhd, m_curr);
  res1_const_c(mhd);
}

static void
calc_resis_nl1_c(struct ggcm_mhd *mhd, int m_curr)
{
  // used to zero _RESIS field, but that's not needed.
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

static inline void
calc_dz_By(struct ggcm_mhd *mhd, struct mrc_fld *f, int m_curr, int XX, int YY, int ZZ,
	   int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2)
{
  float *bd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD1);
  float *bd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD1);
  float *bd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD1);

  // d_z B_y, d_y B_z
  mrc_fld_foreach(f, ix,iy,iz, 2, 1) {
    float bd1[3] = { bd1x[ix], bd1y[iy], bd1z[iz] };

    F3(f, _TMP1, ix,iy,iz) = bd1[ZZ] * 
      (F3(f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2) - F3(f, m_curr + _B1X + YY, ix,iy,iz));
    F3(f, _TMP2, ix,iy,iz) = bd1[YY] * 
      (F3(f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1) - F3(f, m_curr + _B1X + ZZ, ix,iy,iz));
  } mrc_fld_foreach_end;
}

static void
bcthy3z_NL1(struct ggcm_mhd *mhd, int XX, int YY, int ZZ, int IX, int IY, int IZ,
	    int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	    int FF, float dt, int m_curr)
{
  struct mrc_fld *f = mhd->fld;

  float *bd2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD2);
  float *bd2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD2);
  float *bd2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD2);

  calc_dz_By(mhd, f, m_curr, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  // harmonic average if same sign
  mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
    float s1, s2;
    s1 = F3(f, _TMP1, ix,iy,iz) * F3(f, _TMP1, ix-JX2,iy-JY2,iz-JZ2);
    s2 = F3(f, _TMP1, ix,iy,iz) + F3(f, _TMP1, ix-JX2,iy-JY2,iz-JZ2);
    F3(f, _TMP3, ix,iy,iz) = bcthy3f(s1, s2);
    s1 = F3(f, _TMP2, ix,iy,iz) * F3(f, _TMP2, ix-JX1,iy-JY1,iz-JZ1);
    s2 = F3(f, _TMP2, ix,iy,iz) + F3(f, _TMP2, ix-JX1,iy-JY1,iz-JZ1);
    F3(f, _TMP4, ix,iy,iz) = bcthy3f(s1, s2);
  } mrc_fld_foreach_end;

  // edge centered velocity
  mrc_fld_foreach(f, ix,iy,iz, 1, 0) {
    float vvYY = .25f * (F3(f, _VX + YY, ix   ,iy   ,iz   ) + 
			 F3(f, _VX + YY, ix   ,iy+IY,iz+IZ) +
			 F3(f, _VX + YY, ix+IX,iy   ,iz+IZ) +
			 F3(f, _VX + YY, ix+IX,iy+IY,iz   ));
    F3(f, _TMP1, ix,iy,iz) = vvYY; /* - d_i * vcurrYY */
    float vvZZ = .25f * (F3(f, _VX + ZZ, ix   ,iy   ,iz   ) + 
			 F3(f, _VX + ZZ, ix   ,iy+IY,iz+IZ) +
			 F3(f, _VX + ZZ, ix+IX,iy   ,iz+IZ) +
			 F3(f, _VX + ZZ, ix+IX,iy+IY,iz   ));
    F3(f, _TMP2, ix,iy,iz) = vvZZ; /* - d_i * vcurrZZ */
  } mrc_fld_foreach_end;

  float diffmul=1.0;
  if (mhd->time < mhd->par.diff_timelo) { // no anomalous res at startup
    diffmul = 0.f;
  }

  // edge centered E = - v x B (+ dissipation)
  mrc_fld_foreach(f, ix,iy,iz, 1, 0) {
    float bd2[3] = { bd2x[ix], bd2y[iy], bd2z[iz] };
    float bd2p[3] = { bd2x[ix+1], bd2y[iy+1], bd2z[iz+1] };
    float e1, vv;
    vv = F3(f, _TMP1, ix,iy,iz);
    if (vv > 0.f) {
      e1 = F3(f, m_curr + _B1X + ZZ, ix,iy,iz) +
	F3(f, _TMP4, ix,iy,iz) * (bd2[YY] - dt*vv);
    } else {
      e1 = F3(f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1) -
	F3(f, _TMP4, ix+JX1,iy+JY1,iz+JZ1) * (bd2p[YY] + dt*vv);
    }
    float ttmp1 = e1 * vv;

    vv = F3(f, _TMP2, ix,iy,iz);
    if (vv > 0.f) {
      e1 = F3(f, m_curr + _B1X + YY, ix,iy,iz) +
	F3(f, _TMP3, ix,iy,iz) * (bd2[ZZ] - dt*vv);
    } else {
      e1 = F3(f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2) -
	F3(f, _TMP3, ix+JX2,iy+JY2,iz+JZ2) * (bd2p[ZZ] + dt*vv);
    }
    float ttmp2 = e1 * vv;

    float t1m = F3(f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1) - F3(f, m_curr + _B1X + ZZ, ix,iy,iz);
    float t1p = fabsf(F3(f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1)) + fabsf(F3(f, m_curr + _B1X + ZZ, ix,iy,iz));
    float t2m = F3(f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2) - F3(f, m_curr + _B1X + YY, ix,iy,iz);
    float t2p = fabsf(F3(f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2)) + fabsf(F3(f, m_curr + _B1X + YY, ix,iy,iz));
    float tp = t1p + t2p + REPS;
    float tpi = diffmul / tp;
    float d1 = sqr(t1m * tpi);
    float d2 = sqr(t2m * tpi);
    if (d1 < mhd->par.diffth) d1 = 0.;
    if (d2 < mhd->par.diffth) d2 = 0.;
    ttmp1 -= d1 * t1m * F3(f, _RMASK, ix,iy,iz);
    ttmp2 -= d2 * t2m * F3(f, _RMASK, ix,iy,iz);
    F3(f, _RESIS, ix,iy,iz) += fabsf(d1+d2) * F3(f, _ZMASK, ix,iy,iz);
    F3(f, FF, ix,iy,iz) = ttmp1 - ttmp2;
  } mrc_fld_foreach_end;
}

static void
bcthy3z_const(struct ggcm_mhd *mhd, int XX, int YY, int ZZ, int IX, int IY, int IZ,
	      int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	      int FF, float dt, int m_curr)
{
  struct mrc_fld *f = mhd->fld;

  float *bd2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD2);
  float *bd2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD2);
  float *bd2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD2);

  calc_dz_By(mhd, f, m_curr, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  // harmonic average if same sign
  mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
    float s1, s2;
    s1 = F3(f, _TMP1, ix,iy,iz) * F3(f, _TMP1, ix-JX2,iy-JY2,iz-JZ2);
    s2 = F3(f, _TMP1, ix,iy,iz) + F3(f, _TMP1, ix-JX2,iy-JY2,iz-JZ2);
    F3(f, _TMP3, ix,iy,iz) = bcthy3f(s1, s2);
    s1 = F3(f, _TMP2, ix,iy,iz) * F3(f, _TMP2, ix-JX1,iy-JY1,iz-JZ1);
    s2 = F3(f, _TMP2, ix,iy,iz) + F3(f, _TMP2, ix-JX1,iy-JY1,iz-JZ1);
    F3(f, _TMP4, ix,iy,iz) = bcthy3f(s1, s2);
  } mrc_fld_foreach_end;

  // edge centered velocity
  mrc_fld_foreach(f, ix,iy,iz, 1, 0) {
    float vvYY = .25f * (F3(f, _VX + YY, ix   ,iy   ,iz   ) + 
			 F3(f, _VX + YY, ix   ,iy+IY,iz+IZ) +
			 F3(f, _VX + YY, ix+IX,iy   ,iz+IZ) +
			 F3(f, _VX + YY, ix+IX,iy+IY,iz   ));
    F3(f, _TMP1, ix,iy,iz) = vvYY; /* - d_i * vcurrYY */
    float vvZZ = .25f * (F3(f, _VX + ZZ, ix   ,iy   ,iz   ) + 
			 F3(f, _VX + ZZ, ix   ,iy+IY,iz+IZ) +
			 F3(f, _VX + ZZ, ix+IX,iy   ,iz+IZ) +
			 F3(f, _VX + ZZ, ix+IX,iy+IY,iz   ));
    F3(f, _TMP2, ix,iy,iz) = vvZZ; /* - d_i * vcurrZZ */
  } mrc_fld_foreach_end;

  // edge centered E = - v x B (+ dissipation)
  mrc_fld_foreach(f, ix,iy,iz, 1, 0) {
    float bd2[3] = { bd2x[ix], bd2y[iy], bd2z[iz] };
    float bd2p[3] = { bd2x[ix+1], bd2y[iy+1], bd2z[iz+1] };
    float e1, vv;
    vv = F3(f, _TMP1, ix,iy,iz);
    if (vv > 0.f) {
      e1 = F3(f, m_curr + _B1X + ZZ, ix,iy,iz) +
	F3(f, _TMP4, ix,iy,iz) * (bd2[YY] - dt*vv);
    } else {
      e1 = F3(f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1) -
	F3(f, _TMP4, ix+JX1,iy+JY1,iz+JZ1) * (bd2p[YY] + dt*vv);
    }
    float ttmp1 = e1 * vv;

    vv = F3(f, _TMP2, ix,iy,iz);
    if (vv > 0.f) {
      e1 = F3(f, m_curr + _B1X + YY, ix,iy,iz) +
	F3(f, _TMP3, ix,iy,iz) * (bd2[ZZ] - dt*vv);
    } else {
      e1 = F3(f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2) -
	F3(f, _TMP3, ix+JX2,iy+JY2,iz+JZ2) * (bd2p[ZZ] + dt*vv);
    }
    float ttmp2 = e1 * vv;

    float vcurrXX = .25f * (F3(f, _CURRX + XX, ix   ,iy   ,iz   ) + 
			    F3(f, _CURRX + XX, ix   ,iy+IY,iz+IZ) +
			    F3(f, _CURRX + XX, ix+IX,iy   ,iz+IZ) +
			    F3(f, _CURRX + XX, ix+IX,iy+IY,iz   ));
    float vresis = .25f * (F3(f, _RESIS, ix   ,iy   ,iz   ) + 
			   F3(f, _RESIS, ix   ,iy+IY,iz+IZ) +
			   F3(f, _RESIS, ix+IX,iy   ,iz+IZ) +
			   F3(f, _RESIS, ix+IX,iy+IY,iz   ));
    F3(f, FF, ix,iy,iz) = ttmp1 - ttmp2 - vresis * vcurrXX;
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
calce_const_c(struct ggcm_mhd *mhd, float dt, int m_curr)
{
  bcthy3z_const(mhd, 0,1,2, 0,1,1, 0,1,0, 0,0,1, _FLX, dt, m_curr);
  bcthy3z_const(mhd, 1,2,0, 1,0,1, 0,0,1, 1,0,0, _FLY, dt, m_curr);
  bcthy3z_const(mhd, 2,0,1, 1,1,0, 1,0,0, 0,1,0, _FLZ, dt, m_curr);
}

static void
calce_c(struct ggcm_mhd *mhd, float dt, int m_curr)
{
  switch (mhd->par.magdiffu) {
  case MAGDIFFU_NL1:
    return calce_nl1_c(mhd, dt, m_curr);
  case MAGDIFFU_CONST:
    return calce_const_c(mhd, dt, m_curr);
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
    F3(f, m_next + _B1X, ix,iy,iz) = F3(f, m_prev + _B1X, ix,iy,iz) +
      dt * (bd3y[iy] * (F3(f,_FLZ, ix,iy,iz) - F3(f,_FLZ, ix,iy-1,iz)) -
	    bd3z[iz] * (F3(f,_FLY, ix,iy,iz) - F3(f,_FLY, ix,iy,iz-1)));
    F3(f, m_next + _B1Y, ix,iy,iz) = F3(f, m_prev + _B1Y, ix,iy,iz) +
      dt * (bd3z[iz] * (F3(f,_FLX, ix,iy,iz) - F3(f,_FLX, ix,iy,iz-1)) -
	    bd3x[ix] * (F3(f,_FLZ, ix,iy,iz) - F3(f,_FLZ, ix-1,iy,iz)));
    F3(f, m_next + _B1Z, ix,iy,iz) = F3(f, m_prev + _B1Z, ix,iy,iz) +
      dt * (bd3x[ix] * (F3(f,_FLY, ix,iy,iz) - F3(f,_FLY, ix-1,iy,iz)) -
	    bd3y[iy] * (F3(f,_FLX, ix,iy,iz) - F3(f,_FLX, ix,iy-1,iz)));
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

  switch (mhd->par.magdiffu) {
  case MAGDIFFU_NL1:
    calc_resis_nl1_c(mhd, m_curr);
    break;
  case MAGDIFFU_CONST:
    calc_resis_const_c(mhd, m_curr);
    break;
  default:
    assert(0);
  }

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
  static int PR;
  if (!PR) {
    PR = prof_register("pred_c", 1., 0, 0);
  }
  prof_start(PR);
  pushstage_c(step->mhd, dth, _RR1, _RR1, _RR2, LIMIT_NONE);
  prof_stop(PR);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_corr

static void
ggcm_mhd_step_c_corr(struct ggcm_mhd_step *step)
{
  primvar_c(step->mhd, _RR2);
  primbb_c(step->mhd, _RR2);
  //  zmaskn_c(step->mhd);

  pushstage_c(step->mhd, step->mhd->dt, _RR1, _RR2, _RR1, LIMIT_1);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "c"

struct ggcm_mhd_step_ops ggcm_mhd_step_c_ops = {
  .name        = "c",
  .mhd_type    = MT_SEMI_CONSERVATIVE_GGCM,
  .fld_type    = FLD_TYPE,
  .nr_ghosts   = 2,
  .pred        = ggcm_mhd_step_c_pred,
  .corr        = ggcm_mhd_step_c_corr,
  .run         = ggcm_mhd_step_run_predcorr,
};
