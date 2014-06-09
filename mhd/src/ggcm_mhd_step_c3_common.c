
#include "ggcm_mhd_step_private.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_crds.h"

#include <mrc_domain.h>
#include <mrc_profile.h>

#include <math.h>

// FIXME, mv to right place
static void
mrc_fld_copy_range(struct mrc_fld *to, struct mrc_fld *from, int mb, int me)
{
  assert(to->_nr_ghosts == from->_nr_ghosts);
  int bnd = to->_nr_ghosts;

   mrc_fld_foreach(to, i,j,k, bnd, bnd) {
    for (int m = mb; m < me; m++) {
      F3(to, m, i,j,k) = F3(from, m, i,j,k);
    }
  } mrc_fld_foreach_end;
}

// ======================================================================
// ggcm_mhd_step subclass "c3"

struct ggcm_mhd_step_c3 {
  struct mrc_fld *x_half;
  struct mrc_fld *prim;
  struct mrc_fld *fluxes[3];

  struct mrc_fld *masks;
  struct mrc_fld *b;
  struct mrc_fld *c;
  struct mrc_fld *tmp;
  struct mrc_fld *curr;
  struct mrc_fld *resis;
  struct mrc_fld *E;
};

#define ggcm_mhd_step_c3(step) mrc_to_subobj(step, struct ggcm_mhd_step_c3)

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
setup_mrc_fld_3d(struct mrc_fld *x, struct mrc_fld *tmpl, int nr_comps)
{
  mrc_fld_set_type(x, FLD_TYPE);
  mrc_fld_set_param_obj(x, "domain", tmpl->_domain);
  mrc_fld_set_param_int(x, "nr_spatial_dims", 3);
  mrc_fld_set_param_int(x, "nr_comps", nr_comps);
  mrc_fld_set_param_int(x, "nr_ghosts", tmpl->_nr_ghosts);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_setup

static void
ggcm_mhd_step_c_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;

  assert(mhd);

  setup_mrc_fld_3d(sub->x_half, mhd->fld, 8);
  mrc_fld_dict_add_int(sub->x_half, "mhd_type", ggcm_mhd_step_mhd_type(step));
  setup_mrc_fld_3d(sub->prim, mhd->fld, _VZ + 1);
  for (int d = 0; d < 3; d++) {
    setup_mrc_fld_3d(sub->fluxes[d], mhd->fld, 5);
  }
  setup_mrc_fld_3d(sub->tmp , mhd->fld, 4);
  setup_mrc_fld_3d(sub->b   , mhd->fld, 3);
  setup_mrc_fld_3d(sub->c   , mhd->fld, 3);
  setup_mrc_fld_3d(sub->E   , mhd->fld, 3);
  setup_mrc_fld_3d(sub->curr, mhd->fld, 3);
  setup_mrc_fld_3d(sub->resis,mhd->fld, 1);

  sub->masks = mhd->fld;

  ggcm_mhd_step_setup_member_objs_sub(step);
  ggcm_mhd_step_setup_super(step);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_primvar

static void
ggcm_mhd_step_c_primvar(struct ggcm_mhd_step *step, struct mrc_fld *prim,
			struct mrc_fld *x)
{
  mrc_fld_data_t gamm = step->mhd->par.gamm;
  mrc_fld_data_t s = gamm - 1.f;

  mrc_fld_foreach(x, i,j,k, 2, 2) {
    F3(prim,_RR, i,j,k) = RR1(x, i,j,k);
    mrc_fld_data_t rri = 1.f / RR1(x, i,j,k);
    F3(prim,_VX, i,j,k) = rri * RV1X(x, i,j,k);
    F3(prim,_VY, i,j,k) = rri * RV1Y(x, i,j,k);
    F3(prim,_VZ, i,j,k) = rri * RV1Z(x, i,j,k);
    mrc_fld_data_t rvv =
      F3(prim,_VX, i,j,k) * RV1X(x, i,j,k) +
      F3(prim,_VY, i,j,k) * RV1Y(x, i,j,k) +
      F3(prim,_VZ, i,j,k) * RV1Z(x, i,j,k);
    F3(prim,_PP, i,j,k) = s * (UU1(x, i,j,k) - .5f * rvv);
    mrc_fld_data_t cs2 = mrc_fld_max(gamm * F3(prim,_PP, i,j,k) * rri, 0.f);
    F3(prim,_CMSV, i,j,k) = sqrtf(rvv * rri) + sqrtf(cs2);
  } mrc_fld_foreach_end;
}

// ======================================================================

static void
rmaskn_c(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *masks = sub->masks;

  mrc_fld_data_t diffco = mhd->par.diffco;
  mrc_fld_data_t diff_swbnd = mhd->par.diff_swbnd;
  int diff_obnd = mhd->par.diff_obnd;
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, 0, &info);

  float *fx1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FX1);

  mrc_fld_foreach(masks, i,j,k, 2, 2) {
    F3(masks, _RMASK, i,j,k) = 0.f;
    mrc_fld_data_t xxx = fx1x[i];
    if (xxx < diff_swbnd)
      continue;
    if (j + info.off[1] < diff_obnd)
      continue;
    if (k + info.off[2] < diff_obnd)
      continue;
    if (i + info.off[0] >= gdims[0] - diff_obnd)
      continue;
    if (j + info.off[1] >= gdims[1] - diff_obnd)
      continue;
    if (k + info.off[2] >= gdims[2] - diff_obnd)
      continue;
    F3(masks, _RMASK, i,j,k) = diffco * F3(masks, _ZMASK, i,j,k);
  } mrc_fld_foreach_end;
}

static void
vgflrr_c(struct ggcm_mhd *mhd, struct mrc_fld *tmp,
	 struct mrc_fld *prim)
{
  mrc_fld_foreach(tmp, i,j,k, 2, 2) {
    mrc_fld_data_t a = F3(prim, _RR, i,j,k);
    F3(tmp, 0, i,j,k) = a * F3(prim, _VX, i,j,k);
    F3(tmp, 1, i,j,k) = a * F3(prim, _VY, i,j,k);
    F3(tmp, 2, i,j,k) = a * F3(prim, _VZ, i,j,k);
  } mrc_fld_foreach_end;
}

static void
vgflrvx_c(struct ggcm_mhd *mhd, struct mrc_fld *tmp,
	 struct mrc_fld *prim)
{
  mrc_fld_foreach(tmp, i,j,k, 2, 2) {
    mrc_fld_data_t a = F3(prim, _RR, i,j,k) * F3(prim, _VX, i,j,k);
    F3(tmp, 0, i,j,k) = a * F3(prim, _VX, i,j,k);
    F3(tmp, 1, i,j,k) = a * F3(prim, _VY, i,j,k);
    F3(tmp, 2, i,j,k) = a * F3(prim, _VZ, i,j,k);
  } mrc_fld_foreach_end;
}

static void
vgflrvy_c(struct ggcm_mhd *mhd, struct mrc_fld *tmp,
	 struct mrc_fld *prim)
{
  mrc_fld_foreach(tmp, i,j,k, 2, 2) {
    mrc_fld_data_t a = F3(prim, _RR, i,j,k) * F3(prim, _VY, i,j,k);
    F3(tmp, 0, i,j,k) = a * F3(prim, _VX, i,j,k);
    F3(tmp, 1, i,j,k) = a * F3(prim, _VY, i,j,k);
    F3(tmp, 2, i,j,k) = a * F3(prim, _VZ, i,j,k);
  } mrc_fld_foreach_end;
}

static void
vgflrvz_c(struct ggcm_mhd *mhd, struct mrc_fld *tmp,
	 struct mrc_fld *prim)
{
  mrc_fld_foreach(tmp, i,j,k, 2, 2) {
    mrc_fld_data_t a = F3(prim, _RR, i,j,k) * F3(prim, _VZ, i,j,k);
    F3(tmp, 0, i,j,k) = a * F3(prim, _VX, i,j,k);
    F3(tmp, 1, i,j,k) = a * F3(prim, _VY, i,j,k);
    F3(tmp, 2, i,j,k) = a * F3(prim, _VZ, i,j,k);
  } mrc_fld_foreach_end;
}

static void
vgfluu_c(struct ggcm_mhd *mhd, struct mrc_fld *tmp,
	 struct mrc_fld *prim)
{
  mrc_fld_data_t gamma = mhd->par.gamm;
  mrc_fld_data_t s = gamma / (gamma - 1.f);
  mrc_fld_foreach(tmp, i,j,k, 2, 2) {
    mrc_fld_data_t ep = s * F3(prim, _PP, i,j,k) +
      .5f * F3(prim, _RR, i,j,k) * (sqr(F3(prim, _VX, i,j,k)) + 
				    sqr(F3(prim, _VY, i,j,k)) + 
				    sqr(F3(prim, _VZ, i,j,k)));
    F3(tmp, 0, i,j,k) = ep * F3(prim, _VX, i,j,k);
    F3(tmp, 1, i,j,k) = ep * F3(prim, _VY, i,j,k);
    F3(tmp, 2, i,j,k) = ep * F3(prim, _VZ, i,j,k);
  } mrc_fld_foreach_end;
}

static void
fluxl_c(struct ggcm_mhd *mhd, struct mrc_fld **fluxes, struct mrc_fld *tmp,
	struct mrc_fld *x, int m, struct mrc_fld *prim)
{
  mrc_fld_foreach(fluxes[0], i,j,k, 1, 0) {
    mrc_fld_data_t aa = F3(x, m, i,j,k);
    mrc_fld_data_t cmsv = F3(prim, _CMSV, i,j,k);
    F3(fluxes[0], m, i,j,k) =
      .5f * ((F3(tmp, 0, i  ,j,k) + F3(tmp, 0, i+1,j,k)) -
	     .5f * (F3(prim, _CMSV, i+1,j,k) + cmsv) * (F3(x, m, i+1,j,k) - aa));
    F3(fluxes[1], m, i,j,k) =
      .5f * ((F3(tmp, 1, i,j  ,k) + F3(tmp, 1, i,j+1,k)) -
	     .5f * (F3(prim, _CMSV, i,j+1,k) + cmsv) * (F3(x, m, i,j+1,k) - aa));
    F3(fluxes[2], m, i,j,k) =
      .5f * ((F3(tmp, 2, i,j,k  ) + F3(tmp, 2, i,j,k+1)) -
	     .5f * (F3(prim, _CMSV, i,j,k+1) + cmsv) * (F3(x, m, i,j,k+1) - aa));
  } mrc_fld_foreach_end;
}

static void
fluxb_c(struct ggcm_mhd *mhd, struct mrc_fld **fluxes, struct mrc_fld *tmp,
	struct mrc_fld *x, int m, struct mrc_fld *prim, struct mrc_fld *c)
{
  mrc_fld_data_t s1 = 1.f/12.f;
  mrc_fld_data_t s7 = 7.f * s1;

  mrc_fld_foreach(fluxes[0], i,j,k, 1, 0) {
    mrc_fld_data_t fhx = (s7 * (F3(tmp, 0, i  ,j,k) + F3(tmp, 0, i+1,j,k)) -
			  s1 * (F3(tmp, 0, i-1,j,k) + F3(tmp, 0, i+2,j,k)));
    mrc_fld_data_t fhy = (s7 * (F3(tmp, 1, i,j  ,k) + F3(tmp, 1, i,j+1,k)) -
			  s1 * (F3(tmp, 1, i,j-1,k) + F3(tmp, 1, i,j+2,k)));
    mrc_fld_data_t fhz = (s7 * (F3(tmp, 2, i,j,k  ) + F3(tmp, 2, i,j,k+1)) -
			  s1 * (F3(tmp, 2, i,j,k-1) + F3(tmp, 2, i,j,k+2)));

    mrc_fld_data_t aa = F3(x, m, i,j,k);
    mrc_fld_data_t cmsv = F3(prim, _CMSV, i,j,k);
    mrc_fld_data_t flx =
      .5f * ((F3(tmp, 0, i  ,j,k) + F3(tmp, 0, i+1,j,k)) -
	     .5f * (F3(prim, _CMSV, i+1,j,k) + cmsv) * (F3(x, m, i+1,j,k) - aa));
    mrc_fld_data_t fly =
      .5f * ((F3(tmp, 1, i,j  ,k) + F3(tmp, 1, i,j+1,k)) -
	     .5f * (F3(prim, _CMSV, i,j+1,k) + cmsv) * (F3(x, m, i,j+1,k) - aa));
    mrc_fld_data_t flz = 
      .5f * ((F3(tmp, 2, i,j,k  ) + F3(tmp, 2, i,j,k+1)) -
	     .5f * (F3(prim, _CMSV, i,j,k+1) + cmsv) * (F3(x, m, i,j,k+1) - aa));

    mrc_fld_data_t cx = F3(c, 0, i,j,k);
    F3(fluxes[0], m, i,j,k) = cx * flx + (1.f - cx) * fhx;
    mrc_fld_data_t cy = F3(c, 1, i,j,k);
    F3(fluxes[1], m, i,j,k) = cy * fly + (1.f - cy) * fhy;
    mrc_fld_data_t cz = F3(c, 2, i,j,k);
    F3(fluxes[2], m, i,j,k) = cz * flz + (1.f - cz) * fhz;
  } mrc_fld_foreach_end;
}

static void
update_finite_volume(struct ggcm_mhd *mhd, struct mrc_fld *x, struct mrc_fld **fluxes,
		     struct mrc_fld *masks, mrc_fld_data_t dt)
{
  float *fd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FD1);
  float *fd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FD1);
  float *fd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FD1);

  mrc_fld_foreach(x, i,j,k, 0, 0) {
    mrc_fld_data_t s = dt * F3(masks, _YMASK, i,j,k);
    for (int m = 0; m < 5; m++) {
      F3(x, m, i,j,k) +=
	- s * (fd1x[i] * (F3(fluxes[0], m, i,j,k) - F3(fluxes[0], m, i-1,j,k)) +
	       fd1y[j] * (F3(fluxes[1], m, i,j,k) - F3(fluxes[1], m, i,j-1,k)) +
	       fd1z[k] * (F3(fluxes[2], m, i,j,k) - F3(fluxes[2], m, i,j,k-1)));
    }
  } mrc_fld_foreach_end;
}

static void
pushpp_c(struct ggcm_mhd_step *step, mrc_fld_data_t dt, struct mrc_fld *x,
	 struct mrc_fld *prim)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *masks = sub->masks;
  float *fd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FD1);
  float *fd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FD1);
  float *fd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FD1);

  mrc_fld_data_t dth = -.5f * dt;
  mrc_fld_foreach(x, i,j,k, 0, 0) {
    mrc_fld_data_t fpx = fd1x[i] * (F3(prim, _PP, i+1,j,k) - F3(prim, _PP, i-1,j,k));
    mrc_fld_data_t fpy = fd1y[j] * (F3(prim, _PP, i,j+1,k) - F3(prim, _PP, i,j-1,k));
    mrc_fld_data_t fpz = fd1z[k] * (F3(prim, _PP, i,j,k+1) - F3(prim, _PP, i,j,k-1));
    mrc_fld_data_t z = dth * F3(masks, _ZMASK, i,j,k);
    F3(x, _RV1X, i,j,k) += z * fpx;
    F3(x, _RV1Y, i,j,k) += z * fpy;
    F3(x, _RV1Z, i,j,k) += z * fpz;
  } mrc_fld_foreach_end;
}

static void
vgrs(struct mrc_fld *f, int m, mrc_fld_data_t s)
{
  mrc_fld_foreach(f, i,j,k, 2, 2) {
    F3(f, m, i,j,k) = s;
  } mrc_fld_foreach_end;
}

static inline void
limit1a(struct mrc_fld *x, int m, int i, int j, int k, int I, int J, int K,
	struct mrc_fld *c, int m_c)
{
  const mrc_fld_data_t reps = 0.003;
  const mrc_fld_data_t seps = -0.001;
  const mrc_fld_data_t teps = 1.e-25;

  // Harten/Zwas type switch
  mrc_fld_data_t aa = F3(x, m, i,j,k);
  mrc_fld_data_t a1 = F3(x, m, i+I,j+J,k+K);
  mrc_fld_data_t a2 = F3(x, m, i-I,j-J,k-K);
  mrc_fld_data_t d1 = aa - a2;
  mrc_fld_data_t d2 = a1 - aa;
  mrc_fld_data_t s1 = fabsf(d1);
  mrc_fld_data_t s2 = fabsf(d2);
  mrc_fld_data_t f1 = fabsf(a1) + fabsf(a2) + fabsf(aa);
  mrc_fld_data_t s5 = s1 + s2 + reps*f1 + teps;
  mrc_fld_data_t r3 = fabsf(s1 - s2) / s5; // edge condition
  mrc_fld_data_t f2 = seps * f1 * f1;
  if (d1 * d2 < f2) {
    r3 = 1.f;
  }
  r3 = r3 * r3;
  r3 = r3 * r3;
  r3 = fminf(2.f * r3, 1.);
  F3(c, m_c, i   ,j   ,k   ) = fmaxf(F3(c, m_c, i   ,j   ,k   ), r3);
  F3(c, m_c, i-I,j-J,k-K) = fmaxf(F3(c, m_c, i-I,j-J,k-K), r3);
}

static void
limit1_c(struct mrc_fld *x, int m, mrc_fld_data_t time, mrc_fld_data_t timelo,
	 struct mrc_fld *c, int m_c)
{
  if (time < timelo) {
    vgrs(c, m_c + 0, 1.f);
    vgrs(c, m_c + 1, 1.f);
    vgrs(c, m_c + 2, 1.f);
    return;
  }

  mrc_fld_foreach(c, i,j,k, 1, 1) {
/* .if (limit_aspect_low) then */
/* .call lowmask(0,0,0,tl1) */
    limit1a(x, m, i,j,k, 1,0,0, c, m_c + 0);
    limit1a(x, m, i,j,k, 0,1,0, c, m_c + 1);
    limit1a(x, m, i,j,k, 0,0,1, c, m_c + 2);
  } mrc_fld_foreach_end;
}

static void
vgfl_c(struct ggcm_mhd *mhd, int m, struct mrc_fld *tmp, struct mrc_fld *prim)
{
  switch (m) {
  case _RR1:  return vgflrr_c(mhd, tmp, prim);
  case _RV1X: return vgflrvx_c(mhd, tmp, prim);
  case _RV1Y: return vgflrvy_c(mhd, tmp, prim);
  case _RV1Z: return vgflrvz_c(mhd, tmp, prim);
  case _UU1:  return vgfluu_c(mhd, tmp, prim);
  default: assert(0);
  }
}

static void
pushfv_c(struct ggcm_mhd_step *step, struct mrc_fld **fluxes,
	 mrc_fld_data_t dt, struct mrc_fld *x_curr,
	 struct mrc_fld *x_next, int limit)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *prim = sub->prim, *b = sub->b, *c = sub->c, *tmp = sub->tmp;

  if (limit == LIMIT_NONE) {
    for (int m = 0; m < 5; m++) {
      vgfl_c(mhd, m, tmp, prim);
      fluxl_c(mhd, fluxes, tmp, x_curr, m, prim);
    }
  } else {
    vgrs(b, 0, 0.f); vgrs(b, 1, 0.f); vgrs(b, 2, 0.f);
    limit1_c(prim, _PP, mhd->time, mhd->par.timelo, b, 0);
    // limit2, 3

    for (int m = 0; m < 5; m++) {
      vgfl_c(mhd, m, tmp, prim);
      mrc_fld_foreach(c, i,j,k, 2,2) {
	F3(c, 0, i,j,k) = F3(b, 0, i,j,k);
	F3(c, 1, i,j,k) = F3(b, 1, i,j,k);
	F3(c, 2, i,j,k) = F3(b, 2, i,j,k);
      } mrc_fld_foreach_end;
      
      limit1_c(x_curr, m, mhd->time, mhd->par.timelo, c, 0);
      fluxb_c(mhd, fluxes, tmp, x_curr, m, prim, c);
    }
  }
}

// ----------------------------------------------------------------------
// curr_c
//
// edge centered current density

static void
curr_c(struct ggcm_mhd *mhd, struct mrc_fld *j_ec, struct mrc_fld *x)
{
  float *bd4x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD4) - 1;
  float *bd4y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD4) - 1;
  float *bd4z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD4) - 1;

  mrc_fld_foreach(j_ec, i,j,k, 1, 2) {
    F3(j_ec, 0, i,j,k) =
      (F3(x, _B1Z, i,j,k) - F3(x, _B1Z, i,j-1,k)) * bd4y[j] -
      (F3(x, _B1Y, i,j,k) - F3(x, _B1Y, i,j,k-1)) * bd4z[k];
    F3(j_ec, 1, i,j,k) =
      (F3(x, _B1X, i,j,k) - F3(x, _B1X, i,j,k-1)) * bd4z[k] -
      (F3(x, _B1Z, i,j,k) - F3(x, _B1Z, i-1,j,k)) * bd4x[i];
    F3(j_ec, 2, i,j,k) =
      (F3(x, _B1Y, i,j,k) - F3(x, _B1Y, i-1,j,k)) * bd4x[i] -
      (F3(x, _B1X, i,j,k) - F3(x, _B1X, i,j-1,k)) * bd4y[j];
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// compute_B_cc
//
// cell-averaged B

static void
compute_B_cc(struct ggcm_mhd *mhd, struct mrc_fld *b_cc, struct mrc_fld *x)
{
  mrc_fld_foreach(b_cc, i,j,k, 1, 1) {
    F3(b_cc, 0, i,j,k) = .5f * (F3(x, _B1X, i  ,j,k) +
				   F3(x, _B1X, i+1,j,k));
    F3(b_cc, 1, i,j,k) = .5f * (F3(x, _B1Y, i,j  ,k) +
				   F3(x, _B1Y, i,j+1,k));
    F3(b_cc, 2, i,j,k) = .5f * (F3(x, _B1Z, i,j,k  ) +
				   F3(x, _B1Z, i,j,k+1));
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// curbc_c
//
// cell-centered j

static void
curbc_c(struct ggcm_mhd_step *step, struct mrc_fld *j_cc,
	struct mrc_fld *x)
{ 
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *masks = sub->masks;

  // get j on edges
  struct mrc_fld *j_ec = sub->tmp;
  curr_c(mhd, j_ec, x);

  // then average to cell centers
  mrc_fld_foreach(j_cc, i,j,k, 1, 1) {
    mrc_fld_data_t s = .25f * F3(masks, _ZMASK, i, j, k);
    F3(j_cc, 0, i,j,k) = s * (F3(j_ec, 0, i,j+1,k+1) + F3(j_ec, 0, i,j,k+1) +
			      F3(j_ec, 0, i,j+1,k  ) + F3(j_ec, 0, i,j,k  ));
    F3(j_cc, 1, i,j,k) = s * (F3(j_ec, 1, i+1,j,k+1) + F3(j_ec, 1, i,j,k+1) +
			      F3(j_ec, 1, i+1,j,k  ) + F3(j_ec, 1, i,j,k  ));
    F3(j_cc, 2, i,j,k) = s * (F3(j_ec, 2, i+1,j+1,k) + F3(j_ec, 2, i,j+1,k) +
			      F3(j_ec, 2, i+1,j  ,k) + F3(j_ec, 2, i,j  ,k));
  } mrc_fld_foreach_end;
}

static void
push_ej_c(struct ggcm_mhd_step *step, mrc_fld_data_t dt, struct mrc_fld *x_curr,
	  struct mrc_fld *x_next)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *prim = sub->prim;
  struct mrc_fld *j_ec = sub->c, *b_cc = sub->b, *masks = sub->masks;

  curr_c(mhd, j_ec, x_curr);
  compute_B_cc(mhd, b_cc, x_curr);
	
  mrc_fld_data_t s1 = .25f * dt;
  mrc_fld_foreach(x_next, i,j,k, 0, 0) {
    mrc_fld_data_t z = F3(masks, _ZMASK, i,j,k);
    mrc_fld_data_t s2 = s1 * z;
    mrc_fld_data_t cx = (F3(j_ec, 0, i  ,j+1,k+1) + F3(j_ec, 0, i  ,j  ,k+1) +
			 F3(j_ec, 0, i  ,j+1,k  ) + F3(j_ec, 0, i  ,j  ,k  ));
    mrc_fld_data_t cy = (F3(j_ec, 1, i+1,j  ,k+1) + F3(j_ec, 1, i  ,j  ,k+1) +
			 F3(j_ec, 1, i+1,j  ,k  ) + F3(j_ec, 1, i  ,j  ,k  ));
    mrc_fld_data_t cz = (F3(j_ec, 2, i+1,j+1,k  ) + F3(j_ec, 2, i  ,j+1,k  ) +
			 F3(j_ec, 2, i+1,j  ,k  ) + F3(j_ec, 2, i  ,j  ,k  ));
    mrc_fld_data_t ffx = s2 * (cy * F3(b_cc, 2, i,j,k) - cz * F3(b_cc, 1, i,j,k));
    mrc_fld_data_t ffy = s2 * (cz * F3(b_cc, 0, i,j,k) - cx * F3(b_cc, 2, i,j,k));
    mrc_fld_data_t ffz = s2 * (cx * F3(b_cc, 1, i,j,k) - cy * F3(b_cc, 0, i,j,k));
    mrc_fld_data_t duu = (ffx * F3(prim, _VX, i,j,k) +
			  ffy * F3(prim, _VY, i,j,k) +
			  ffz * F3(prim, _VZ, i,j,k));

    F3(x_next, _RV1X, i,j,k) += ffx;
    F3(x_next, _RV1Y, i,j,k) += ffy;
    F3(x_next, _RV1Z, i,j,k) += ffz;
    F3(x_next, _UU1 , i,j,k) += duu;
  } mrc_fld_foreach_end;
}

static void
res1_const_c(struct ggcm_mhd *mhd, struct mrc_fld *resis)
{
  // resistivity comes in ohm*m
  int diff_obnd = mhd->par.diff_obnd;
  mrc_fld_data_t eta0i = 1.0/53.5848e6;
  mrc_fld_data_t diffsphere2 = sqr(mhd->par.diffsphere);
  mrc_fld_data_t diff = mhd->par.diffco * eta0i;

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, 0, &info);

  float *fx2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FX2);
  float *fx2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FX2);
  float *fx2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FX2);

  mrc_fld_foreach(resis, i,j,k, 1, 1) {
    F3(resis, 0, i,j,k) = 0.f;
    mrc_fld_data_t r2 = fx2x[i] + fx2y[j] + fx2z[k];
    if (r2 < diffsphere2)
      continue;
    if (j + info.off[1] < diff_obnd)
      continue;
    if (k + info.off[2] < diff_obnd)
      continue;
    if (i + info.off[0] >= gdims[0] - diff_obnd)
      continue;
    if (j + info.off[1] >= gdims[1] - diff_obnd)
      continue;
    if (k + info.off[2] >= gdims[2] - diff_obnd)
      continue;

    F3(resis, 0, i,j,k) = diff;
  } mrc_fld_foreach_end;
}

static void
calc_resis_const_c(struct ggcm_mhd_step *step, struct mrc_fld *curr,
		   struct mrc_fld *resis, struct mrc_fld *x)
{
  struct ggcm_mhd *mhd = step->mhd;

  curbc_c(step, curr, x);
  res1_const_c(mhd, resis);
}

static void
calc_resis_nl1_c(struct ggcm_mhd *mhd)
{
  // used to zero _RESIS field, but that's not needed.
}

static inline mrc_fld_data_t
bcthy3f(mrc_fld_data_t s1, mrc_fld_data_t s2)
{
  if (s1 > 0.f && fabsf(s2) > REPS) {
/* .if(calce_aspect_low) then */
/* .call lowmask(I, 0, 0,tl1) */
/* .call lowmask( 0,J, 0,tl2) */
/* .call lowmask( 0, 0,K,tl3) */
/* .call lowmask(I,J,K,tl4) */
/*       tt=tt*(1.0-max(tl1,tl2,tl3,tl4)) */
    return s1 / s2;
  }
  return 0.f;
}

static inline void
calc_avg_dz_By(struct ggcm_mhd_step *step, struct mrc_fld *x, int XX, int YY, int ZZ,
	       int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *tmp = sub->tmp;

  float *bd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD1);
  float *bd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD1);
  float *bd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD1);

  // d_z B_y, d_y B_z on x edges
  mrc_fld_foreach(tmp, i,j,k, 1, 2) {
    mrc_fld_data_t bd1[3] = { bd1x[i-1], bd1y[j-1], bd1z[k-1] };

    F3(tmp, 0, i,j,k) = bd1[ZZ] * 
      (F3(x, _B1X + YY, i,j,k) - F3(x, _B1X + YY, i-JX2,j-JY2,k-JZ2));
    F3(tmp, 1, i,j,k) = bd1[YY] * 
      (F3(x, _B1X + ZZ, i,j,k) - F3(x, _B1X + ZZ, i-JX1,j-JY1,k-JZ1));
  } mrc_fld_foreach_end;

  // .5 * harmonic average if same sign
  mrc_fld_foreach(tmp, i,j,k, 1, 1) {
    mrc_fld_data_t s1, s2;
    // dz_By on y face
    s1 = F3(tmp, 0, i+JX2,j+JY2,k+JZ2) * F3(tmp, 0, i,j,k);
    s2 = F3(tmp, 0, i+JX2,j+JY2,k+JZ2) + F3(tmp, 0, i,j,k);
    F3(tmp, 2, i,j,k) = bcthy3f(s1, s2);
    // dy_Bz on z face
    s1 = F3(tmp, 1, i+JX1,j+JY1,k+JZ1) * F3(tmp, 1, i,j,k);
    s2 = F3(tmp, 1, i+JX1,j+JY1,k+JZ1) + F3(tmp, 1, i,j,k);
    F3(tmp, 3, i,j,k) = bcthy3f(s1, s2);
  } mrc_fld_foreach_end;
}

#define CC_TO_EC(f, m, i,j,k, I,J,K) \
  (.25f * (F3(f, m, i-I,j-J,k-K) +  \
	   F3(f, m, i-I,j   ,k   ) +  \
	   F3(f, m, i   ,j-J,k   ) +  \
	   F3(f, m, i   ,j   ,k-K)))

static inline void
calc_v_x_B(mrc_fld_data_t ttmp[2], struct mrc_fld *x,
	   struct mrc_fld *prim, struct mrc_fld *tmp,
	   int i, int j, int k,
	   int XX, int YY, int ZZ, int I, int J, int K,
	   int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	   float *bd2x, float *bd2y, float *bd2z, mrc_fld_data_t dt)
{
    mrc_fld_data_t bd2m[3] = { bd2x[i-1], bd2y[j-1], bd2z[k-1] };
    mrc_fld_data_t bd2[3] = { bd2x[i], bd2y[j], bd2z[k] };
    mrc_fld_data_t vbZZ;
    // edge centered velocity
    mrc_fld_data_t vvYY = CC_TO_EC(prim, _VX + YY, i,j,k, I,J,K) /* - d_i * vcurrYY */;
    if (vvYY > 0.f) {
      vbZZ = F3(x, _B1X + ZZ, i-JX1,j-JY1,k-JZ1) +
	F3(tmp, 3, i-JX1,j-JY1,k-JZ1) * (bd2m[YY] - dt*vvYY);
    } else {
      vbZZ = F3(x, _B1X + ZZ, i,j,k) -
	F3(tmp, 3, i,j,k) * (bd2[YY] + dt*vvYY);
    }
    ttmp[0] = vbZZ * vvYY;

    mrc_fld_data_t vbYY;
    // edge centered velocity
    mrc_fld_data_t vvZZ = CC_TO_EC(prim, _VX + ZZ, i,j,k, I,J,K) /* - d_i * vcurrZZ */;
    if (vvZZ > 0.f) {
      vbYY = F3(x, _B1X + YY, i-JX2,j-JY2,k-JZ2) +
	F3(tmp, 2, i-JX2,j-JY2,k-JZ2) * (bd2m[ZZ] - dt*vvZZ);
    } else {
      vbYY = F3(x, _B1X + YY, i,j,k) -
	F3(tmp, 2, i,j,k) * (bd2[ZZ] + dt*vvZZ);
    }
    ttmp[1] = vbYY * vvZZ;
}

static void
bcthy3z_NL1(struct ggcm_mhd_step *step, int XX, int YY, int ZZ, int I, int J, int K,
	    int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	    struct mrc_fld *E, mrc_fld_data_t dt, struct mrc_fld *x)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *prim = sub->prim, *tmp = sub->tmp, *masks = sub->masks;

  float *bd2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD2);
  float *bd2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD2);
  float *bd2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD2);

  calc_avg_dz_By(step, x, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  mrc_fld_data_t diffmul=1.0;
  if (mhd->time < mhd->par.diff_timelo) { // no anomalous res at startup
    diffmul = 0.f;
  }

  // edge centered E = - v x B (+ dissipation)
  mrc_fld_foreach(E, i,j,k, 1, 0) {
    mrc_fld_data_t ttmp[2];
    calc_v_x_B(ttmp, x, prim, tmp, i, j, k, XX, YY, ZZ, I, J, K,
	       JX1, JY1, JZ1, JX2, JY2, JZ2, bd2x, bd2y, bd2z, dt);
    
    mrc_fld_data_t t1m = F3(x, _B1X + ZZ, i+JX1,j+JY1,k+JZ1) - F3(x, _B1X + ZZ, i,j,k);
    mrc_fld_data_t t1p = fabsf(F3(x, _B1X + ZZ, i+JX1,j+JY1,k+JZ1)) + fabsf(F3(x, _B1X + ZZ, i,j,k));
    mrc_fld_data_t t2m = F3(x, _B1X + YY, i+JX2,j+JY2,k+JZ2) - F3(x, _B1X + YY, i,j,k);
    mrc_fld_data_t t2p = fabsf(F3(x, _B1X + YY, i+JX2,j+JY2,k+JZ2)) + fabsf(F3(x, _B1X + YY, i,j,k));
    mrc_fld_data_t tp = t1p + t2p + REPS;
    mrc_fld_data_t tpi = diffmul / tp;
    mrc_fld_data_t d1 = sqr(t1m * tpi);
    mrc_fld_data_t d2 = sqr(t2m * tpi);
    if (d1 < mhd->par.diffth) d1 = 0.;
    if (d2 < mhd->par.diffth) d2 = 0.;
    ttmp[0] -= d1 * t1m * F3(masks, _RMASK, i,j,k);
    ttmp[1] -= d2 * t2m * F3(masks, _RMASK, i,j,k);
    //    F3(f, _RESIS, i,j,k) += fabsf(d1+d2) * F3(masks, _ZMASK, i,j,k);
    F3(E, XX, i,j,k) = - (ttmp[0] - ttmp[1]);
  } mrc_fld_foreach_end;
}

static void
bcthy3z_const(struct ggcm_mhd_step *step, int XX, int YY, int ZZ, int I, int J, int K,
	      int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	      struct mrc_fld *E, mrc_fld_data_t dt, struct mrc_fld *x,
	      struct mrc_fld *curr, struct mrc_fld *resis)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *prim = sub->prim, *tmp = sub->tmp;

  float *bd2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD2);
  float *bd2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD2);
  float *bd2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD2);

  calc_avg_dz_By(step, x, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  // edge centered E = - v x B (+ dissipation)
  mrc_fld_foreach(E, i,j,k, 0, 1) {
    mrc_fld_data_t ttmp[2];
    calc_v_x_B(ttmp, x, prim, tmp, i, j, k, XX, YY, ZZ, I, J, K,
	       JX1, JY1, JZ1, JX2, JY2, JZ2, bd2x, bd2y, bd2z, dt);

    mrc_fld_data_t vcurrXX = CC_TO_EC(curr, XX, i,j,k, I,J,K);
    mrc_fld_data_t vresis = CC_TO_EC(resis, 0, i,j,k, I,J,K);
    F3(E, XX, i,j,k) = - (ttmp[0] - ttmp[1]) + vresis * vcurrXX;
  } mrc_fld_foreach_end;
}

static void
calce_nl1_c(struct ggcm_mhd_step *step, struct mrc_fld *E,
	    mrc_fld_data_t dt, struct mrc_fld *x)
{
  bcthy3z_NL1(step, 0,1,2, 0,1,1, 0,1,0, 0,0,1, E, dt, x);
  bcthy3z_NL1(step, 1,2,0, 1,0,1, 0,0,1, 1,0,0, E, dt, x);
  bcthy3z_NL1(step, 2,0,1, 1,1,0, 1,0,0, 0,1,0, E, dt, x);
}

static void
calce_const_c(struct ggcm_mhd_step *step, struct mrc_fld *E,
	      mrc_fld_data_t dt, struct mrc_fld *x,
	      struct mrc_fld *curr, struct mrc_fld *resis)
{
  bcthy3z_const(step, 0,1,2, 0,1,1, 0,1,0, 0,0,1, E, dt, x, curr, resis);
  bcthy3z_const(step, 1,2,0, 1,0,1, 0,0,1, 1,0,0, E, dt, x, curr, resis);
  bcthy3z_const(step, 2,0,1, 1,1,0, 1,0,0, 0,1,0, E, dt, x, curr, resis);
}

static void
calce_c(struct ggcm_mhd_step *step, struct mrc_fld *E,
	struct mrc_fld *x, mrc_fld_data_t dt)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;

  switch (mhd->par.magdiffu) {
  case MAGDIFFU_NL1:
    calc_resis_nl1_c(mhd);
    calce_nl1_c(step, E, dt, x);
    break;
  case MAGDIFFU_CONST: {
    struct mrc_fld *curr = sub->curr;
    struct mrc_fld *resis = sub->resis;
    calc_resis_const_c(step, curr, resis, x);
    calce_const_c(step, E, dt, x, curr, resis);
    break;
  }
  default:
    assert(0);
  }
}

static void
update_ct(struct ggcm_mhd *mhd, struct mrc_fld *x, struct mrc_fld *E_ec,
	  mrc_fld_data_t dt)
{
  float *bd3x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD3);
  float *bd3y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD3);
  float *bd3z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD3);

  mrc_fld_foreach(x, i,j,k, 0, 0) {
    F3(x, _B1X, i,j,k) -=
      dt * (bd3y[j] * (F3(E_ec, 2, i,j+1,k) - F3(E_ec, 2, i,j,k)) -
	    bd3z[k] * (F3(E_ec, 1, i,j,k+1) - F3(E_ec, 1, i,j,k)));
    F3(x, _B1Y, i,j,k) -=
      dt * (bd3z[k] * (F3(E_ec, 0, i,j,k+1) - F3(E_ec, 0, i,j,k)) -
	    bd3x[i] * (F3(E_ec, 2, i+1,j,k) - F3(E_ec, 2, i,j,k)));
    F3(x, _B1Z, i,j,k) -=
      dt * (bd3x[i] * (F3(E_ec, 1, i+1,j,k) - F3(E_ec, 1, i,j,k)) -
	    bd3y[j] * (F3(E_ec, 0, i,j+1,k) - F3(E_ec, 0, i,j,k)));
  } mrc_fld_foreach_end;
}

static void
pushstage_c(struct ggcm_mhd_step *step, mrc_fld_data_t dt,
	    struct mrc_fld *x_curr, struct mrc_fld *x_next,
	    struct mrc_fld *prim, int limit)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld **fluxes = sub->fluxes;
  struct mrc_fld *E = sub->E;
  struct mrc_fld *masks = sub->masks;

  rmaskn_c(step);

  pushfv_c(step, fluxes, dt, x_curr, x_next, limit);
  update_finite_volume(mhd, x_next, fluxes, masks, dt);
  pushpp_c(step, dt, x_next, prim);

  push_ej_c(step, dt, x_curr, x_next);
  calce_c(step, E, x_curr, dt);
  update_ct(mhd, x_next, E, dt);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_pred

static void
ggcm_mhd_step_c_pred(struct ggcm_mhd_step *step,
		     struct mrc_fld *x_half, struct mrc_fld *x)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct mrc_fld *prim = sub->prim;

  ggcm_mhd_step_c_primvar(step, prim, x);
  primbb_c2_c(step->mhd, _RR1);
  zmaskn_c(step->mhd);

  mrc_fld_data_t dt = .5f * step->mhd->dt;

  // set x_half = x^n, then advance to n+1/2
  mrc_fld_copy_range(x_half, x, 0, 8);

  pushstage_c(step, dt, x, x_half, prim, LIMIT_NONE);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_corr

static void
ggcm_mhd_step_c_corr(struct ggcm_mhd_step *step,
		     struct mrc_fld *x, struct mrc_fld *x_half)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct mrc_fld *prim = sub->prim;

  ggcm_mhd_step_c_primvar(step, prim, x_half);
  //  primbb_c2_c(step->mhd, _RR2);
  //  zmaskn_c(step->mhd);

  pushstage_c(step, step->mhd->dt, x_half, x, prim, LIMIT_1);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_run

static void
ggcm_mhd_step_c_run(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd *mhd = step->mhd;
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct mrc_fld *x_half = sub->x_half;

  float dtn;
  if (step->do_nwst) {
    newstep(mhd, &dtn);
  }

  ggcm_mhd_fill_ghosts(mhd, x, _RR1, mhd->time);
  ggcm_mhd_step_c_pred(step, x_half, x);

  ggcm_mhd_fill_ghosts(mhd, x_half, 0, mhd->time + mhd->bndt);
  ggcm_mhd_step_c_corr(step, x, x_half);

  if (step->do_nwst) {
    dtn = fminf(1., dtn); // FIXME, only kept for compatibility

    if (dtn > 1.02 * mhd->dt || dtn < mhd->dt / 1.01) {
      mpi_printf(ggcm_mhd_comm(mhd), "switched dt %g <- %g\n",
		 dtn, mhd->dt);
      mhd->dt = dtn;
    }
  }
}

// ----------------------------------------------------------------------
// subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_step_c3, x)
static struct param ggcm_mhd_step_c_descr[] = {
  { "x_half"          , VAR(x_half)          , MRC_VAR_OBJ(mrc_fld)           },
  { "prim"            , VAR(prim)            , MRC_VAR_OBJ(mrc_fld)           },
  { "tmp"             , VAR(tmp)             , MRC_VAR_OBJ(mrc_fld)           },
  { "fluxes[0]"       , VAR(fluxes[0])       , MRC_VAR_OBJ(mrc_fld)           },
  { "fluxes[1]"       , VAR(fluxes[1])       , MRC_VAR_OBJ(mrc_fld)           },
  { "fluxes[2]"       , VAR(fluxes[2])       , MRC_VAR_OBJ(mrc_fld)           },
  { "b"               , VAR(b)               , MRC_VAR_OBJ(mrc_fld)           },
  { "c"               , VAR(c)               , MRC_VAR_OBJ(mrc_fld)           },
  { "E"               , VAR(E)               , MRC_VAR_OBJ(mrc_fld)           },
  { "curr"            , VAR(curr)            , MRC_VAR_OBJ(mrc_fld)           },
  { "resis"           , VAR(resis)           , MRC_VAR_OBJ(mrc_fld)           },

  {},
};
#undef VAR

