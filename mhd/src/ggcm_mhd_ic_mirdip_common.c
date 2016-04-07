#include <stdlib.h>
#include <math.h>

#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_bndsw.h"
#include "ggcm_mhd_dipole.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_dipole_private.h" // FIXME

#include "mrc_domain.h"

// ======================================================================
// ggcm_mhd_ic subclass "mirdip"

struct ggcm_mhd_ic_mirdip {
  float xxx1;
  float xxx2;
  float xmir;
  float rrini;
  float prat; // determines initial interior pressure as a fraction of solar wind pressure
  float stretch_tail; // vx == 0 is stretched by this factor on the tail side

  float dipole_moment[3];

  // solar wind values to use if there is not "bndsw" object around
  double bnvals[SW_NR];
};

#define ggcm_mhd_ic_mirdip(ic) mrc_to_subobj(ic, struct ggcm_mhd_ic_mirdip)

// ----------------------------------------------------------------------
// get_solar_wind

static void
get_solar_wind(struct ggcm_mhd_ic *ic, float vals[])
{
  struct ggcm_mhd_ic_mirdip *sub = ggcm_mhd_ic_mirdip(ic);
  struct ggcm_mhd *mhd = ic->mhd;

  vals[SW_RR] = sub->bnvals[SW_RR] / mhd->rrnorm;
  vals[SW_VX] = sub->bnvals[SW_VX] / mhd->vvnorm;
  vals[SW_VY] = sub->bnvals[SW_VY] / mhd->vvnorm;
  vals[SW_VZ] = sub->bnvals[SW_VZ] / mhd->vvnorm;
  vals[SW_PP] = sub->bnvals[SW_PP] / mhd->ppnorm;
  vals[SW_BX] = sub->bnvals[SW_BX] / mhd->bbnorm;
  vals[SW_BY] = sub->bnvals[SW_BY] / mhd->bbnorm;
  vals[SW_BZ] = sub->bnvals[SW_BZ] / mhd->bbnorm;
}

// ----------------------------------------------------------------------
// lmbda

static inline mrc_fld_data_t
lmbda(mrc_fld_data_t XX, mrc_fld_data_t R1, mrc_fld_data_t R2)
{
  mrc_fld_data_t LL = (XX - R2) / (R1 - R2);
  LL = fmax(LL, 0.);
  LL = fmin(LL, 1.);
  return LL;
}

// ----------------------------------------------------------------------
// vxsta1

static inline mrc_fld_data_t
vxsta1(mrc_fld_data_t x, mrc_fld_data_t y, mrc_fld_data_t z, mrc_fld_data_t v0,
       mrc_fld_data_t r1, mrc_fld_data_t r2, mrc_fld_data_t str, mrc_fld_data_t xmir)
{
  mrc_fld_data_t xx = x;
  if (xx > 0.) {
    xx = x / str;
  }
  mrc_fld_data_t r = sqrt(sqr(xx) + sqr(y) + sqr(z));
  mrc_fld_data_t s = (r - r1) / (r2 - r1);
  s = fmax(s, 0.);
  s = fmin(s, 1.);
  if (xmir != 0. && x < xmir) {
    return v0;
  } else {
    return v0 * s;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_ini1

static void
ggcm_mhd_ic_mirdip_ini1(struct ggcm_mhd_ic *ic, float vals[])
{
  struct ggcm_mhd_ic_mirdip *sub = ggcm_mhd_ic_mirdip(ic);
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);

  mrc_fld_data_t xxx1 = sub->xxx1, xxx2 = sub->xxx2, xmir = sub->xmir;
  mrc_fld_data_t rrini = sub->rrini, stretch_tail = sub->stretch_tail;

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 2, 2) {
      mrc_fld_data_t tmplam = lmbda(MRC_MCRDX(crds, ix, p), -xxx1, -xxx2);
      RR_(fld, ix,iy,iz, p) = 
	tmplam * vals[SW_RR] + (1. - tmplam) * (rrini / mhd->rrnorm);
      
      VX_(fld, ix,iy,iz, p) = vxsta1(MRC_MCRDX(crds, ix, p), MRC_MCRDY(crds, iy, p), MRC_MCRDZ(crds, iz, p),
				     vals[SW_VX], xxx2, xxx1, stretch_tail, xmir);
      VY_(fld, ix,iy,iz, p) = 0.;
      VZ_(fld, ix,iy,iz, p) = 0.;
      
      mrc_fld_data_t ppmin = sub->prat * vals[SW_PP];
      PP_(fld, ix,iy,iz, p) = vals[SW_PP] * tmplam + (1. - tmplam) * ppmin;
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(fld, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_get_mhd_dipole

static struct ggcm_mhd_dipole *
ggcm_mhd_ic_mirdip_get_mhd_dipole(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_dipole *mhd_dipole = ggcm_mhd_get_var_obj(ic->mhd, "mhd_dipole");
  if (mhd_dipole) {
    ggcm_mhd_dipole_get(mhd_dipole);
  } else {
    mhd_dipole = ggcm_mhd_dipole_create(ggcm_mhd_ic_comm(ic));
    ggcm_mhd_dipole_set_type(mhd_dipole, FLD_TYPE);
    ggcm_mhd_dipole_set_from_options(mhd_dipole);
    ggcm_mhd_dipole_set_param_obj(mhd_dipole, "mhd", ic->mhd);
    ggcm_mhd_dipole_setup(mhd_dipole);
  }

  return mhd_dipole;
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_ini_b

static void
ggcm_mhd_ic_mirdip_ini_b(struct ggcm_mhd_ic *ic, float b_sw[3])
{
  struct ggcm_mhd_ic_mirdip *sub = ggcm_mhd_ic_mirdip(ic);
  struct ggcm_mhd *mhd = ic->mhd;
  struct ggcm_mhd_dipole *mhd_dipole = ggcm_mhd_ic_mirdip_get_mhd_dipole(ic);

  float x0[3] = {0.0, 0.0, 0.0};

  mprintf("mirdip_ini_b dipole moment: %f %f %f\n", sub->dipole_moment[0],
  	  sub->dipole_moment[1], sub->dipole_moment[2]);

  struct mrc_fld *b_base = mrc_fld_make_view(mhd->fld, BX, BX + 3);
  ggcm_mhd_dipole_add_dipole(mhd_dipole, b_base, x0, sub->dipole_moment, sub->xmir, 0.);

  if (sub->xmir != 0.0) {
    x0[0] = 2.0 * sub->xmir;
    sub->dipole_moment[0] *= -1.0;
    ggcm_mhd_dipole_add_dipole(mhd_dipole, b_base, x0, sub->dipole_moment, sub->xmir, 1.);
    sub->dipole_moment[0] *= -1.0;
  }

  mrc_fld_destroy(b_base);

  struct mrc_fld *f = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  struct mrc_fld *b0 = mrc_fld_get_as(mhd->b0, FLD_TYPE);

  // finish up
  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    mrc_fld_foreach(f, ix,iy,iz, 1, 1) {
      for (int d = 0; d < 3; d++){
	// add B_IMF
	M3(f, BX + d, ix,iy,iz, p) += b_sw[d];

	// subtract previously calculated background dipole
	M3(f, BX + d, ix,iy,iz, p) -= M3(b0, d, ix,iy,iz, p);
      }
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(f, mhd->fld);
  mrc_fld_put_as(b0, mhd->b0);

  ggcm_mhd_dipole_put(mhd_dipole);
}

// ----------------------------------------------------------------------
// ggcm_mhd_mirdip_ic_run

static void
ggcm_mhd_ic_mirdip_run(struct ggcm_mhd_ic *ic)
{
  // ini_b is done when we still have the primitive variables stored in
  // mhd->fld, so it needs to use the right ("regular") staggering because
  // the convert_from_primitive() will fix up the staggering to whatever
  // we really need in the end.
  int mhd_type_save;
  mrc_fld_get_param_int(ic->mhd->fld, "mhd_type", &mhd_type_save);
  int mhd_type = MT_PRIMITIVE;
  if (mhd_type_save == MT_FULLY_CONSERVATIVE_CC) {
    mhd_type = MT_PRIMITIVE_CC;
  }
  mrc_fld_set_param_int(ic->mhd->fld, "mhd_type", mhd_type);

  float vals[SW_NR];
  get_solar_wind(ic, vals);

  ggcm_mhd_ic_mirdip_ini1(ic, vals);
  ggcm_mhd_ic_mirdip_ini_b(ic, &vals[SW_BX]);

  mrc_fld_set_param_int(ic->mhd->fld, "mhd_type", mhd_type_save);
  if (mhd_type == MT_PRIMITIVE) {
    ggcm_mhd_convert_from_primitive(ic->mhd, ic->mhd->fld);
  } else if (mhd_type == MT_PRIMITIVE_CC) {
    ggcm_mhd_convert_from_primitive_cc(ic->mhd, ic->mhd->fld);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_vector_potential_b0

static double
ggcm_mhd_ic_mirdip_vector_potential_b0(struct ggcm_mhd_ic *ic, int m, double x[3])
{
  struct ggcm_mhd_ic_mirdip *sub = ggcm_mhd_ic_mirdip(ic);
  static struct ggcm_mhd_dipole *mhd_dipole;
  static bool first_time = true;
  if (first_time) {
    mhd_dipole = ggcm_mhd_ic_mirdip_get_mhd_dipole(ic);
    first_time = false;
  }

  float x0[3] = { 0.f, 0.f, 0.f };
  float *moment = sub->dipole_moment;
  float xmir = 0.f;

  // get main dipole vector potential
  double A = ggcm_mhd_dipole_vector_potential(mhd_dipole, m, x, x0, moment, xmir);

  return A;
}

// ----------------------------------------------------------------------
// ggcm_mhd_mirdip_ic_init_b0

static void
ggcm_mhd_ic_mirdip_init_b0(struct ggcm_mhd_ic *ic, struct mrc_fld *b0_base)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct ggcm_mhd_dipole *mhd_dipole = ggcm_mhd_ic_mirdip_get_mhd_dipole(ic);

  float vals[SW_NR];
  get_solar_wind(ic, vals);

  int mhd_type;
  mrc_fld_get_param_int(mhd->fld, "mhd_type", &mhd_type);

  struct mrc_fld *a_base = ggcm_mhd_get_3d_fld(mhd, 3);
  struct mrc_fld *a = mrc_fld_get_as(a_base, FLD_TYPE);
  struct mrc_fld *b0 = mrc_fld_get_as(b0_base, FLD_TYPE);

  // calculate A first, then take its curl
  for (int p = 0; p < mrc_fld_nr_patches(a); p++) {
    mrc_fld_foreach(a, ix,iy,iz, 2, 2) {
      for (int m = 0; m < 3; m++) {
	float crd[3];
	if (mhd_type == MT_PRIMITIVE_CC ||
	    mhd_type == MT_FULLY_CONSERVATIVE_CC) { // cell-centered B
	  ggcm_mhd_get_crds_cc(mhd, ix,iy,iz, p, crd);
	} else {
	  ggcm_mhd_get_crds_ec(mhd, ix,iy,iz, p, m, crd);
	}
	double x[3] = { crd[0], crd[1], crd[2] };
	M3(a, m, ix,iy,iz, p) = ggcm_mhd_ic_mirdip_vector_potential_b0(ic, m, x);
      }
    } mrc_fld_foreach_end;
  }

  // B = keep * B + curl A

  mrc_fld_data_t curl_a[3];

  if (mhd_type == MT_PRIMITIVE_CC ||
      mhd_type == MT_FULLY_CONSERVATIVE_CC) { // cell-centered B
    for (int p = 0; p < mrc_fld_nr_patches(b0); p++) {
      float *fd1x = ggcm_mhd_crds_get_crd_p(mhd->crds, 0, FD1, p);
      float *fd1y = ggcm_mhd_crds_get_crd_p(mhd->crds, 1, FD1, p);
      float *fd1z = ggcm_mhd_crds_get_crd_p(mhd->crds, 2, FD1, p);
      
      mrc_fld_foreach(b0, ix,iy,iz, 1, 1) {
	curl_a[0] = ((M3(a, 2, ix,iy+1,iz, p) - M3(a, 2, ix,iy-1,iz, p)) * .5f * fd1y[iy] -
		     (M3(a, 1, ix,iy,iz+1, p) - M3(a, 1, ix,iy,iz-1, p)) * .5f * fd1z[iz]);
	curl_a[1] = ((M3(a, 0, ix,iy,iz+1, p) - M3(a, 0, ix,iy,iz-1, p)) * .5f * fd1z[iz] -
		     (M3(a, 2, ix+1,iy,iz, p) - M3(a, 2, ix-1,iy,iz, p)) * .5f * fd1x[ix]);
	curl_a[2] = ((M3(a, 1, ix+1,iy,iz, p) - M3(a, 1, ix-1,iy,iz, p)) * .5f * fd1x[ix] -
		     (M3(a, 0, ix,iy+1,iz, p) - M3(a, 0, ix,iy-1,iz, p)) * .5f * fd1y[iy]);
	
	float crd_cc[3];
	ggcm_mhd_get_crds_cc(mhd, ix,iy,iz, p, crd_cc);
	float r = sqrtf(sqr(crd_cc[0]) + sqr(crd_cc[1]) + sqr(crd_cc[2]));
	// only set B outside of r1lim
	if (r >= mhd_dipole->r1lim) {
	  for (int d = 0; d < 3; d++){
	    M3(b0, d, ix,iy,iz, p) = curl_a[d];
	  }
	}
      } mrc_fld_foreach_end;
    }
  } else { // face-centered B
    // FIXME, this doesn't fill Bnormal one ghost cell out on the right (high) side
    for (int p = 0; p < mrc_fld_nr_patches(b0); p++) {
      float *bd3x = ggcm_mhd_crds_get_crd_p(mhd->crds, 0, BD3, p);
      float *bd3y = ggcm_mhd_crds_get_crd_p(mhd->crds, 1, BD3, p);
      float *bd3z = ggcm_mhd_crds_get_crd_p(mhd->crds, 2, BD3, p);
      
      mrc_fld_foreach(b0, ix,iy,iz, 2, 1) {
	curl_a[0] = ((M3(a, 2, ix,iy+1,iz, p) - M3(a, 2, ix,iy,iz, p)) * bd3y[iy] -
		     (M3(a, 1, ix,iy,iz+1, p) - M3(a, 1, ix,iy,iz, p)) * bd3z[iz]);
	curl_a[1] = ((M3(a, 0, ix,iy,iz+1, p) - M3(a, 0, ix,iy,iz, p)) * bd3z[iz] -
		     (M3(a, 2, ix+1,iy,iz, p) - M3(a, 2, ix,iy,iz, p)) * bd3x[ix]);
	curl_a[2] = ((M3(a, 1, ix+1,iy,iz, p) - M3(a, 1, ix,iy,iz, p)) * bd3x[ix] -
		     (M3(a, 0, ix,iy+1,iz, p) - M3(a, 0, ix,iy,iz, p)) * bd3y[iy]);
	
	switch (mhd_type) {
	case MT_SEMI_CONSERVATIVE_GGCM:
	  M3(b0, 0, ix-1,iy,iz, p) = curl_a[0];
	  M3(b0, 1, ix,iy-1,iz, p) = curl_a[1];
	  M3(b0, 2, ix,iy,iz-1, p) = curl_a[2];
	  break;
	case MT_PRIMITIVE:
	case MT_SEMI_CONSERVATIVE:
	case MT_FULLY_CONSERVATIVE:
	  for (int d = 0; d < 3; d++){
	    float crd_fc[3];
	    ggcm_mhd_get_crds_fc(mhd, ix,iy,iz, p, d, crd_fc);
	    
	    // only set B outside of r1lim
	    float r = sqrtf(sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]));
	    if (r >= mhd_dipole->r1lim) {
	      M3(b0, d, ix,iy,iz, p) = curl_a[d];
	    }
	  }
	  break;
	default:
	  assert(0);
	}
      } mrc_fld_foreach_end;
    }
  }

  mrc_fld_put_as(a, a_base);
  ggcm_mhd_put_3d_fld(mhd, a_base);

  //  ggcm_mhd_dipole_add_dipole(mhd_dipole, b0, x0, sub->dipole_moment, 0., 0.);

  // finish up
  for (int p = 0; p < mrc_fld_nr_patches(b0); p++) {
    mrc_fld_foreach(b0, ix,iy,iz, 2, 2) {
      for (int d = 0; d < 3; d++){
	// add B_IMF
	M3(b0, d, ix,iy,iz, p) += vals[SW_BX + d];
      }
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(b0, b0_base);

  ggcm_mhd_dipole_put(mhd_dipole);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mirdip_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_mirdip, x)
static struct param ggcm_mhd_ic_mirdip_descr[] = {
  { "xxx1"         , VAR(xxx1)         , PARAM_FLOAT(14.0)         },
  { "xxx2"         , VAR(xxx2)         , PARAM_FLOAT(12.0)         },
  { "xmir"         , VAR(xmir)         , PARAM_FLOAT(-15.0)        },  // off if == 0.0
  { "rrini"        , VAR(rrini)        , PARAM_FLOAT(3.0)          },
  { "prat"         , VAR(prat)         , PARAM_FLOAT(.5)           },
  { "stretch_tail" , VAR(stretch_tail) , PARAM_FLOAT(2.)           },

  { "dipole_moment", VAR(dipole_moment), PARAM_FLOAT3(0., 0., -1.) },

  { "rr"           , VAR(bnvals[SW_RR]), PARAM_DOUBLE(1.)          },
  { "pp"           , VAR(bnvals[SW_PP]), PARAM_DOUBLE(1.)          },
  { "vx"           , VAR(bnvals[SW_VX]), PARAM_DOUBLE(0.)          },
  { "vy"           , VAR(bnvals[SW_VY]), PARAM_DOUBLE(0.)          },
  { "vz"           , VAR(bnvals[SW_VZ]), PARAM_DOUBLE(0.)          },
  { "bx"           , VAR(bnvals[SW_BX]), PARAM_DOUBLE(0.)          },
  { "by"           , VAR(bnvals[SW_BY]), PARAM_DOUBLE(0.)          },
  { "bz"           , VAR(bnvals[SW_BZ]), PARAM_DOUBLE(0.)          },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic subclass "mirdip"

struct ggcm_mhd_ic_ops ggcm_mhd_ic_mirdip_ops = {
  .name             = ggcm_mhd_ic_mirdip_name,
  .size             = sizeof(struct ggcm_mhd_ic_mirdip),
  .param_descr      = ggcm_mhd_ic_mirdip_descr,
  .run              = ggcm_mhd_ic_mirdip_run,
  .init_b0          = ggcm_mhd_ic_mirdip_init_b0,
};
