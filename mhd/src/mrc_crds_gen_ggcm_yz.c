
#include "mrc_crds_gen_private.h"

#include <mrc_params.h>
#include <mrc_domain.h>
#include <math.h>
#include <assert.h>

// ======================================================================
// mrc_crds_gen subclass "ggcm_yz"

struct mrc_crds_gen_ggcm_yz {
  float dx0;
  float xn;
  float xm;
  float xshift;
};

#define mrc_crds_gen_ggcm_yz(gen) mrc_to_subobj(gen, struct mrc_crds_gen_ggcm_yz)

// ----------------------------------------------------------------------
// acoff

static float
acoff(int n, float y, float xm, float xn, float d0)
{
  float x = n - .5;
  float yy = y;
  yy /= d0 * x;
  yy = pow(yy, 1./xm);
  yy -= 1.;
  yy /= pow(x, 2.*xn);
  return yy;
}

// ----------------------------------------------------------------------
// mrc_crds_gen_ggcm_yz_run

static void
mrc_crds_gen_ggcm_yz_run(struct mrc_crds_gen *gen, float *xx, float *dx)
{
  struct mrc_crds_gen_ggcm_yz *sub = mrc_crds_gen_ggcm_yz(gen);

  int gdims[3];
  mrc_domain_get_global_dims(gen->crds->domain, gdims);
  int n = gdims[gen->d];
  
  float xl[3], xh[3];
  mrc_crds_get_param_float3(gen->crds, "l", xl);
  mrc_crds_get_param_float3(gen->crds, "h", xh);
  
  // FIXME, maybe we should calculate xx2, xshift from this...
  assert(xl[gen->d] == -xh[gen->d]);
  float xx2 = xh[gen->d];

  int nx2 = n / 2;
  int nx1 = 1 - n / 2;
  float a = acoff(nx2, xx2, sub->xm, sub->xn, sub->dx0);
  //  printf("gridyz: n = %d nx12 = %d, %d a = %g\n", gen->n, nx1, nx2, a);

  for (int i = -2; i < n + 2; i++) {
    float x = i + nx1 - .5;
    float d0 = sub->dx0;
    float xn = sub->xn;
    float xm = sub->xm;
    float s = 1 + a*(pow(x, (2.*xn)));
    float sm = pow(s, xm);
    float dg = d0 * (sm + xm*x*2.*xn*a*(pow(x, (2.*xn-1.))) * sm / s);
    float g = d0 * x * sm - sub->xshift;
    xx[i] = g;
    dx[i] = dg;
  }
}

// ----------------------------------------------------------------------
// mrc_crds_gen "ggcm_yz" description

#define VAR(x) (void *)offsetof(struct mrc_crds_gen_ggcm_yz, x)
static struct param mrc_crds_gen_ggcm_yz_descr[] = {
  { "center_spacing"  , VAR(dx0)             , PARAM_FLOAT(.4)       },
  { "center_shift"    , VAR(xshift)          , PARAM_FLOAT(0.)       },
  { "xn"              , VAR(xn)              , PARAM_FLOAT(2.)       },
  { "xm"              , VAR(xm)              , PARAM_FLOAT(.5)       },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mrc_crds_gen subclass "ggcm_yz"

struct mrc_crds_gen_ops mrc_crds_gen_ggcm_yz_ops = {
  .name        = "ggcm_yz",
  .size        = sizeof(struct mrc_crds_gen_ggcm_yz),
  .param_descr = mrc_crds_gen_ggcm_yz_descr,
  .run         = mrc_crds_gen_ggcm_yz_run,
};


