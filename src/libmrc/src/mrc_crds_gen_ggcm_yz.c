
#include "mrc_crds_gen_private.h"

#include <mrc_params.h>
#include <mrc_domain.h>
#include <math.h>
#include <assert.h>

// ======================================================================
// mrc_crds_gen subclass "ggcm_yz"

struct mrc_crds_gen_ggcm_yz {
  double dx0;
  double xn;
  double xm;
  double xshift;
};

#define mrc_crds_gen_ggcm_yz(gen) mrc_to_subobj(gen, struct mrc_crds_gen_ggcm_yz)

// ----------------------------------------------------------------------
// acoff

static double
acoff(int n, double y, double xm, double xn, double d0)
{
  double x = n - .5;
  double yy = y;
  yy /= d0 * x;
  yy = pow(yy, 1./xm);
  yy -= 1.;
  yy /= pow(x, 2.*xn);
  return yy;
}

// ----------------------------------------------------------------------
// mrc_crds_gen_ggcm_yz_run

static void
mrc_crds_gen_ggcm_yz_run(struct mrc_crds_gen *gen, double *xx, double *dx)
{
  struct mrc_crds_gen_ggcm_yz *sub = mrc_crds_gen_ggcm_yz(gen);
  double xshift_code = sub->xshift / gen->crds->xnorm;
  double dx0_code = sub->dx0 / gen->crds->xnorm;

  // FIXME, maybe we should calculate xx2, xshift from this...
  assert(gen->xl == -gen->xh);
  double xx2 = gen->xh;

  int nx2 = gen->n / 2;
  int nx1 = 1 - gen->n / 2;
  double a = acoff(nx2, xx2, sub->xm, sub->xn, dx0_code);
  //  printf("gridyz: n = %d nx12 = %d, %d a = %g\n", gen->n, nx1, nx2, a);

  for (int i = -gen->sw; i < gen->n + gen->sw; i++) {
    double x = i + nx1 - .5;
    double xn = sub->xn;
    double xm = sub->xm;
    double s = 1 + a*(pow(x, (2.*xn)));
    double sm = pow(s, xm);
    double dg = dx0_code * (sm + xm*x*2.*xn*a*(pow(x, (2.*xn-1.))) * sm / s);
    double g = dx0_code * x * sm - xshift_code;
    xx[i] = g;
    dx[i] = dg;
  }
}

// ----------------------------------------------------------------------
// mrc_crds_gen "ggcm_yz" description

#define VAR(x) (void *)offsetof(struct mrc_crds_gen_ggcm_yz, x)
static struct param mrc_crds_gen_ggcm_yz_descr[] = {
  { "center_spacing"  , VAR(dx0)             , PARAM_DOUBLE(.4)       },
  { "center_shift"    , VAR(xshift)          , PARAM_DOUBLE(0.)       },
  { "xn"              , VAR(xn)              , PARAM_DOUBLE(2.)       },
  { "xm"              , VAR(xm)              , PARAM_DOUBLE(.5)       },
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


