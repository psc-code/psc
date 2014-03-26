
#include "mrc_crds_gen_private.h"

#include <mrc_domain.h>
#include <mrc_params.h>
#include <math.h>
#include <assert.h>

#define ITERMAX 2000

// ======================================================================
// mrc_crds_gen subclass "ggcm_x"

struct mrc_crds_gen_ggcm_x {
  double x1;
  double x3;
  double dmm;
  double x5;
  double h0;
  double hn;
  double hm;
  double hmm;
  double b1;
  double b2;
  double b3;

  // these aren't parameters, but also needed in the context
  double fak;
  double dxi;
  int n;
  double x0;
  double xn;
};

#define mrc_crds_gen_ggcm_x(gen) mrc_to_subobj(gen, struct mrc_crds_gen_ggcm_x)

static double
tanhp(double x)
{
  return .5 * (1. + tanh(x));
}

static double
tanhm(double x)
{
  return .5 * (1. - tanh(x));
}

//
// h0  |  hm  |  hmm  |  hm  |  hn
//     x1   x3-dmm  x3+dmm   x5
//     b1     b3      b3     b2

static double 
f(double x, struct mrc_crds_gen_ggcm_x *par)
{
  return par->fak * (par->hm +
		     (par->h0 - par->hm) * tanhm(par->b1 * (x - par->x1)) +
		     (par->hn - par->hm) * tanhp(par->b2 * (x - par->x5)) +
		     (par->hmm- par->hm) * (1. - tanhm(par->b3 * (x - (par->x3 - par->dmm)))
 					       - tanhp(par->b3 * (x - (par->x3 + par->dmm)))));
}

static double
ff(double xi, struct mrc_crds_gen_ggcm_x *par)
{
  int ns = fabsf((xi - 1.) / par->dxi + .5) + 5;
  double ddxi = (xi - 1.) / ns;

  double x = par->x0;
  for (int i = 0; i < ns; i++) {
    // 4th order runge kutta
    double xk1 = ddxi * f(x, par);
    double xk2 = ddxi * f(x + .5 * xk1, par);
    double xk3 = ddxi * f(x + .5 * xk2, par);
    double xk4 = ddxi * f(x + .5 * xk3, par);
    x += (1./6.) * (xk1 + xk4) + (1./3.) * (xk2 + xk3);
  }

  return x;
}

// ----------------------------------------------------------------------
// mrc_crds_gen_ggcm_x_run

static void
mrc_crds_gen_ggcm_x_run(struct mrc_crds_gen *gen, double *xx, double *dx)
{
  struct mrc_crds_gen_ggcm_x *sub = mrc_crds_gen_ggcm_x(gen);

  sub->n = gen->dims[gen->d];
  
  float xl[3], xh[3];
  mrc_crds_get_param_float3(gen->crds, "l", xl);
  mrc_crds_get_param_float3(gen->crds, "h", xh);
  int sw;
  mrc_crds_get_param_int(gen->crds, "sw", &sw);
  
  sub->x0  = xl[gen->d];
  sub->xn  = xh[gen->d];
  //  printf("gridx: n = %d\n", sub->n);
  sub->dxi = .1;
  double xn1, xn1o = 0.;
  sub->fak = 1.;
  for (int k = 0; k < ITERMAX; k++) {
    xn1 = ff(sub->n, sub);
    double rel = 200. * fabsf(xn1 - xn1o) / (fabsf(xn1) + fabs(xn1o));
    //    printf("%d: convergence test: %g, %g, %g\n", k, dxi, xn1, rel);
    if (rel < .02)
      break;
    xn1o = xn1;
    sub->dxi *= .95;
  }
  //  printf("gridx: convergence test: %g, %g\n", sub->dxi, xn1);

  sub->fak = .1;
  double s = .1;
  for (int k = 0; k < ITERMAX; k++) {
    xn1 = ff(sub->n, sub);
    //    printf("%d: fak=%g xn1=%g xn=%g\n", k, sub->fak, xn1, sub->xn);

    if (xn1 > sub->xn) {
      if (fabsf(xn1 - sub->xn) < 1e-2)
	break;
      sub->fak -= s;
      s *= .1;
    } else {
      sub->fak += s;
    }
  }
  //  printf("gridx: fak=%g xn1=%g xn=%g\n", sub->fak, xn1, sub->xn);

  int step_out = sub->n / 20;
  if (step_out > 20) step_out = 20;

  for (int i = -sw; i <= sub->n + sw; i++) {
    xx[i] = ff(i+1, sub);
    dx[i] = f(xx[i], sub);

#if 0
    if (i % step_out == 1) {
      printf("%d: xx=%g dx=%g\n", i, xx[i], dx[i]);
    }
#endif
  }
}

// ----------------------------------------------------------------------
// mrc_crds_gen "ggcm_yz" description

#define VAR(x) (void *)offsetof(struct mrc_crds_gen_ggcm_x, x)
static struct param mrc_crds_gen_ggcm_x_descr[] = {
  { "x1"              , VAR(x1)              , PARAM_DOUBLE(-26.)     },
  { "x3"              , VAR(x3)              , PARAM_DOUBLE(10.)      },
  { "dmm"             , VAR(dmm)             , PARAM_DOUBLE(8.)       },
  { "x5"              , VAR(x5)              , PARAM_DOUBLE(80.)      },
  { "h0"              , VAR(h0)              , PARAM_DOUBLE(2.)       },
  { "hn"              , VAR(hn)              , PARAM_DOUBLE(43.)      },
  { "hm"              , VAR(hm)              , PARAM_DOUBLE(1.7)      },
  { "hmm"             , VAR(hmm)             , PARAM_DOUBLE(.7)       },
  { "b1"              , VAR(b1)              , PARAM_DOUBLE(.15)      },
  { "b2"              , VAR(b2)              , PARAM_DOUBLE(.025)     },
  { "b3"              , VAR(b3)              , PARAM_DOUBLE(.3)       },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mrc_crds_gen subclass "ggcm_x"

struct mrc_crds_gen_ops mrc_crds_gen_ggcm_x_ops = {
  .name        = "ggcm_x",
  .size        = sizeof(struct mrc_crds_gen_ggcm_x),
  .param_descr = mrc_crds_gen_ggcm_x_descr,
  .run         = mrc_crds_gen_ggcm_x_run,
};


