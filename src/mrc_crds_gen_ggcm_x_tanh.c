
#include "mrc_ggcm_gridx_gen.h"

#include <mrc_params.h>
#include <math.h>
#include <assert.h>

// ======================================================================
// mrc_crds_gen subclass "ggcm_x_tanh"

struct mrc_crds_gen_ggcm_x_tanh {
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
};

#define mrc_crds_gen_ggcm_x_tanh(gen) mrc_to_subobj(gen, struct mrc_crds_gen_ggcm_x_tanh)

static double
tanhp(double x)
{
  return 0.5 * (1.0 + tanh(x));
}

static double
tanhm(double x)
{
  return 0.5 * (1.0 - tanh(x));
}

static double
f_tanh(struct mrc_crds_gen *gen, double x, double fak)
{
  struct mrc_crds_gen_ggcm_x_tanh *par = mrc_crds_gen_ggcm_x_tanh(gen);

  return fak * (par->hm +
         (par->h0 - par->hm) * tanhm(par->b1 * (x - par->x1)) +
         (par->hn - par->hm) * tanhp(par->b2 * (x - par->x5)) +
         (par->hmm - par->hm) * (1.0 - tanhm(par->b3 * (x - (par->x3 - par->dmm)))
                 - tanhp(par->b3 * (x - (par->x3 + par->dmm)))));
}

static void
mrc_crds_gen_ggcm_x_tanh_run(struct mrc_crds_gen *gen, double *xx, double *dx)
{
  generate_ggcm_x_grid(gen, xx, dx, f_tanh);
}

// ----------------------------------------------------------------------
// mrc_crds_gen "ggcm_x_tanh" description

#define VAR(x) (void *)offsetof(struct mrc_crds_gen_ggcm_x_tanh, x)
static struct param mrc_crds_gen_ggcm_x_tanh_descr[] = {
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
// mrc_crds_gen subclass "ggcm_x_tanh"

struct mrc_crds_gen_ops mrc_crds_gen_ggcm_x_tanh_ops = {
  .name        = "ggcm_x_tanh",
  .size        = sizeof(struct mrc_crds_gen_ggcm_x_tanh),
  .param_descr = mrc_crds_gen_ggcm_x_tanh_descr,
  .run         = mrc_crds_gen_ggcm_x_tanh_run,
};
