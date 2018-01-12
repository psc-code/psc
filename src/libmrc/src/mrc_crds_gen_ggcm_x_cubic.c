
#include "mrc_ggcm_gridx_gen.h"

#include <mrc_params.h>
#include <mrc_crds.h>

#include <math.h>
#include <assert.h>

// ======================================================================
// mrc_crds_gen subclass "ggcm_x_cubic"

struct mrc_crds_gen_ggcm_x_cubic {
  double w0;  // overall constant offset?
  double w1;  // weight of 1st fsm function
  double a1;  // a parameter for first fsm
  double b1;  // b parameter for first fsm
  double w2;
  double a2;
  double b2;
};

#define mrc_crds_gen_ggcm_x_cubic(gen) mrc_to_subobj(gen, struct mrc_crds_gen_ggcm_x_cubic)

static double
fsm3(double xx, double a, double b)
{
  double x, y;

  x = xx - 0.5 * (b + a);
  x = 2.0 * x / (b - a);

  if (x > 0.0) {
    y = 1.0 + pow(x - 1.0, 3.0);
  } else {
    y = -1.0 + pow(x + 1.0, 3.0);
  }

  y = fmax(-1.0, fmin(1.0, y));
  return 0.5 * (1.0 + y);
}

static double
f_fsm(struct mrc_crds_gen *gen, double x, double fak)
{
  struct mrc_crds_gen_ggcm_x_cubic *sub = mrc_crds_gen_ggcm_x_cubic(gen);
  double xnorm = gen->crds->xnorm;
  double w0 = sub->w0 / xnorm;
  double w1 = sub->w1 / xnorm, a1 = sub->a1 / xnorm, b1 = sub->b1 / xnorm;
  double w2 = sub->w2 / xnorm, a2 = sub->a2 / xnorm, b2 = sub->b2 / xnorm;

  double dx = w0 + w1 * fsm3(x, a1, b1) + w2 * fsm3(x, a2, b2);
  return fak * dx;
}

static void
mrc_crds_gen_ggcm_x_cubic_run(struct mrc_crds_gen *gen, double *xx, double *dx)
{
  generate_ggcm_x_grid(gen, xx, dx, f_fsm);
}

// ----------------------------------------------------------------------
// mrc_crds_gen "ggcm_x_cubic" description

#define VAR(x) (void *)offsetof(struct mrc_crds_gen_ggcm_x_cubic, x)
static struct param mrc_crds_gen_ggcm_x_cubic_descr[] = {
  { "w0", VAR(w0), PARAM_DOUBLE(1.0)   },
  { "w1", VAR(w1), PARAM_DOUBLE(150.0) },
  { "a1", VAR(a1), PARAM_DOUBLE(4.0)   },
  { "b1", VAR(b1), PARAM_DOUBLE(400.0) },
  { "w2", VAR(w2), PARAM_DOUBLE(2.0)   },
  { "a2", VAR(a2), PARAM_DOUBLE(-8.0)  },
  { "b2", VAR(b2), PARAM_DOUBLE(-30.0) },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mrc_crds_gen subclass "ggcm_x_cubic"

struct mrc_crds_gen_ops mrc_crds_gen_ggcm_x_cubic_ops = {
  .name        = "ggcm_x_cubic",
  .size        = sizeof(struct mrc_crds_gen_ggcm_x_cubic),
  .param_descr = mrc_crds_gen_ggcm_x_cubic_descr,
  .run         = mrc_crds_gen_ggcm_x_cubic_run,
};
