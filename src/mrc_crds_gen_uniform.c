
#include "mrc_crds_gen_private.h"

#include <mrc_params.h>
#include <mrc_domain.h>
#include <math.h>
#include <assert.h>

// ======================================================================
// mrc_crds_gen subclass "uniform"

// ----------------------------------------------------------------------
// mrc_crds_gen_uniform_run

static void
mrc_crds_gen_uniform_run(struct mrc_crds_gen *gen, double *xx, double *dx)
{
  int d = gen->d;
  int n = gen->dims[d];

  float xl[3], xh[3];
  mrc_crds_get_param_float3(gen->crds, "l", xl);
  mrc_crds_get_param_float3(gen->crds, "h", xh);
  int sw;
  mrc_crds_get_param_int(gen->crds, "sw", &sw);

  double h = (xh[d] - xl[d]) / n;
  
  for (int i = -sw; i < n + sw; i++) {
    xx[i] = xl[d] + (i + .5) * h;
    dx[i] = h;
  }
}

// ----------------------------------------------------------------------
// mrc_crds_gen subclass "uniform"

struct mrc_crds_gen_ops mrc_crds_gen_uniform_ops = {
  .name        = "uniform",
  .run         = mrc_crds_gen_uniform_run,
};


