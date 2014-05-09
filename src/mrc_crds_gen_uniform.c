
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
  double h = (gen->xh - gen->xl) / gen->n;

  for (int i = -gen->sw; i < gen->n + gen->sw; i++) {
    xx[i] = gen->xl + (i + .5) * h;
    dx[i] = h;
  }
}

// ----------------------------------------------------------------------
// mrc_crds_gen subclass "uniform"

struct mrc_crds_gen_ops mrc_crds_gen_uniform_ops = {
  .name        = "uniform",
  .run         = mrc_crds_gen_uniform_run,
};


