
#include "mrc_crds_gen_private.h"

#include <mrc_crds.h>
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_io.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

//#include "libmrc_extensions.h"

#define sqr(x) (x)*(x)

// ======================================================================
// mrc_crd_gen subclass "gaussian"
//
// Nonuniform coordinates initiated via a Gaussian function,
// rather than the funky sort used by jimmy_crd
// This is similar to the default version in the original 
// libmrc (ie, the one found in mrc-v3)


struct mrc_crds_gen_gaussian {
  double gc_x0; // Center (ie point of smallest dx)
  double gc_r; // approximate ratio from peak to min crd
  double gc_w; // Width, i.e. sigma of the gaussian
  bool bc_cyl; // Hack for cylindrical coord bc
};

#define mrc_crds_gen_gaussian(gen) mrc_to_subobj(gen, struct mrc_crds_gen_gaussian)

static double
f_dx(struct mrc_crds_gen_gaussian *mbc, double x)
{
  double r = mbc->gc_r, w = mbc->gc_w, x0 = mbc->gc_x0;
  return r - (r - 1.) * exp(-sqr(x-x0)/(2.*sqr(w)));
}


static void
mrc_crds_gen_gaussian_run(struct mrc_crds_gen *gen, double *xx, double *dx)
{
  struct mrc_crds_gen_gaussian *sub = mrc_crds_gen_gaussian(gen);

  // Calculate a set of grid points
  double *_xx, *ncxx;
  _xx = (double *)calloc(gen->n+1+2*gen->sw, sizeof(*xx));
  ncxx = _xx + gen->sw;
  for (int jx = 0; jx < gen->n + gen->sw; jx++) {
    ncxx[jx+1] = ncxx[jx] + f_dx(sub, (jx + .5)/gen->n);
  }

  if (sub->bc_cyl) {
    for (int jx = -gen->sw; jx < 0; jx++) {
    ncxx[jx] = -ncxx[-jx];
    }
  } else { 
    for (int jx = 0; jx > -gen->sw; jx--) {
      ncxx[jx-1] = ncxx[jx] - f_dx(sub, (jx - .5)/gen->n);
    }
  }
  double fac = 1. / ncxx[gen->n];
  for (int jx = -gen->sw; jx<= gen->n + gen->sw; jx++){
    ncxx[jx] *= fac;
  }

  // average them down to the actual coordinates (basically NC to CC)
  for (int ii = -gen->sw; ii < gen->n + gen->sw; ii++) {
    xx[ii] = gen->xl + 0.5*(ncxx[ii] + ncxx[ii+1]) * (gen->xh - gen->xl);
  }
  free(_xx);
}


#define VAR(x) (void *)offsetof(struct mrc_crds_gen_gaussian, x)
static struct param mrc_crds_gen_gaussian_param_descr[] = {
  { "gc_x0"             , VAR(gc_x0)            , PARAM_DOUBLE(0.0)       },
  { "gc_r"              , VAR(gc_r)             , PARAM_DOUBLE(1.0)       },
  { "gc_w"              , VAR(gc_w)             , PARAM_DOUBLE(1.0)       },
  { "bc_cyl"            , VAR(bc_cyl)           , PARAM_BOOL(false)       },
  {}
};
#undef VAR

struct mrc_crds_gen_ops mrc_crds_gen_gaussian_ops = {
  .name        = "gaussian",
  .size        = sizeof(struct mrc_crds_gen_gaussian),
  .param_descr = mrc_crds_gen_gaussian_param_descr,
  .run         = mrc_crds_gen_gaussian_run,
};

#undef mrc_crds_gen_gaussian


// ======================================================================
// mrc_crds_gen subclass "two_gaussian"
//
// Nonuniform coordinates initiated via the sum of two Gaussian functions.
// Usefull for double tearing modes.

struct mrc_crds_gen_two_gaussian {
  double gc_x0; // Center (ie point of smallest dx)
  double gc_r; // approximate ratio from peak to min crd
  double gc_w; // Width, i.e. sigma of the gaussian
  double gc_x1;  
  double gc_rx; 
  double gc_wx; 
  bool bc_cyl; // Hack for cylindrical coord bc
};


#define mrc_crds_gen_two_gaussian(gen) mrc_to_subobj(gen, struct mrc_crds_gen_two_gaussian)

static double
f2_dx(struct mrc_crds_gen_two_gaussian *mbc, double x)
{
  double r = mbc->gc_r, w = mbc->gc_w, x0 = mbc->gc_x0, x1 = mbc->gc_x1 ;   
  return r - (r - 1.) * ( exp(-sqr(x-x0)/(2.*sqr(w))) + exp(-sqr(x-x1)/(2.*sqr(w)))   )   ;
}



static void
mrc_crds_gen_two_gaussian_run(struct mrc_crds_gen *gen, double *xx, double *dx)
{
  struct mrc_crds_gen_two_gaussian *sub = mrc_crds_gen_two_gaussian(gen);

  // Calculate a set of grid points
  double *_xx, *ncxx;
  _xx = (double *)calloc(gen->n+1+2*gen->sw, sizeof(*xx));
  ncxx = _xx + gen->sw;
  for (int jx = 0; jx < gen->n + gen->sw; jx++) {
    ncxx[jx+1] = ncxx[jx] + f2_dx(sub, (jx + .5)/gen->n);
  }

  if (sub->bc_cyl) {
    for (int jx = -gen->sw; jx < 0; jx++) {
    ncxx[jx] = -ncxx[-jx];
    }
  } else { 
    for (int jx = 0; jx > -gen->sw; jx--) {
      ncxx[jx-1] = ncxx[jx] - f2_dx(sub, (jx - .5)/gen->n);
    }
  }

  // normalize the domain to [0,1]
  double fac = 1. / ncxx[gen->n];
  for (int jx = -gen->sw; jx<= gen->n + gen->sw; jx++){
    ncxx[jx] *= fac;
  }

  // average them down to the actual coordinates (basically NC to CC, I thnk)
  for (int ii = -gen->sw; ii < gen->n + gen->sw; ii++) {
    xx[ii] = gen->xl + 0.5*(ncxx[ii] + ncxx[ii+1]) * (gen->xh - gen->xl);
  }
  free(_xx);
}


#define VAR(x) (void *)offsetof(struct mrc_crds_gen_two_gaussian, x)
static struct param mrc_crds_gen_two_gaussian_param_descr[] = {
  { "gc_x0"             , VAR(gc_x0)            , PARAM_DOUBLE(0.0)       },
  { "gc_r"              , VAR(gc_r)             , PARAM_DOUBLE(1.0)       },
  { "gc_w"              , VAR(gc_w)             , PARAM_DOUBLE(1.0)       },
  { "gc_x1"             , VAR(gc_x1)            , PARAM_DOUBLE(1.0)       },  
  { "gc_rx"             , VAR(gc_rx)            , PARAM_DOUBLE(1.0)       },
  { "gc_wx"             , VAR(gc_wx)            , PARAM_DOUBLE(1.0)       },
  { "bc_cyl"            , VAR(bc_cyl)           , PARAM_BOOL(false)       },
  {}
};
#undef VAR


// ======================================================================
// mrc_crd_gen subclass "gaussian"

struct mrc_crds_gen_ops mrc_crds_gen_two_gaussian_ops = {
  .name        = "two_gaussian",
  .size        = sizeof(struct mrc_crds_gen_two_gaussian),
  .param_descr = mrc_crds_gen_two_gaussian_param_descr,
  .run         = mrc_crds_gen_two_gaussian_run,
};

#undef mrc_crds_gen_two_gaussian
