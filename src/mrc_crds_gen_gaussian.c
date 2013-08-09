
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
  float gc_x0; // Center (ie point of smallest dx)
  float gc_r; // approximate ratio from peak to min crd
  float gc_w; // Width, i.e. sigma of the gaussian
};

#define mrc_crds_gen_gaussian(gen) mrc_to_subobj(gen, struct mrc_crds_gen_gaussian)

static double
f_dx(struct mrc_crds_gen_gaussian *mbc, double x)
{
  double r = mbc->gc_r, w = mbc->gc_w, x0 = mbc->gc_x0;
  return r - (r - 1.) * exp(-sqr(x-x0)/(2.*sqr(w)));
}


static void
mrc_crds_gen_gaussian_run(struct mrc_crds_gen *gen, float *xx, float *dx)
{
  struct mrc_crds_gen_gaussian *sub = mrc_crds_gen_gaussian(gen);

  int gdims[3];
  mrc_domain_get_global_dims(gen->crds->domain, gdims);

  float xl[3], xh[3];
  mrc_crds_get_param_float3(gen->crds, "l", xl);
  mrc_crds_get_param_float3(gen->crds, "h", xh);

  int sw;
  mrc_crds_get_param_int(gen->crds, "sw", &sw);

  int d = gen->d;

  // Calculate a set of grid points
  float *_xx, *ncxx;
  _xx = (float *)calloc(gdims[d]+1+2*sw, sizeof(*xx));
  ncxx = _xx + sw;
  for (int jx = 0; jx < gdims[d] + sw; jx++) {
    ncxx[jx+1] = ncxx[jx] + f_dx(sub, (jx + .5)/gdims[d]);
  }
  for (int jx = 0; jx > -sw; jx--) {
    ncxx[jx-1] = ncxx[jx] - f_dx(sub, (jx - .5)/gdims[d]);
  }
  double fac = 1. / ncxx[gdims[d]];
  for (int jx = -sw; jx<= gdims[d] + sw; jx++){
    ncxx[jx] *= fac;
  }

  // average them down to the actual coordinates (basically NC to CC)
  for (int ii = -sw; ii < gdims[d] + sw; ii++) {
    xx[ii] = xl[d] + 0.5*(ncxx[ii] + ncxx[ii+1]) * (xh[d] - xl[d]);
  }
  free(_xx);
}


#define VAR(x) (void *)offsetof(struct mrc_crds_gen_gaussian, x)
static struct param mrc_crds_gen_gaussian_param_descr[] = {
  { "gc_x0"             , VAR(gc_x0)            , PARAM_FLOAT(0.0)       },
  { "gc_r"              , VAR(gc_r)             , PARAM_FLOAT(1.0)       },
  { "gc_w"              , VAR(gc_w)             , PARAM_FLOAT(1.0)       },
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
  float gc_x0; // Center (ie point of smallest dx)
  float gc_r; // approximate ratio from peak to min crd
  float gc_w; // Width, i.e. sigma of the gaussian
  float gc_x1;  
  float gc_rx; 
  float gc_wx; 
};


#define mrc_crds_gen_two_gaussian(gen) mrc_to_subobj(gen, struct mrc_crds_gen_two_gaussian)

static double
f2_dx(struct mrc_crds_gen_two_gaussian *mbc, double x)
{
  double r = mbc->gc_r, w = mbc->gc_w, x0 = mbc->gc_x0, x1 = mbc->gc_x1 ;   
  return r - (r - 1.) * ( exp(-sqr(x-x0)/(2.*sqr(w))) + exp(-sqr(x-x1)/(2.*sqr(w)))   )   ;
}



static void
mrc_crds_gen_two_gaussian_run(struct mrc_crds_gen *gen, float *xx, float *dx)
{
  struct mrc_crds_gen_two_gaussian *sub = mrc_crds_gen_two_gaussian(gen);

  int gdims[3];
  mrc_domain_get_global_dims(gen->crds->domain, gdims);

  float xl[3], xh[3];
  mrc_crds_get_param_float3(gen->crds, "l", xl);
  mrc_crds_get_param_float3(gen->crds, "h", xh);

  int sw;
  mrc_crds_get_param_int(gen->crds, "sw", &sw);

  int d = gen->d;
  
  // Calculate a set of grid points
  float *_xx, *ncxx;
  _xx = (float *)calloc(gdims[d]+1+2*sw, sizeof(*xx));
  ncxx = _xx + sw;
  for (int jx = 0; jx < gdims[d] + sw; jx++) {
    ncxx[jx+1] = ncxx[jx] + f2_dx(sub, (jx + .5)/gdims[d]);
  }
  for (int jx = 0; jx > -sw; jx--) {
    ncxx[jx-1] = ncxx[jx] - f2_dx(sub, (jx - .5)/gdims[d]);
  }

  // normalize the domain to [0,1]
  double fac = 1. / ncxx[gdims[d]];
  for (int jx = -sw; jx<= gdims[d] + sw; jx++){
    ncxx[jx] *= fac;
  }

  // average them down to the actual coordinates (basically NC to CC, I thnk)
  for (int ii = -sw; ii < gdims[d] + sw; ii++) {
    xx[ii] = xl[d] + 0.5*(ncxx[ii] + ncxx[ii+1]) * (xh[d] - xl[d]);
  }
  free(_xx);
}


#define VAR(x) (void *)offsetof(struct mrc_crds_gen_two_gaussian, x)
static struct param mrc_crds_gen_two_gaussian_param_descr[] = {
  { "gc_x0"             , VAR(gc_x0)            , PARAM_FLOAT(0.0)       },
  { "gc_r"              , VAR(gc_r)             , PARAM_FLOAT(1.0)       },
  { "gc_w"              , VAR(gc_w)             , PARAM_FLOAT(1.0)       },
  { "gc_x1"             , VAR(gc_x1)            , PARAM_FLOAT(1.0)       },  
  { "gc_rx"             , VAR(gc_rx)            , PARAM_FLOAT(1.0)       },
  { "gc_wx"             , VAR(gc_wx)            , PARAM_FLOAT(1.0)       },
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
