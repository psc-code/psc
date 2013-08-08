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

static void
mrc_crds_alloc(struct mrc_crds *crds, int d, int dim, int sw)
{
  mrc_fld_destroy(crds->crd[d]);
  crds->crd[d] = mrc_domain_f1_create(crds->domain);
  char s[5]; sprintf(s, "crd%d", d);
  mrc_fld_set_name(crds->crd[d], s);
  mrc_fld_set_sw(crds->crd[d], sw);
  mrc_fld_set_param_int(crds->crd[d], "dim", d);
  mrc_fld_setup(crds->crd[d]);
  mrc_fld_set_comp_name(crds->crd[d], 0, s);
}


static inline
struct mrc_crds_ops *mrc_crds_ops(struct mrc_crds *crds)
{
  return (struct mrc_crds_ops *) crds->obj.ops;
}

// ======================================================================
// mrc_crds_gaussian

// Nonuniform coordinates initiated via a Gaussian function,
// rather than the funky sort used by jimmy_crd
// This is similar to the default version in the original 
// libmrc (ie, the one found in mrc-v3)

struct mrc_crds_gaussian {
  float gc_x0; // Center (ie point of smallest dx)
  float gc_r; // approximate ratio from peak to min crd
  float gc_w; // Width, i.e. sigma of the gaussian
  float gc_x1;  
  float gc_rx; 
  float gc_wx; 
  float gc_y0;
  float gc_ry; 
  float gc_wy;
};

#define VAR(x) (void *)offsetof(struct mrc_crds_gaussian, x)
static struct param mrc_crds_gaussian_param_descr[] = {
  { "gc_x0"             , VAR(gc_x0)            , PARAM_FLOAT(0.0)       },
  { "gc_r"              , VAR(gc_r)             , PARAM_FLOAT(1.0)       },
  { "gc_w"              , VAR(gc_w)             , PARAM_FLOAT(1.0)       },
  { "gc_x1"             , VAR(gc_x1)            , PARAM_FLOAT(1.0)       },  
  { "gc_rx"              , VAR(gc_rx)             , PARAM_FLOAT(1.0)       },
  { "gc_wx"              , VAR(gc_wx)             , PARAM_FLOAT(1.0)       },
  { "gc_y0"             , VAR(gc_y0)            , PARAM_FLOAT(0.0)       },
  { "gc_ry"              , VAR(gc_ry)             , PARAM_FLOAT(1.0)       },
  { "gc_wy"              , VAR(gc_wy)             , PARAM_FLOAT(1.0)       },
  {}
};
#undef VAR

#define to_mrc_crds_gaussian(crds) \
  mrc_to_subobj(crds, struct mrc_crds_gaussian)

static double
f_dx(struct mrc_crds_gaussian *mbc, double x)
{
  double r = mbc->gc_r, w = mbc->gc_w, x0 = mbc->gc_x0;
  return r - (r - 1.) * exp(-sqr(x-x0)/(2.*sqr(w)));
}


static double
f_2dx(struct mrc_crds_gaussian *mbc, double x)
{
  double r = mbc->gc_rx, w = mbc->gc_wx, x0 = mbc->gc_x0;
  return r - (r - 1.) * exp(-sqr(x-x0)/(2.*sqr(w)));
}

static double
f_2dy(struct mrc_crds_gaussian *mbc, double x)
{
  double r = mbc->gc_ry, w = mbc->gc_wy, x0 = mbc->gc_y0;
  return r - (r - 1.) * exp(-sqr(x-x0)/(2.*sqr(w)));
}


static double
f2_dx(struct mrc_crds_gaussian *mbc, double x)
{
  double r = mbc->gc_r, w = mbc->gc_w, x0 = mbc->gc_x0, x1 = mbc->gc_x1 ;   
  return r - (r - 1.) * ( exp(-sqr(x-x0)/(2.*sqr(w))) + exp(-sqr(x-x1)/(2.*sqr(w)))   )   ;
}

static void
mrc_crds_gaussian_setup(struct mrc_crds *crds)
{
  struct mrc_crds_gaussian *gauss = to_mrc_crds_gaussian(crds);
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;
  int sw = crds->sw;
  float *xl = crds->xl, *xh = crds->xh;


  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(crds->domain, &nr_patches);
  assert(nr_patches == 1);

  float *_xx, *xx;
  _xx = (float *)calloc(gdims[0]+1+2*sw, sizeof(*xx));
  xx = _xx + sw;
  for (int jx = 0; jx < gdims[0] + sw; jx++) {
    xx[jx+1] = xx[jx] + f_dx(gauss, (jx + .5)/gdims[0]);
  }
  for (int jx = 0; jx > -sw; jx--) {
    xx[jx-1] = xx[jx] - f_dx(gauss, (jx - .5)/gdims[0]);
  }
  double fac = 1. / xx[gdims[0]];
  for (int jx = -sw; jx<= gdims[0] + sw; jx++){
    xx[jx] *= fac;
  }
  for (int d = 1; d < 3; d++) {
    mrc_crds_alloc(crds, d, patches[0].ldims[d], sw);    
    for (int i = -sw; i < patches[0].ldims[d] +  sw; i++) {
      int ii =  i+patches[0].off[d];
      // accumulate grid and make sure to normalize to intended domain size.
      // fine tuning required at user end to ensure that grid refinement profile 
      // suitable for physically relevant profile..
      MRC_CRD(crds, d, i) = xl[d] + ((float) ii+0.5) / ((float) gdims[d]) * (xh[d] - xl[d]);
    }
  }
  mrc_crds_alloc(crds, 0, patches[0].ldims[0], sw);
  for (int i = -sw; i < patches[0].ldims[0] + sw; i++) {
    // +patches[0].off[0] added to 
    int ii =  i+patches[0].off[0];
    MRC_CRD(crds, 0, i) = xl[0] + 0.5*(xx[ii] + xx[ii+1]) * (xh[0] - xl[0]);
  }
}


static void
mrc_crds_two_gaussian_setup(struct mrc_crds *crds)
{
  struct mrc_crds_gaussian *gauss = to_mrc_crds_gaussian(crds);
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;
  int sw = crds->sw;
  float *xl = crds->xl, *xh = crds->xh;
  

  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(crds->domain, &nr_patches);
  assert(nr_patches == 1);
  // FIXME, parallel (look at version mrc-v3, already parallel there)
  // Also, this is just gaussian in x, and uniform in the other three directions
  // it would probably nice to fix that somehow.
  float *_xx, *xx;
  _xx = (float *)calloc(gdims[0]+1+2*sw, sizeof(*xx));
  xx = _xx + sw;
  for (int jx = 0; jx < gdims[0] + sw; jx++) {
    xx[jx+1] = xx[jx] + f2_dx(gauss, (jx + .5)/gdims[0]);
  }
  for (int jx = 0; jx > -sw; jx--) {
    xx[jx-1] = xx[jx] - f2_dx(gauss, (jx - .5)/gdims[0]);
  }
  double fac = 1. / xx[gdims[0]];
  for (int jx = -sw; jx<= gdims[0] + sw; jx++){
    xx[jx] *= fac;
  }
  for (int d = 1; d < 3; d++) {
    mrc_crds_alloc(crds, d, patches[0].ldims[d], sw);    
    for (int i = -sw; i < patches[0].ldims[d] +  sw; i++) {
      int ii =  i+patches[0].off[d];
      MRC_CRD(crds, d, i) = xl[d] + ((float) ii+0.5) / ((float) gdims[d]) * (xh[d] - xl[d]);
    }
  }
  mrc_crds_alloc(crds, 0, patches[0].ldims[0], sw);
  for (int i = -sw; i < patches[0].ldims[0] + sw; i++) {
    // patches --> 
    
    int ii =  i+patches[0].off[0];
    MRC_CRD(crds, 0, i) = xl[0] + 0.5*(xx[ii] + xx[ii+1]) * (xh[0] - xl[0]);
  }
}


static void
mrc_crds_gaussian_2D_setup(struct mrc_crds *crds)
{
  struct mrc_crds_gaussian *gauss = to_mrc_crds_gaussian(crds);
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;
  int sw = crds->sw;
  float *xl = crds->xl, *xh = crds->xh;

  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(crds->domain, &nr_patches);
  assert(nr_patches == 1);

  // populate grid coordinate arrays
  float *_xx, *xx;
  _xx = (float *)calloc(gdims[0]+1+2*sw, sizeof(*xx));
  xx = _xx + sw;
  for (int jx = 0; jx < gdims[0] + sw; jx++) {
    xx[jx+1] = xx[jx] + f_2dx(gauss, (jx + .5)/gdims[0]);
  }
  for (int jx = 0; jx > -sw; jx--) {
    xx[jx-1] = xx[jx] - f_2dx(gauss, (jx - .5)/gdims[0]);
  }
  double facx = 1. / xx[gdims[0]];
  for (int jx = -sw; jx<= gdims[0] + sw; jx++){
    xx[jx] *= facx;
  }
  float *_yy, *yy;
  _yy = (float *)calloc(gdims[1]+1+2*sw, sizeof(*yy));
  yy = _yy + sw;
  for (int jx = 0; jx < gdims[1] + sw; jx++) {
    yy[jx+1] = yy[jx] + f_2dy(gauss, (jx + .5)/gdims[1]);
  }
  for (int jx = 0; jx > -sw; jx--) {
    yy[jx-1] = yy[jx] - f_2dy(gauss, (jx - .5)/gdims[1]);
  }
  double facy = 1. / yy[gdims[1]];
  for (int jx = -sw; jx<= gdims[1] + sw; jx++){
    yy[jx] *= facy;
  }
  
  
  for (int d = 1; d < 3; d++) {
    mrc_crds_alloc(crds, d, patches[0].ldims[d], sw);    
    for (int i = -sw; i < patches[0].ldims[d] +  sw; i++) {
      int ii =  i+patches[0].off[d];
      // accumulate grid and make sure to normalize to intended domain size.
      // fine tuning required at user end to ensure that grid refinement profile 
      // suitable for physically relevant profile..
      MRC_CRD(crds, d, i) = xl[d] + ((float) ii+0.5) / ((float) gdims[d]) * (xh[d] - xl[d]);
    }
  }
  mrc_crds_alloc(crds, 0, patches[0].ldims[0], sw);
  for (int i = -sw; i < patches[0].ldims[0] + sw; i++) {
    // +patches[0].off[0] added to 
    int ii =  i+patches[0].off[0];
    MRC_CRD(crds, 0, i) = xl[0] + 0.5*(xx[ii] + xx[ii+1]) * (xh[0] - xl[0]);
  }
  
  mrc_crds_alloc(crds, 1, patches[0].ldims[1], sw);
  for (int i = -sw; i < patches[0].ldims[1] + sw; i++) {
    // +patches[0].off[0] added to 
    int ii =  i+patches[0].off[1];
    MRC_CRD(crds, 1, i) = xl[1] + 0.5*(yy[ii] + yy[ii+1]) * (xh[1] - xl[1]);
  }
  

}




struct mrc_crds_ops mrc_crds_gaussian_ops = {
  .name        = "gaussian",
  .size        = sizeof(struct mrc_crds_gaussian),
  .param_descr = mrc_crds_gaussian_param_descr,
  .setup       = mrc_crds_gaussian_setup,
};

struct mrc_crds_ops mrc_crds_two_gaussian_ops = {
  .name        = "two_gaussian",
  .size        = sizeof(struct mrc_crds_gaussian),
  .param_descr = mrc_crds_gaussian_param_descr,
  .setup       = mrc_crds_two_gaussian_setup,
};

struct mrc_crds_ops mrc_crds_gaussian_2D_ops = {
  .name        = "gaussian_2D",
  .size        = sizeof(struct mrc_crds_gaussian),
  .param_descr = mrc_crds_gaussian_param_descr,
  .setup       = mrc_crds_gaussian_2D_setup,
};
