
#include <mrc_crds.h>
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_io.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

static inline
struct mrc_crds_ops *mrc_crds_ops(struct mrc_crds *crds)
{
  return (struct mrc_crds_ops *) crds->obj.ops;
}

// ----------------------------------------------------------------------
// mrc_crds_* wrappers

static void
_mrc_crds_destroy(struct mrc_crds *crds)
{
  for (int d = 0; d < 3; d++) {
    mrc_f1_destroy(crds->crd[d]);
    mrc_m1_destroy(crds->mcrd[d]);
  }
}

static void
_mrc_crds_read(struct mrc_crds *crds, struct mrc_io *io)
{
  crds->domain = (struct mrc_domain *)
    mrc_io_read_obj_ref(io, mrc_crds_name(crds), "domain", &mrc_class_mrc_domain);
  if (strcmp(mrc_crds_type(crds), "multi_uniform") == 0 ||
      strcmp(mrc_crds_type(crds), "multi_rectilinear") == 0) {
    for (int d = 0; d < 3; d++) {
      char s[6];
      sprintf(s, "mcrd%d", d);
      crds->mcrd[d] = (struct mrc_m1 *)
	mrc_io_read_obj_ref(io, mrc_crds_name(crds), s, &mrc_class_mrc_m1);
    }
  } else {
    for (int d = 0; d < 3; d++) {
      char s[5];
      sprintf(s, "crd%d", d);
      crds->crd[d] = (struct mrc_f1 *)
	mrc_io_read_obj_ref(io, mrc_crds_name(crds), s, &mrc_class_mrc_f1);
    }
  }
  mrc_crds_setup(crds);
}

static void
_mrc_crds_write(struct mrc_crds *crds, struct mrc_io *io)
{
  mrc_io_write_obj_ref(io, mrc_crds_name(crds), "domain",
		       (struct mrc_obj *) crds->domain);
  for (int d = 0; d < 3; d++) {
    if (crds->crd[d]) {
      char s[5];
      sprintf(s, "crd%d", d);
      mrc_io_write_obj_ref(io, mrc_crds_name(crds), s, (struct mrc_obj *) crds->crd[d]);
    }
    if (crds->mcrd[d]) {
      char s[6];
      sprintf(s, "mcrd%d", d);
      mrc_io_write_obj_ref(io, mrc_crds_name(crds), s, (struct mrc_obj *) crds->mcrd[d]);
    }
  }
}

void
mrc_crds_set_domain(struct mrc_crds *crds, struct mrc_domain *domain)
{
  crds->domain = domain;
}

void
mrc_crds_set_values(struct mrc_crds *crds, float *crdx, int mx,
		    float *crdy, int my, float *crdz, int mz)
{
  struct mrc_crds_ops *ops = mrc_crds_ops(crds);
  if (ops->set_values) {
    ops->set_values(crds, crdx, mx, crdy, my, crdz, mz);
  }
}

void
mrc_crds_get_xl_xh(struct mrc_crds *crds, float xl[3], float xh[3])
{
  if (xl) {
    for (int d = 0; d < 3; d++) {
      xl[d] = crds->par.xl[d];
    }
  }
  if (xh) {
    for (int d = 0; d < 3; d++) {
      xh[d] = crds->par.xh[d];
    }
  }
}

void
mrc_crds_get_dx(struct mrc_crds *crds, float dx[3])
{
  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  // FIXME, only makes sense for uniform coords, should be dispatched!!!
  for (int d = 0; d < 3; d++) {
    dx[d] = (crds->par.xh[d] - crds->par.xl[d]) / gdims[d];
  }
}

static void
mrc_crds_alloc(struct mrc_crds *crds, int d, int dim, int sw)
{
  mrc_f1_destroy(crds->crd[d]);
  crds->crd[d] = mrc_domain_f1_create(crds->domain);
  char s[5]; sprintf(s, "crd%d", d);
  mrc_f1_set_name(crds->crd[d], s);
  mrc_f1_set_param_int(crds->crd[d], "sw", sw);
  mrc_f1_set_param_int(crds->crd[d], "dim", d);
  mrc_f1_setup(crds->crd[d]);
  mrc_f1_set_comp_name(crds->crd[d], 0, s);
}

static void
mrc_crds_multi_alloc(struct mrc_crds *crds, int d)
{
  mrc_m1_destroy(crds->mcrd[d]);
  crds->mcrd[d] = mrc_domain_m1_create(crds->domain);
  char s[5]; sprintf(s, "crd%d", d);
  mrc_m1_set_name(crds->mcrd[d], s);
  mrc_m1_set_param_int(crds->mcrd[d], "sw", crds->par.sw);
  mrc_m1_set_param_int(crds->mcrd[d], "dim", d);
  mrc_m1_setup(crds->mcrd[d]);
  mrc_m1_set_comp_name(crds->mcrd[d], 0, s);
}

void
mrc_crds_patch_get(struct mrc_crds *crds, int p)
{
  for (int d = 0; d < 3; d++) {
    crds->mcrd_p[d] = mrc_m1_patch_get(crds->mcrd[d], p);
  }
}

void
mrc_crds_patch_put(struct mrc_crds *crds)
{
  for (int d = 0; d < 3; d++) {
    crds->mcrd_p[d] = NULL;
  }
}

// ======================================================================
// mrc_crds_uniform

static void
mrc_crds_uniform_setup(struct mrc_crds *crds)
{
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;

  int sw = crds->par.sw;
  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(crds->domain, &nr_patches);
  assert(nr_patches == 1);
  float *xl = crds->par.xl, *xh = crds->par.xh;
  for (int d = 0; d < 3; d++) {
    mrc_crds_alloc(crds, d, patches[0].ldims[d], sw);
    for (int i = -sw; i < patches[0].ldims[d] +  sw; i++) {
      MRC_CRD(crds, d, i) = xl[d] + (i + patches[0].off[d] + .5) / gdims[d] * (xh[d] - xl[d]);
    }
  }
}

static struct mrc_crds_ops mrc_crds_uniform_ops = {
  .name  = "uniform",
  .setup = mrc_crds_uniform_setup,
};

// ======================================================================
// mrc_crds_rectilinear

static void
mrc_crds_rectilinear_set_values(struct mrc_crds *crds, float *crdx, int mx,
				float *crdy, int my, float *crdz, int mz)
{
  float *crd[3] = { crdx, crdy, crdz };
  int m[3] = { mx, my, mz };
  for (int d = 0; d < 3; d++) {
    mrc_crds_alloc(crds, d, m[d], crds->par.sw);
    memcpy(crds->crd[d]->arr, crd[d], crds->crd[d]->len * sizeof(*crd[d]));
  }
}

static void
mrc_crds_rectilinear_setup(struct mrc_crds *crds)
{
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;

  int sw = crds->par.sw;
  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(crds->domain, &nr_patches);
  assert(nr_patches == 1);
  for (int d = 0; d < 3; d++) {
    if (!crds->crd[d]) {
      mrc_crds_alloc(crds, d, patches[0].ldims[d], sw);
    }
  }
}

static struct mrc_crds_ops mrc_crds_rectilinear_ops = {
  .name       = "rectilinear",
  .setup      = mrc_crds_rectilinear_setup,
  .set_values = mrc_crds_rectilinear_set_values,
};

// ======================================================================
// mrc_crds_rectilinear_jr2

struct mrc_crds_rectilinear_jr2 {
  float xm;
  float xn;
  float dx0;
};

#define VAR(x) (void *)offsetof(struct mrc_crds_rectilinear_jr2, x)
static struct param mrc_crds_rectilinear_jr2_param_descr[] = {
  { "dx0"             , VAR(dx0)            , PARAM_FLOAT(.1)       },
  { "xn"              , VAR(xn)             , PARAM_FLOAT(2.)       },
  { "xm"              , VAR(xm)             , PARAM_FLOAT(.5)       },
  {},
};
#undef VAR

#define to_mrc_crds_rectilinear_jr2(crds) \
  mrc_to_subobj(crds, struct mrc_crds_rectilinear_jr2)

static void
mrc_crds_rectilinear_jr2_setup(struct mrc_crds *crds)
{
  struct mrc_crds_rectilinear_jr2 *jr2 = to_mrc_crds_rectilinear_jr2(crds);
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;

  int sw = crds->par.sw;
  float *xl = crds->par.xl, *xh = crds->par.xh;

  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(crds->domain, &nr_patches);
  assert(nr_patches == 1);

  // FIXME, parallel
  float xm = jr2->xm, xn = jr2->xn, dx0 = jr2->dx0;
  for (int d = 0; d < 3; d++) {
    assert(-xl[d] == xh[d]);
    float Lx2 = xh[d];
    float xc = gdims[d] / 2 - .5;
    float a = (pow((Lx2) / (dx0 * xc), 1./xm) - 1.) / pow(xc, 2.*xn);
    mrc_crds_alloc(crds, d, patches[0].ldims[d], sw);
    for (int i = -sw; i < patches[0].ldims[d] +  sw; i++) {
      float xi = i - xc;
      float s = 1 + a*(pow(xi, (2. * xn)));
      float sm = pow(s, xm);
      float g = dx0 * xi * sm;
      //    float dg = rmhd->dx0 * (sm + rmhd->xm*xi*2.*rmhd->xn*a*(pow(xi, (2.*rmhd->xn-1.))) * sm / s);
      MRC_CRD(crds, d, i) = g;
    }
  }
}

static struct mrc_crds_ops mrc_crds_rectilinear_jr2_ops = {
  .name        = "rectilinear_jr2",
  .size        = sizeof(struct mrc_crds_rectilinear_jr2),
  .param_descr = mrc_crds_rectilinear_jr2_param_descr,
  .setup       = mrc_crds_rectilinear_jr2_setup,
};

// ======================================================================
// mrc_crds_multi_uniform

static void
mrc_crds_multi_uniform_setup(struct mrc_crds *crds)
{
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;

  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  float *xl = crds->par.xl, *xh = crds->par.xh;

  struct mrc_patch *patches = mrc_domain_get_patches(crds->domain, NULL);
  for (int d = 0; d < 3; d++) {
    mrc_crds_multi_alloc(crds, d);
    struct mrc_m1 *mcrd = crds->mcrd[d];
    mrc_m1_foreach_patch(mcrd, p) {
      struct mrc_m1_patch *mcrd_p = mrc_m1_patch_get(mcrd, p);
      mrc_m1_foreach_bnd(mcrd_p, i) {
	MRC_M1(mcrd_p,0, i) = xl[d] + (i + patches[p].off[d] + .5) / gdims[d] * (xh[d] - xl[d]);
      } mrc_m1_foreach_end;
      mrc_m1_patch_put(mcrd);
    }
  }
}

static struct mrc_crds_ops mrc_crds_multi_uniform_ops = {
  .name  = "multi_uniform",
  .setup = mrc_crds_multi_uniform_setup,
};

// ======================================================================
// mrc_crds_multi_rectilinear

static void
mrc_crds_multi_rectilinear_setup(struct mrc_crds *crds)
{
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;

  for (int d = 0; d < 3; d++) {
    if (!crds->mcrd[d]) {
      mrc_crds_multi_alloc(crds, d);
    }
  }
}

static struct mrc_crds_ops mrc_crds_multi_rectilinear_ops = {
  .name  = "multi_rectilinear",
  .setup = mrc_crds_multi_rectilinear_setup,
};

// ======================================================================
// mrc_crds_init

static void
mrc_crds_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_uniform_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_rectilinear_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_rectilinear_jr2_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_multi_uniform_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_multi_rectilinear_ops);
}

// ======================================================================
// mrc_crds class

#define VAR(x) (void *)offsetof(struct mrc_crds_params, x)
static struct param mrc_crds_params_descr[] = {
  { "l"              , VAR(xl)            , PARAM_FLOAT3(0., 0., 0.) },
  { "h"              , VAR(xh)            , PARAM_FLOAT3(1., 1., 1.) },
  { "sw"             , VAR(sw)            , PARAM_INT(0)             },
  {},
};
#undef VAR

struct mrc_class_mrc_crds mrc_class_mrc_crds = {
  .name         = "mrc_crds",
  .size         = sizeof(struct mrc_crds),
  .param_descr  = mrc_crds_params_descr,
  .param_offset = offsetof(struct mrc_crds, par),
  .init         = mrc_crds_init,
  .destroy      = _mrc_crds_destroy,
  .write        = _mrc_crds_write,
  .read         = _mrc_crds_read,
};

