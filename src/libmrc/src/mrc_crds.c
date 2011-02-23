
#include <mrc_crds.h>
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_io.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

// ----------------------------------------------------------------------

static inline
struct mrc_crds_ops *mrc_crds_ops(struct mrc_crds *crds)
{
  return (struct mrc_crds_ops *) crds->obj.ops;
}

static inline
struct mrc_crds *to_mrc_crds(struct mrc_obj *obj)
{
  assert(obj->class == &mrc_class_mrc_crds);
  return container_of(obj, struct mrc_crds, obj);
}

// ----------------------------------------------------------------------
// mrc_crds_* wrappers

static void
_mrc_crds_destroy(struct mrc_obj *obj)
{
  struct mrc_crds *crds = to_mrc_crds(obj);
  for (int d = 0; d < 3; d++) {
    mrc_f1_destroy(crds->crd[d]);
  }
}

static void
_mrc_crds_read(struct mrc_obj *obj, struct mrc_io *io)
{
  struct mrc_crds *crds = to_mrc_crds(obj);
  crds->domain = (struct mrc_domain *)
    mrc_io_read_obj_ref(io, mrc_obj_name(obj), "domain", &mrc_class_mrc_domain);
  for (int d = 0; d < 3; d++) {
    char s[5];
    sprintf(s, "crd%d", d);
    crds->crd[d] = (struct mrc_f1 *)
      mrc_io_read_obj_ref(io, mrc_obj_name(obj), s, &mrc_class_mrc_f1);
  }
}

static void
_mrc_crds_write(struct mrc_obj *obj, struct mrc_io *io)
{
  struct mrc_crds *crds = to_mrc_crds(obj);
  mrc_io_write_obj_ref(io, mrc_obj_name(obj), "domain",
		       (struct mrc_obj *) crds->domain);
  for (int d = 0; d < 3; d++) {
    char s[5];
    sprintf(s, "crd%d", d);
    mrc_io_write_obj_ref(io, mrc_obj_name(obj), s, (struct mrc_obj *) crds->crd[d]);
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
  for (int d = 0; d < 3; d++) {
    dx[d] = MRC_CRD(crds, d, 1) - MRC_CRD(crds, d, 0);
  }
}

static void
mrc_crds_alloc(struct mrc_crds *crds, int d, int dim)
{
  mrc_f1_destroy(crds->crd[d]);
  crds->crd[d] = mrc_f1_create(mrc_crds_comm(crds));
  char s[5]; sprintf(s, "crd%d", d);
  mrc_f1_set_name(crds->crd[d], s);
  crds->crd[d]->name[0] = strdup(s);
  mrc_f1_set_param_int(crds->crd[d], "imx", dim);
  mrc_f1_setup(crds->crd[d]);
}

// ======================================================================
// mrc_crds_uniform

static void
mrc_crds_uniform_setup(struct mrc_obj *obj)
{
  struct mrc_crds *crds = to_mrc_crds(obj);
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;

  int sw = crds->par.sw;
  int off[3], ldims[3], gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  mrc_domain_get_local_offset_dims(crds->domain, off, ldims);
  float *xl = crds->par.xl, *xh = crds->par.xh;
  for (int d = 0; d < 3; d++) {
    mrc_crds_alloc(crds, d, ldims[d] + 2 * sw);
    for (int i = -sw; i < ldims[d] +  sw; i++) {
      MRC_CRD(crds, d, i + sw) = xl[d] + (i + off[d] + .5) / gdims[d] * (xh[d] - xl[d]);
    }
  }
}

static struct mrc_crds_ops crds_ops_uniform = {
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
    mrc_crds_alloc(crds, d, m[d]);
    memcpy(crds->crd[d]->arr, crd[d], m[d] * sizeof(*crd[d]));
  }
}

static void
mrc_crds_rectilinear_setup(struct mrc_obj *obj)
{
  struct mrc_crds *crds = to_mrc_crds(obj);
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;

  int sw = crds->par.sw;
  int off[3], ldims[3], gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  mrc_domain_get_local_offset_dims(crds->domain, off, ldims);
  for (int d = 0; d < 3; d++) {
    if (!crds->crd[d]) {
      mrc_crds_alloc(crds, d, ldims[d] + 2 * sw);
    }
  }
}

static struct mrc_crds_ops crds_ops_rectilinear = {
  .name       = "rectilinear",
  .setup      = mrc_crds_rectilinear_setup,
  .set_values = mrc_crds_rectilinear_set_values,
};

// ======================================================================
// mrc_crds class

static LIST_HEAD(mrc_crds_subs);

void
mrc_crds_register(struct mrc_crds_ops *ops)
{
  list_add_tail(&ops->list, &mrc_crds_subs);
}

// ----------------------------------------------------------------------
// mrc_crds_init

static void
mrc_crds_init()
{
  mrc_crds_register(&crds_ops_uniform);
  mrc_crds_register(&crds_ops_rectilinear);
}

#define VAR(x) (void *)offsetof(struct mrc_crds_params, x)
static struct param mrc_crds_params_descr[] = {
  { "l"              , VAR(xl)            , PARAM_FLOAT3(0., 0., 0.) },
  { "h"              , VAR(xh)            , PARAM_FLOAT3(1., 1., 1.) },
  { "sw"             , VAR(sw)            , PARAM_INT(0)             },
  {},
};
#undef VAR

struct mrc_class mrc_class_mrc_crds = {
  .name         = "mrc_crds",
  .size         = sizeof(struct mrc_crds),
  .subclasses   = &mrc_crds_subs,
  .param_descr  = mrc_crds_params_descr,
  .param_offset = offsetof(struct mrc_crds, par),
  .init         = mrc_crds_init,
  .destroy      = _mrc_crds_destroy,
  .write        = _mrc_crds_write,
  .read         = _mrc_crds_read,
};

