
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
    mrc_fld_destroy(crds->crd[d]);
  }
}

static void
_mrc_crds_read(struct mrc_crds *crds, struct mrc_io *io)
{
  crds->domain = mrc_io_read_ref(io, crds, "domain", mrc_domain);
  for (int d = 0; d < 3; d++) {
    char s[5];
    sprintf(s, "crd%d", d);
    crds->crd[d] = mrc_io_read_ref(io, crds, s, mrc_fld);
  }
  mrc_crds_setup(crds);
}

static void
_mrc_crds_write(struct mrc_crds *crds, struct mrc_io *io)
{
  mrc_io_write_ref(io, crds, "domain", crds->domain);
  for (int d = 0; d < 3; d++) {
    int slab_off_save[3], slab_dims_save[3];
    if (strcmp(mrc_io_type(io), "xdmf_collective") == 0) { // FIXME
      mrc_io_get_param_int3(io, "slab_off", slab_off_save);
      mrc_io_get_param_int3(io, "slab_dims", slab_dims_save);
      mrc_io_set_param_int3(io, "slab_off", (int[3]) { 0, 0, 0});
      mrc_io_set_param_int3(io, "slab_dims", (int[3]) { 0, 0, 0 });
    }

    struct mrc_fld *crd_cc = crds->crd[d];
    char s[10];
    sprintf(s, "crd%d", d);
    mrc_io_write_ref(io, crds, s, crd_cc);

    if (strcmp(mrc_io_type(io), "xdmf_collective") == 0) { // FIXME
      struct mrc_fld *crd_nc = crds->crd_nc[d];
      if (!crd_nc) {
	crd_nc = mrc_fld_create(mrc_crds_comm(crds)); // FIXME, leaked
	crds->crd_nc[d] = crd_nc;
	char s[10];
	sprintf(s, "crd%d_nc", d);
	mrc_fld_set_name(crd_nc, s);
	mrc_fld_set_param_obj(crd_nc, "domain", crds->domain);
	mrc_fld_set_param_int(crd_nc, "dim", d);
	mrc_fld_set_param_int_array(crd_nc, "dims", 3, NULL);
	mrc_fld_set_sw(crd_nc, 1);
	mrc_fld_set_nr_comps(crd_nc, 1);
	mrc_fld_setup(crd_nc);
	mrc_fld_set_comp_name(crd_nc, 0, s);

	mrc_m1_foreach_patch(crd_nc, p) {
	  if (crds->sw > 0) {
	    mrc_m1_foreach(crd_cc, i, 0, 1) {
	      MRC_M1(crd_nc,0, i, p) = .5 * (MRC_M1(crd_cc,0, i-1, p) + MRC_M1(crd_cc,0, i, p));
	    } mrc_m1_foreach_end;
	  } else {
	    mrc_m1_foreach(crd_cc, i, -1, 0) {
	      MRC_M1(crd_nc,0, i, p) = .5 * (MRC_M1(crd_cc,0, i-1, p) + MRC_M1(crd_cc,0, i, p));
	    } mrc_m1_foreach_end;
	    int ld = mrc_fld_dims(crd_nc)[0];
	    // extrapolate
	    MRC_M1(crd_nc,0, 0 , p) = MRC_M1(crd_cc,0, 0   , p) 
	      - .5 * (MRC_M1(crd_cc,0,    1, p) - MRC_M1(crd_cc,0, 0   , p));
	    MRC_M1(crd_nc,0, ld, p) = MRC_M1(crd_cc,0, ld-1, p)
	      + .5 * (MRC_M1(crd_cc,0, ld-1, p) - MRC_M1(crd_cc,0, ld-2, p));
	  }
	}
      }
      int gdims[3];
      mrc_domain_get_global_dims(crds->domain, gdims);
      // FIXME, this is really too hacky... should per m1 / m3, not per mrc_io
      mrc_io_set_param_int3(io, "slab_off", (int[3]) { 0, 0, 0});
      mrc_io_set_param_int3(io, "slab_dims", (int[3]) { gdims[d] + 1, 0, 0 });
      mrc_fld_write(crd_nc, io);
    }
    if (strcmp(mrc_io_type(io), "xdmf_collective") == 0) { // FIXME
      mrc_io_set_param_int3(io, "slab_off", slab_off_save);
      mrc_io_set_param_int3(io, "slab_dims", slab_dims_save);
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
  assert(mrc_crds_is_setup(crds));
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
      xl[d] = crds->xl[d];
    }
  }
  if (xh) {
    for (int d = 0; d < 3; d++) {
      xh[d] = crds->xh[d];
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
    dx[d] = (crds->xh[d] - crds->xl[d]) / gdims[d];
  }
}

static void
mrc_crds_alloc(struct mrc_crds *crds)
{
  for (int d = 0; d < 3; d++) {
    mrc_fld_destroy(crds->crd[d]);
    crds->crd[d] = mrc_domain_m1_create(crds->domain);
    char s[5]; sprintf(s, "crd%d", d);
    mrc_fld_set_name(crds->crd[d], s);
    mrc_fld_set_sw(crds->crd[d], crds->sw);
    mrc_fld_set_param_int(crds->crd[d], "dim", d);
    mrc_fld_setup(crds->crd[d]);
    mrc_fld_set_comp_name(crds->crd[d], 0, s);
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

  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  float *xl = crds->xl, *xh = crds->xh;

  mrc_crds_alloc(crds);
  struct mrc_patch *patches = mrc_domain_get_patches(crds->domain, NULL);
  for (int d = 0; d < 3; d++) {
    struct mrc_fld *mcrd = crds->crd[d];
    mrc_m1_foreach_patch(mcrd, p) {
      mrc_m1_foreach_bnd(mcrd, i) {
	MRC_M1(mcrd,0, i, p) = xl[d] + (i + patches[p].off[d] + .5) / gdims[d] * (xh[d] - xl[d]);
      } mrc_m1_foreach_end;
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
  for (int d = 0; d < 3; d++) {
    memcpy(crds->crd[d]->_arr, crd[d] - crds->sw, crds->crd[d]->_len * sizeof(*crd[d]));
  }
}

static void
mrc_crds_rectilinear_setup(struct mrc_crds *crds)
{
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;

  assert(!crds->crd[0]);
  mrc_crds_alloc(crds);
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

  int sw = crds->sw;
  float *xl = crds->xl, *xh = crds->xh;

  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(crds->domain, &nr_patches);
  assert(nr_patches == 1);

  // FIXME, parallel
  float xm = jr2->xm, xn = jr2->xn, dx0 = jr2->dx0;
  mrc_crds_alloc(crds);
  for (int d = 0; d < 3; d++) {
    assert(-xl[d] == xh[d]);
    float Lx2 = xh[d];
    float xc = gdims[d] / 2 - .5;
    float a = (pow((Lx2) / (dx0 * xc), 1./xm) - 1.) / pow(xc, 2.*xn);
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
// mrc_crds_amr_uniform

// FIXME, this should use mrc_a1 not mrc_m1

static void
mrc_crds_amr_uniform_setup(struct mrc_crds *crds)
{
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;

  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  float *xl = crds->xl, *xh = crds->xh;

  mrc_crds_alloc(crds);
  for (int d = 0; d < 3; d++) {
    struct mrc_fld *mcrd = crds->crd[d];
    mrc_m1_foreach_patch(mcrd, p) {
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(crds->domain, p, &info);
      float xb = (float) info.off[d] / (1 << info.level);
      float xe = (float) (info.off[d] + info.ldims[d]) / (1 << info.level);
      float dx = (xe - xb) / info.ldims[d];

      mrc_m1_foreach_bnd(mcrd, i) {
	MRC_M1(mcrd,0, i, p) = xl[d] + (xb + (i + .5) * dx) / gdims[d] * (xh[d] - xl[d]);
      } mrc_m1_foreach_end;
    }
  }
}

static struct mrc_crds_ops mrc_crds_amr_uniform_ops = {
  .name  = "amr_uniform",
  .setup = mrc_crds_amr_uniform_setup,
};

// ======================================================================
// mrc_crds_multi_rectilinear

static void
mrc_crds_multi_rectilinear_setup(struct mrc_crds *crds)
{
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;

  assert(!crds->crd[0]);
  mrc_crds_alloc(crds);
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
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_multi_rectilinear_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_amr_uniform_ops);
}

// ======================================================================
// mrc_crds class

#define VAR(x) (void *)offsetof(struct mrc_crds, x)
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
  .init         = mrc_crds_init,
  .destroy      = _mrc_crds_destroy,
  .write        = _mrc_crds_write,
  .read         = _mrc_crds_read,
};

