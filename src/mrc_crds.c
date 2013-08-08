
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
_mrc_crds_read(struct mrc_crds *crds, struct mrc_io *io)
{
  mrc_crds_read_member_objs(crds, io);
}

static void
_mrc_crds_write(struct mrc_crds *crds, struct mrc_io *io)
{
  int slab_off_save[3], slab_dims_save[3];
  if (strcmp(mrc_io_type(io), "xdmf_collective") == 0) { // FIXME
    mrc_io_get_param_int3(io, "slab_off", slab_off_save);
    mrc_io_get_param_int3(io, "slab_dims", slab_dims_save);
    mrc_io_set_param_int3(io, "slab_off", (int[3]) { 0, 0, 0});
    mrc_io_set_param_int3(io, "slab_dims", (int[3]) { 0, 0, 0 });
  }

  for (int d = 0; d < 3; d++) {
    struct mrc_fld *crd_cc = crds->crd[d];
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
  }

  if (strcmp(mrc_io_type(io), "xdmf_collective") == 0) { // FIXME
    mrc_io_set_param_int3(io, "slab_off", slab_off_save);
    mrc_io_set_param_int3(io, "slab_dims", slab_dims_save);
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
    char s[7]; sprintf(s, "crd[%d]", d);
    mrc_fld_set_name(crds->crd[d], s);
    mrc_fld_set_param_obj(crds->crd[d], "domain", crds->domain);
    mrc_fld_set_param_int_array(crds->crd[d], "dims", 3, NULL);
    mrc_fld_set_param_int(crds->crd[d], "dim", d);
    mrc_fld_set_nr_comps(crds->crd[d], 1);
    mrc_fld_set_sw(crds->crd[d], crds->sw);
    mrc_fld_set_comp_name(crds->crd[d], 0, mrc_fld_name(crds->crd[d]));
    mrc_fld_setup(crds->crd[d]);
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
mrc_crds_rectilinear_setup(struct mrc_crds *crds)
{
  assert(crds->domain);
  if (!mrc_domain_is_setup(crds->domain))
    return;

  mrc_crds_alloc(crds);
}

static struct mrc_crds_ops mrc_crds_rectilinear_ops = {
  .name       = "rectilinear",
  .setup      = mrc_crds_rectilinear_setup,
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
// mrc_crds_init

static void
mrc_crds_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_uniform_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_rectilinear_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_amr_uniform_ops);
}

// ======================================================================
// mrc_crds class

#define VAR(x) (void *)offsetof(struct mrc_crds, x)
static struct param mrc_crds_params_descr[] = {
  { "l"              , VAR(xl)            , PARAM_FLOAT3(0., 0., 0.) },
  { "h"              , VAR(xh)            , PARAM_FLOAT3(1., 1., 1.) },
  { "sw"             , VAR(sw)            , PARAM_INT(0)             },
  { "domain"         , VAR(domain)        , PARAM_OBJ(mrc_domain)    },

  { "crd[0]"         , VAR(crd[0])        , MRC_VAR_OBJ(mrc_fld)     },
  { "crd[1]"         , VAR(crd[1])        , MRC_VAR_OBJ(mrc_fld)     },
  { "crd[2]"         , VAR(crd[2])        , MRC_VAR_OBJ(mrc_fld)     },

  {},
};
#undef VAR

struct mrc_class_mrc_crds mrc_class_mrc_crds = {
  .name         = "mrc_crds",
  .size         = sizeof(struct mrc_crds),
  .param_descr  = mrc_crds_params_descr,
  .init         = mrc_crds_init,
  .write        = _mrc_crds_write,
  .read         = _mrc_crds_read,
};

