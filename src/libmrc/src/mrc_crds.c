
#include "mrc_crds_gen_private.h"
#include <mrc_crds.h>
#include <mrc_crds_gen.h>
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_ndarray.h>
#include <mrc_io.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

// static inline
// struct mrc_crds_ops *mrc_crds_ops(struct mrc_crds *crds)
// {
//   return (struct mrc_crds_ops *) crds->obj.ops;
// }

// ----------------------------------------------------------------------
// mrc_crds_* wrappers

static void
_mrc_crds_create(struct mrc_crds *crds)
{
  for (int d = 0; d < 3; d++) {
    char s[20];

    sprintf(s, "crd[%d]", d);
    mrc_fld_set_name(crds->crd[d], s);

    sprintf(s, "dcrd[%d]", d);
    mrc_fld_set_name(crds->dcrd[d], s);
    
    // These aren't listed as member objects, so we need to create them ourselves here
    // because we do something hacky when writing them for xdmf_collective
    crds->crd_nc[d] = mrc_fld_create(mrc_crds_comm(crds));
    sprintf(s, "crd_nc[%d]", d);
    mrc_fld_set_name(crds->crd_nc[d], s);

    sprintf(s, "dcrd_nc[%d]", d);
    mrc_fld_set_name(crds->dcrd_nc[d], s);

    crds->global_crd[d] = mrc_ndarray_create(mrc_crds_comm(crds));
    sprintf(s, "global_crd[%d]", d);
    mrc_ndarray_set_name(crds->global_crd[d], s);

    sprintf(s, "crds_gen_%c", 'x' + d);
    mrc_crds_gen_set_name(crds->crds_gen[d], s);
    mrc_crds_gen_set_param_int(crds->crds_gen[d], "d", d);
    mrc_crds_gen_set_param_obj(crds->crds_gen[d], "crds", crds);
  }
}

static void
_mrc_crds_read(struct mrc_crds *crds, struct mrc_io *io)
{
  if (strcmp(mrc_io_type(io), "hdf5_serial") != 0 &&
      strcmp(mrc_io_type(io), "xdmf_serial") != 0) {
    // FIXME, but reading back coordinates is broken for everything but hdf5_serial,
    // because for other mrc_io types, we don't write crd_nc (at least not completely),
    // so there's no (easy) way to restore it.
    assert(0);
  }

  mrc_crds_read_member_objs(crds, io);

  // this is a carbon copy on all nodes that run the crd_gen, so this
  // write is only for checkpointing
  if (strcmp(mrc_io_type(io), "hdf5_serial") == 0) {
    crds->crd_nc[0] = mrc_io_read_ref(io, crds, "crd_nc[0]", mrc_fld);
    crds->crd_nc[1] = mrc_io_read_ref(io, crds, "crd_nc[1]", mrc_fld);
    crds->crd_nc[2] = mrc_io_read_ref(io, crds, "crd_nc[2]", mrc_fld);

    crds->global_crd[0] = mrc_io_read_ref(io, crds, "global_crd[0]", mrc_ndarray);
    crds->global_crd[1] = mrc_io_read_ref(io, crds, "global_crd[1]", mrc_ndarray);
    crds->global_crd[2] = mrc_io_read_ref(io, crds, "global_crd[2]", mrc_ndarray);    
  }

  // FIXME, shouldn't the generic read() take care of this?
  crds->obj.is_setup = true;
}

// ----------------------------------------------------------------------
// _mrc_crds_write

static void
_mrc_crds_write(struct mrc_crds *crds, struct mrc_io *io)
{
  // FIXME, hacky to no end...
  // So for xdmf_collective, we need to write crd_nc while chopping of ghost points,
  // which breaks completely reading them back
  // for checkpointing/hdf5_serial, however, we write them including ghost points,
  // just like every other coordinate field...
  // ...except that every other coordinate field is listed as a member object, so it happens
  // automatically, but here we need to do it explicitly.

  if (strcmp(mrc_io_type(io), "xdmf_collective") == 0) { // FIXME
    int slab_off_save[3], slab_dims_save[3];

    mrc_io_get_param_int3(io, "slab_off", slab_off_save);
    mrc_io_get_param_int3(io, "slab_dims", slab_dims_save);

    int gdims[3];
    mrc_domain_get_global_dims(crds->domain, gdims);
    for (int d = 0; d < 3; d++) {
      // FIXME, this is really too hacky... should per m1 / m3, not per mrc_io
      struct mrc_fld *crd_nc = crds->crd_nc[d];
      mrc_io_set_param_int3(io, "slab_off", (int[3]) { 0, 0, 0});
      int slab_dims[3] = {};
      slab_dims[d] = gdims[d] + 1;
      mrc_io_set_param_int3(io, "slab_dims", slab_dims);
      mrc_fld_write(crd_nc, io);
    }

    mrc_io_set_param_int3(io, "slab_off", slab_off_save);
    mrc_io_set_param_int3(io, "slab_dims", slab_dims_save);
  }

  if (strcmp(mrc_io_type(io), "hdf5_serial") == 0) { // FIXME
    for (int d = 0; d < 3; d++) {
      struct mrc_fld *crd_nc = crds->crd_nc[d];
      mrc_io_write_ref(io, crds, mrc_fld_name(crd_nc), crd_nc);
    }
    
    // this is a carbon copy on all nodes that run the crd_gen, so this
    // write is only for checkpointing
    for (int d = 0; d < 3; d++) {
      struct mrc_ndarray *global_crd = crds->global_crd[d];
      mrc_io_write_ref(io, crds, mrc_ndarray_name(global_crd), global_crd);
    }
  }
}

void
mrc_crds_get_dx(struct mrc_crds *crds, int p, double dx[3])
{
  // FIXME, only for uniform crds, should be dispatched!
  dx[0] = MRC_DMCRDX(crds, 1, p) - MRC_DMCRDX(crds, 0, p);
  dx[1] = MRC_DMCRDY(crds, 1, p) - MRC_DMCRDY(crds, 0, p);
  dx[2] = MRC_DMCRDZ(crds, 1, p) - MRC_DMCRDZ(crds, 0, p);
}

// FIXME, should go away / superseded by mrc_crds_get_dx()
// calculate and return 0th level dx for amr
void
mrc_crds_get_dx_base(struct mrc_crds *crds, double dx[3])
{
  if (strcmp(mrc_crds_type(crds), "amr_uniform") == 0) {
    int lm[3];
    const double *lo = mrc_crds_lo(crds), *hi = mrc_crds_hi(crds);
    mrc_domain_get_param_int3(crds->domain, "m", lm);
    for (int d = 0; d < 3; d++) {
      dx[d] = (hi[d] - lo[d]) / lm[d];
    }
  } else {
    assert(strcmp(mrc_crds_type(crds), "uniform") == 0);
    // the only place where this makes sense is if we have one patch / proc, no AMR, anyway
    mrc_crds_get_dx(crds, 0, dx); // the only 
  }
}

// allocate the coordinate fields common to all crds types.
static void
mrc_crds_setup_alloc_only(struct mrc_crds *crds)
{
  assert(crds->domain && mrc_domain_is_setup(crds->domain));

  for (int d = 0; d < 3; d++) {
    mrc_fld_set_param_obj(crds->crd[d], "domain", crds->domain);
    mrc_fld_set_param_int(crds->crd[d], "nr_spatial_dims", 1);
    mrc_fld_set_param_int(crds->crd[d], "dim", d);
    mrc_fld_set_param_int(crds->crd[d], "nr_ghosts", crds->sw);
    mrc_fld_set_comp_name(crds->crd[d], 0, mrc_fld_name(crds->crd[d]));
    mrc_fld_dict_add_double(crds->crd[d], "io_scale", crds->xnorm);
    mrc_fld_setup(crds->crd[d]);

    // double version of coords
    mrc_fld_set_type(crds->dcrd[d], "double");
    mrc_fld_set_param_obj(crds->dcrd[d], "domain", crds->domain);
    mrc_fld_set_param_int(crds->dcrd[d], "nr_spatial_dims", 1);
    mrc_fld_set_param_int(crds->dcrd[d], "dim", d);
    mrc_fld_set_param_int(crds->dcrd[d], "nr_ghosts", crds->sw);
    mrc_fld_set_comp_name(crds->dcrd[d], 0, mrc_fld_name(crds->dcrd[d]));
    mrc_fld_setup(crds->dcrd[d]);

    // node-centered coords
    mrc_fld_set_param_obj(crds->crd_nc[d], "domain", crds->domain);
    mrc_fld_set_param_int(crds->crd_nc[d], "nr_spatial_dims", 1);
    mrc_fld_set_param_int(crds->crd_nc[d], "dim", d);
    mrc_fld_set_param_int(crds->crd_nc[d], "nr_ghosts", crds->sw + 1);
    mrc_fld_set_comp_name(crds->crd_nc[d], 0, mrc_fld_name(crds->crd_nc[d]));
    mrc_fld_dict_add_double(crds->crd_nc[d], "io_scale", crds->xnorm);
    mrc_fld_setup(crds->crd_nc[d]);

    // node-centered coords
    mrc_fld_set_type(crds->dcrd_nc[d], "double");
    mrc_fld_set_param_obj(crds->dcrd_nc[d], "domain", crds->domain);
    mrc_fld_set_param_int(crds->dcrd_nc[d], "nr_spatial_dims", 1);
    mrc_fld_set_param_int(crds->dcrd_nc[d], "dim", d);
    mrc_fld_set_param_int(crds->dcrd_nc[d], "nr_ghosts", crds->sw + 1);
    mrc_fld_set_comp_name(crds->dcrd_nc[d], 0, mrc_fld_name(crds->dcrd_nc[d]));
    mrc_fld_setup(crds->dcrd_nc[d]);
  }
}

// Allocate global coordinate fields (for domains that they make sense)
static void
mrc_crds_setup_alloc_global_array(struct mrc_crds *crds)
{
  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  
  for (int d = 0; d < 3; d++) {
    mrc_ndarray_set_type(crds->global_crd[d], "double");
    mrc_ndarray_set_param_int_array(crds->global_crd[d], "dims", 2, (int[2]) { gdims[d] + 2 * crds->sw, 2 });
    mrc_ndarray_set_param_int_array(crds->global_crd[d], "offs"  , 2, (int[2]) { -crds->sw, 0 });
    mrc_ndarray_setup(crds->global_crd[d]);
  }
}

// ----------------------------------------------------------------------
// _mrc_crds_setup

static void
_mrc_crds_setup(struct mrc_crds *crds)
{
  int gdims[3];
  int nr_patches;
  struct mrc_patch *patches;

  crds->xnorm = crds->norm_length / crds->norm_length_scale;

  for (int d = 0; d < 3; d ++) {
    crds->lo_code[d] = crds->l[d] / crds->xnorm;
    crds->hi_code[d] = crds->h[d] / crds->xnorm;
  }
  
  mrc_crds_setup_alloc_only(crds);
  mrc_crds_setup_alloc_global_array(crds);

  mrc_domain_get_global_dims(crds->domain, gdims);
  patches = mrc_domain_get_patches(crds->domain, &nr_patches);

  int sw = crds->sw;

  for (int d = 0; d < 3; d ++) {
    struct mrc_ndarray *x = crds->global_crd[d];

    struct mrc_crds_gen *gen = crds->crds_gen[d];

    mrc_crds_gen_set_param_int(gen, "n", gdims[d]);
    mrc_crds_gen_set_param_int(gen, "sw", gen->crds->sw);
    mrc_crds_gen_set_param_double(gen, "xl", crds->lo_code[d]);
    mrc_crds_gen_set_param_double(gen, "xh", crds->hi_code[d]);
    mrc_crds_gen_run(gen, &MRC_D2(x, 0, 0), &MRC_D2(x, 0, 1));

    mrc_fld_foreach_patch(crds->crd[d], p) {
      // shift to beginning of local domain
      int off = patches[p].off[d];

      mrc_m1_foreach(crds->crd[d], i, sw, sw) {
        MRC_DMCRD(crds, d, i, p) = MRC_D2(x, i + off, 0);
        MRC_MCRD(crds, d, i, p) = MRC_DMCRD(crds, d, i, p);
      } mrc_m1_foreach_end;

      mrc_m1_foreach(crds->crd[d], i, sw, sw + 1) {
	if (i + off == -sw) { // extrapolate on low side
	  MRC_DMCRD_NC(crds, d, i, p) = MRC_D2(x, i + off, 0)
	    - .5 * (MRC_D2(x, i+1 + off, 0) - MRC_D2(x, i + off, 0));
	} else if (i + off == gdims[d] + sw) { // extrapolate on high side
	  MRC_DMCRD_NC(crds, d, i, p) = MRC_D2(x, i-1 + off, 0)
	    + .5 * (MRC_D2(x, i-1 + off, 0) - MRC_D2(x, i-2 + off, 0));
	} else {
	  MRC_DMCRD_NC(crds, d, i, p) = .5 * (MRC_D2(x, i-1 + off, 0) + MRC_D2(x, i + off, 0));
	}
	MRC_MCRD_NC(crds, d, i, p) = MRC_DMCRD_NC(crds, d, i, p);
      } mrc_m1_foreach_end;
    }
  }
}

static void
_mrc_crds_destroy(struct mrc_crds *crds)
{
  for (int d = 0; d < 3; d++) {
    mrc_ndarray_destroy(crds->global_crd[d]);
    mrc_fld_destroy(crds->crd_nc[d]);
  }
}

// ----------------------------------------------------------------------
// mrc_crds_lo

const double *
mrc_crds_lo(struct mrc_crds *crds)
{
  assert(crds->obj.is_setup);
  return crds->lo_code;
}

// ----------------------------------------------------------------------
// mrc_crds_hi

const double *
mrc_crds_hi(struct mrc_crds *crds)
{
  assert(crds->obj.is_setup);
  return crds->hi_code;
}

// ======================================================================
// mrc_crds_uniform

static void
mrc_crds_uniform_setup(struct mrc_crds *crds)
{
  for (int d = 0; d < 3; d++) {
    assert(strcmp(mrc_crds_gen_type(crds->crds_gen[d]), "uniform") == 0);
  }

  mrc_crds_setup_super(crds);
}

static struct mrc_crds_ops mrc_crds_uniform_ops = {
  .name  = "uniform",
  .setup = mrc_crds_uniform_setup,
};

// ======================================================================
// mrc_crds_rectilinear

static struct mrc_crds_ops mrc_crds_rectilinear_ops = {
  .name       = "rectilinear",
};

// ======================================================================
// mrc_crds_amr_uniform

// FIXME, this should use mrc_a1 not mrc_m1

static void
mrc_crds_amr_uniform_setup(struct mrc_crds *crds)
{
  mrc_crds_setup_alloc_only(crds);

  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);
  const double *lo = mrc_crds_lo(crds), *hi = mrc_crds_hi(crds);

  for (int d = 0; d < 3; d++) {
    struct mrc_fld *mcrd = crds->crd[d];
    struct mrc_fld *dcrd = crds->dcrd[d];
    mrc_m1_foreach_patch(mcrd, p) {
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(crds->domain, p, &info);
      double xb = (double) info.off[d] / (1 << info.level);
      double xe = (double) (info.off[d] + info.ldims[d]) / (1 << info.level);
      double dx = (xe - xb) / info.ldims[d];

      mrc_m1_foreach_bnd(mcrd, i) {
        MRC_D3(dcrd,i, 0, p) = lo[d] + (xb + (i + .5) * dx) / gdims[d] * (hi[d] - lo[d]);
        MRC_M1(mcrd,0, i, p) = (float)MRC_D3(dcrd,i, 0, p);
      } mrc_m1_foreach_end;
    }
  }
}

static struct mrc_crds_ops mrc_crds_amr_uniform_ops = {
  .name  = "amr_uniform",
  .setup = mrc_crds_amr_uniform_setup,
};


// ======================================================================
// mrc_crds_mb
// We need a new type of coords to be able to handle multi-block domains,
// since they don't generally have a global coordinate system. We'll iterate
// through and create each block local coordinate system.


static void
mrc_crds_mb_create(struct mrc_crds *crds)
{
  // Name the standard 3 crds_gen to "UNUSED" to (hopefully) avoid some confusion
  for (int d = 0; d < 3; d++) {
    mrc_crds_gen_set_name(crds->crds_gen[d], "UNUSED");
  }
}

static void
mrc_crds_mb_read(struct mrc_crds *crds, struct mrc_io *io)
{
  // need to make sure subclass create doesn't get called at read, but
  // the superclass read does.
  mrc_crds_read_super(crds, io);
}

static void
mrc_crds_mb_setup(struct mrc_crds *crds)
{
  // this should still work
  mrc_crds_setup_alloc_only(crds);

  typedef void (*dgb_t)(struct mrc_domain *, struct MB_block **pblock, int *nr_blocks);
  dgb_t domain_get_blocks = (dgb_t) mrc_domain_get_method(crds->domain, "get_blocks");

  int nr_blocks;
  struct MB_block *blocks;

  domain_get_blocks(crds->domain, &blocks, &nr_blocks);

  int nr_patches;
  mrc_domain_get_patches(crds->domain, &nr_patches);
  int sw = crds->sw;

  for (int b = 0; b < nr_blocks; b++)
    {
      struct MB_block *block = &(blocks[b]);

      for (int d = 0; d < 3; d ++) {
        struct mrc_fld *x = mrc_fld_create(MPI_COMM_SELF);
        mrc_fld_set_type(x, "double");
        mrc_fld_set_param_int_array(x, "dims", 2, (int[2]) { block->mx[d] + 1, 2 });
        mrc_fld_set_param_int_array(x, "sw"  , 2, (int[2]) { sw, 0 });
        mrc_fld_setup(x);

        // If I had my way, I'd kill off the original crds_gen children to minimize confusion,
        // but I guess I'll just have to write the docs to make it clear what's going on here.
        assert(block->coord_gen[d]);

        mrc_crds_gen_set_param_int(block->coord_gen[d], "n", block->mx[d]);
        mrc_crds_gen_set_param_int(block->coord_gen[d], "d", d);
        mrc_crds_gen_set_param_int(block->coord_gen[d], "sw", sw);
        mrc_crds_gen_set_param_obj(block->coord_gen[d], "crds", crds);
        mrc_crds_gen_set_param_double(block->coord_gen[d], "xl", block->xl[d]);
        mrc_crds_gen_set_param_double(block->coord_gen[d], "xh", block->xh[d]);

        mrc_crds_gen_run(block->coord_gen[d], &MRC_D2(x, 0, 0), &MRC_D2(x, 0, 1));
    
        mrc_m1_foreach_patch(crds->crd[d], p) {
          struct mrc_patch_info info;
          mrc_domain_get_local_patch_info(crds->domain, p, &info);
          if (b == info.p_block) {
            // This is offset of the patch in the block
            int off = info.p_ix[d];
            mrc_m1_foreach_bnd(crds->crd[d], ix) {
              MRC_DMCRD(crds, d, ix, p) = MRC_D2(x, ix + off, 0);
              MRC_MCRD(crds, d, ix, p) = (float)MRC_D2(x, ix + off, 0);
          } mrc_m1_foreach_end;
        }
      }
      mrc_fld_destroy(x);
    }
  }
}


static struct mrc_crds_ops mrc_crds_mb_ops = {
  .name      = "mb",
  .create    = mrc_crds_mb_create,
  .setup     = mrc_crds_mb_setup,
  .read      = mrc_crds_mb_read,
};

// ======================================================================
// mrc_crds_init

extern struct mrc_crds_ops mrc_crds_two_gaussian_ops;
extern struct mrc_crds_ops mrc_crds_gaussian_ops;
extern struct mrc_crds_ops mrc_crds_gaussian_2D_ops;

static void
mrc_crds_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_uniform_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_rectilinear_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_amr_uniform_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_mb_ops);
}

// ======================================================================
// mrc_crds class

#define VAR(x) (void *)offsetof(struct mrc_crds, x)
static struct param mrc_crds_params_descr[] = {
  { "l"                , VAR(l)                , PARAM_DOUBLE3(0., 0., 0.) },
  { "h"                , VAR(h)                , PARAM_DOUBLE3(1., 1., 1.) },
  { "sw"               , VAR(sw)               , PARAM_INT(0)              },
  { "norm_length"      , VAR(norm_length)      , PARAM_DOUBLE(1.)          },
  { "norm_length_scale", VAR(norm_length_scale), PARAM_DOUBLE(1.)          },
  { "domain"           , VAR(domain)           , PARAM_OBJ(mrc_domain)     },

  { "xnorm"          , VAR(xnorm)         , MRC_VAR_DOUBLE           },
  { "lo_code"        , VAR(lo_code)       , MRC_VAR_DOUBLE3          },
  { "hi_code"        , VAR(hi_code)       , MRC_VAR_DOUBLE3          },

  { "crd[0]"         , VAR(crd[0])        , MRC_VAR_OBJ(mrc_fld)     },
  { "crd[1]"         , VAR(crd[1])        , MRC_VAR_OBJ(mrc_fld)     },
  { "crd[2]"         , VAR(crd[2])        , MRC_VAR_OBJ(mrc_fld)     },

  { "dcrd[0]"        , VAR(dcrd[0])       , MRC_VAR_OBJ(mrc_fld)     },
  { "dcrd[1]"        , VAR(dcrd[1])       , MRC_VAR_OBJ(mrc_fld)     },
  { "dcrd[2]"        , VAR(dcrd[2])       , MRC_VAR_OBJ(mrc_fld)     },

  /* { "crd_nc[0]"      , VAR(crd_nc[0])     , MRC_VAR_OBJ(mrc_fld)     }, */
  /* { "crd_nc[1]"      , VAR(crd_nc[1])     , MRC_VAR_OBJ(mrc_fld)     }, */
  /* { "crd_nc[2]"      , VAR(crd_nc[2])     , MRC_VAR_OBJ(mrc_fld)     }, */

  { "dcrd_nc[0]"     , VAR(dcrd_nc[0])    , MRC_VAR_OBJ(mrc_fld)     },
  { "dcrd_nc[1]"     , VAR(dcrd_nc[1])    , MRC_VAR_OBJ(mrc_fld)     },
  { "dcrd_nc[2]"     , VAR(dcrd_nc[2])    , MRC_VAR_OBJ(mrc_fld)     },

  { "crds_gen_x"     , VAR(crds_gen[0])   , MRC_VAR_OBJ(mrc_crds_gen)},
  { "crds_gen_y"     , VAR(crds_gen[1])   , MRC_VAR_OBJ(mrc_crds_gen)},
  { "crds_gen_z"     , VAR(crds_gen[2])   , MRC_VAR_OBJ(mrc_crds_gen)},

  {},
};
#undef VAR

struct mrc_class_mrc_crds mrc_class_mrc_crds = {
  .name         = "mrc_crds",
  .size         = sizeof(struct mrc_crds),
  .param_descr  = mrc_crds_params_descr,
  .init         = mrc_crds_init,
  .destroy      = _mrc_crds_destroy,
  .create       = _mrc_crds_create,
  .write        = _mrc_crds_write,
  .read         = _mrc_crds_read,
  .setup        = _mrc_crds_setup,
};

