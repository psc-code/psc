
#include <mrc_domain_private.h>
#include <mrc_fld.h>
#include <mrc_params.h>
#include <mrc_io.h>
#include <mrc_ddc.h>

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

static inline struct mrc_domain_ops *
mrc_domain_ops(struct mrc_domain *domain)
{
  return (struct mrc_domain_ops *) domain->obj.ops;
}

// ======================================================================
// mrc_domain

static void
_mrc_domain_create(struct mrc_domain *domain)
{
  mrc_crds_set_param_obj(domain->crds, "domain", domain);

  // FIXME, this isn't really general, though ok if we always use
  // multi for whatever domain. Otherwise, we need a call back for set type to
  // updated the sub type accordingly... Similar problem exists with crds, btw.
  mrc_ddc_set_type(domain->ddc, "multi");
  mrc_ddc_set_domain(domain->ddc, domain);
  mrc_ddc_set_param_int(domain->ddc, "size_of_type", sizeof(float));
  mrc_ddc_set_funcs(domain->ddc, &mrc_ddc_funcs_fld);
}

static void
_mrc_domain_setup(struct mrc_domain *domain)
{
  mrc_domain_setup_member_objs(domain);
}

static void
_mrc_domain_read(struct mrc_domain *domain, struct mrc_io *io)
{
  // FIXME, mrc_obj should do this by itself
  domain->obj.is_setup = true;
  mrc_domain_read_member_objs(domain, io);
}

static void
_mrc_domain_write(struct mrc_domain *domain, struct mrc_io *io)
{
  mrc_io_write_int(io, domain, "mpi_size", domain->size);
}

int
mrc_domain_nr_patches(struct mrc_domain *domain)
{
  int nr_patches;
  assert(mrc_domain_is_setup(domain));
  assert(mrc_domain_ops(domain)->get_patches);
  mrc_domain_ops(domain)->get_patches(domain, &nr_patches);
  
  return nr_patches;
}

struct mrc_patch *
mrc_domain_get_patches(struct mrc_domain *domain, int *nr_patches)
{
  assert(mrc_domain_is_setup(domain));
  assert(mrc_domain_ops(domain)->get_patches);
  return mrc_domain_ops(domain)->get_patches(domain, nr_patches);
}

struct mrc_crds *
mrc_domain_get_crds(struct mrc_domain *domain)
{
  return domain->crds;
}

struct mrc_ddc *
mrc_domain_get_ddc(struct mrc_domain *domain)
{
  return domain->ddc;
}

int
mrc_domain_get_neighbor_rank(struct mrc_domain *domain, int shift[3])
{
  assert(mrc_domain_is_setup(domain));
  assert(mrc_domain_ops(domain)->get_neighbor_rank);
  return mrc_domain_ops(domain)->get_neighbor_rank(domain, shift);
}

void
mrc_domain_get_global_dims(struct mrc_domain *domain, int *dims)
{
  assert(mrc_domain_is_setup(domain));
  assert(mrc_domain_ops(domain)->get_global_dims);
  mrc_domain_ops(domain)->get_global_dims(domain, dims);
}

void
mrc_domain_get_nr_procs(struct mrc_domain *domain, int *nr_procs)
{
  assert(mrc_domain_is_setup(domain));
  assert(mrc_domain_ops(domain)->get_nr_procs);
  mrc_domain_ops(domain)->get_nr_procs(domain, nr_procs);
}

void
mrc_domain_get_bc(struct mrc_domain *domain, int *bc)
{
  assert(mrc_domain_is_setup(domain));
  
  for (int d = 0; d < 3; d++) {
    bc[d] = domain->bc[d];
  }  
}

void
mrc_domain_get_nr_global_patches(struct mrc_domain *domain, int *nr_global_patches)
{
  assert(mrc_domain_is_setup(domain));
  assert(mrc_domain_ops(domain)->get_nr_global_patches);
  mrc_domain_ops(domain)->get_nr_global_patches(domain, nr_global_patches);
}

void
mrc_domain_get_global_patch_info(struct mrc_domain *domain, int gp,
				 struct mrc_patch_info *info)
{
  assert(mrc_domain_is_setup(domain));
  assert(mrc_domain_ops(domain)->get_global_patch_info);
  mrc_domain_ops(domain)->get_global_patch_info(domain, gp, info);
}

void
mrc_domain_get_local_patch_info(struct mrc_domain *domain, int p,
				struct mrc_patch_info *info)
{
  assert(mrc_domain_is_setup(domain));
  assert(mrc_domain_ops(domain)->get_local_patch_info);
  mrc_domain_ops(domain)->get_local_patch_info(domain, p, info);
}

void
mrc_domain_get_level_idx3_patch_info(struct mrc_domain *domain, int level,
				     int idx[3], struct mrc_patch_info *info)
{
  assert(mrc_domain_is_setup(domain));
  assert(mrc_domain_ops(domain)->get_level_idx3_patch_info);
  mrc_domain_ops(domain)->get_level_idx3_patch_info(domain, level, idx, info);
}

void
mrc_domain_get_nr_levels(struct mrc_domain *domain, int *p_nr_levels)
{
  assert(mrc_domain_is_setup(domain));
  if (!mrc_domain_ops(domain)->get_nr_levels) {
    *p_nr_levels = 1;
    return;
  }
  mrc_domain_ops(domain)->get_nr_levels(domain, p_nr_levels);
}

void
mrc_domain_get_neighbor_rank_patch(struct mrc_domain *domain, int p, int dir[3],
				   int *nei_rank, int *nei_patch)
{
  struct mrc_domain_ops *ops = mrc_domain_ops(domain);
  assert(ops->get_neighbor_rank_patch);
  ops->get_neighbor_rank_patch(domain, p, dir, nei_rank, nei_patch);
}

void
mrc_domain_plot(struct mrc_domain *domain)
{
  assert(mrc_domain_is_setup(domain));
  assert(mrc_domain_ops(domain)->plot);
  mrc_domain_ops(domain)->plot(domain);
}

void
mrc_domain_add_patch(struct mrc_domain *domain, int l, int idx3[3])
{
  assert(mrc_domain_ops(domain)->add_patch);
  mrc_domain_ops(domain)->add_patch(domain, l, idx3);
}

// ======================================================================

struct mrc_ddc *
mrc_domain_create_ddc(struct mrc_domain *domain)
{
  assert(mrc_domain_ops(domain)->create_ddc);
  return mrc_domain_ops(domain)->create_ddc(domain);
}

// ======================================================================

struct mrc_fld *
mrc_domain_f1_create(struct mrc_domain *domain)
{
  int nr_patches;
  mrc_domain_get_patches(domain, &nr_patches);
  assert(nr_patches == 1);
  return mrc_domain_m1_create(domain);
}

// ======================================================================
// mrc_domain_m1_create

struct mrc_fld *
mrc_domain_m1_create(struct mrc_domain *domain)
{
  struct mrc_fld *m1 = mrc_fld_create(domain->obj.comm);
  mrc_fld_set_param_obj(m1, "domain", domain);
  mrc_fld_set_param_int(m1, "nr_spatial_dims", 1);
  // default direction to 0 == x
  mrc_fld_set_param_int(m1, "dim", 0);
  return m1;
}

// ======================================================================
// mrc_domain_fld_create

struct mrc_fld *
mrc_domain_fld_create(struct mrc_domain *domain, int sw, const char *comps)
{
  struct mrc_fld *fld = mrc_fld_create(mrc_domain_comm(domain));
  mrc_fld_set_param_obj(fld, "domain", domain);
  mrc_fld_set_param_int(fld, "nr_spatial_dims", 3);
  if (comps) {
    mrc_fld_set_comp_names(fld, comps);
  }
  mrc_fld_set_param_int(fld, "nr_ghosts", sw);
  return fld;
}

// ======================================================================
// mrc_domain_m3_create

struct mrc_fld *
mrc_domain_m3_create(struct mrc_domain *domain)
{
  struct mrc_fld *m3 = mrc_fld_create(mrc_domain_comm(domain));
  mrc_fld_set_param_obj(m3, "domain", domain);
  mrc_fld_set_param_int(m3, "nr_spatial_dims", 3);
  return m3;
}

// ======================================================================
// mrc_domain_init

static void
mrc_domain_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_domain, &mrc_domain_simple_ops);
  mrc_class_register_subclass(&mrc_class_mrc_domain, &mrc_domain_multi_ops);
  mrc_class_register_subclass(&mrc_class_mrc_domain, &mrc_domain_amr_ops);
#ifdef HAVE_PETSC
  mrc_class_register_subclass(&mrc_class_mrc_domain, &mrc_domain_mb_ops);
#endif
}

// ======================================================================
// mrc_domain class

static struct mrc_param_select bc_descr[] = {
  { .val = BC_NONE       , .str = "none"     },
  { .val = BC_PERIODIC   , .str = "periodic" },
  {},
};

#define VAR(x) (void *)offsetof(struct mrc_domain, x)
static struct param mrc_domain_descr[] = {
  { "crds"           , VAR(crds)          , MRC_VAR_OBJ(mrc_crds)           },
  { "ddc"            , VAR(ddc)           , MRC_VAR_OBJ(mrc_ddc)            },
  { "bcx"            , VAR(bc[0])         , PARAM_SELECT(BC_NONE, bc_descr) },
  { "bcy"            , VAR(bc[1])         , PARAM_SELECT(BC_NONE, bc_descr) },
  { "bcz"            , VAR(bc[2])         , PARAM_SELECT(BC_NONE, bc_descr) },
  {},
};
#undef VAR

struct mrc_class_mrc_domain mrc_class_mrc_domain = {
  .name             = "mrc_domain",
  .size             = sizeof(struct mrc_domain),
  .param_descr      = mrc_domain_descr,
  .init             = mrc_domain_init,
  .create           = _mrc_domain_create,
  .setup            = _mrc_domain_setup,
  .read             = _mrc_domain_read,
  .write            = _mrc_domain_write,
};
