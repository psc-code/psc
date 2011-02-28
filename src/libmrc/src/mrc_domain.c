
#include <mrc_domain_private.h>
#include <mrc_fld.h>
#include <mrc_params.h>
#include <mrc_io.h>

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

static inline struct mrc_domain_ops *
mrc_domain_ops(struct mrc_domain *domain)
{
  return (struct mrc_domain_ops *) domain->obj.ops;
}

#define check_is_setup(domain) do { assert(domain->is_setup); } while (0)

// ======================================================================
// mrc_domain

static void
_mrc_domain_create(struct mrc_obj *obj)
{
  struct mrc_domain *domain = to_mrc_domain(obj);
  domain->crds = mrc_crds_create(domain->obj.comm);
  mrc_crds_set_domain(domain->crds, domain);
}

static void
_mrc_domain_destroy(struct mrc_obj *obj)
{
  struct mrc_domain *domain = to_mrc_domain(obj);
  mrc_crds_destroy(domain->crds);
}

static void
_mrc_domain_set_from_options(struct mrc_obj *obj)
{
  struct mrc_domain *domain = to_mrc_domain(obj);
  mrc_crds_set_from_options(domain->crds);
}

static void
_mrc_domain_view(struct mrc_obj *obj)
{
  struct mrc_domain *domain = to_mrc_domain(obj);
  mrc_crds_view(domain->crds);
}

static void
_mrc_domain_setup(struct mrc_obj *obj)
{
  struct mrc_domain *domain = to_mrc_domain(obj);
  mrc_obj_setup_sub(obj);
  mrc_crds_setup(domain->crds);
}

static void
_mrc_domain_read(struct mrc_obj *obj, struct mrc_io *io)
{
  // FIXME, the whole ldims business doesn't work if ldims aren't the same everywhere
  struct mrc_domain *domain = to_mrc_domain(obj);
  mrc_crds_destroy(domain->crds);
  char *s;
  mrc_io_read_attr_string(io, mrc_domain_name(domain), "crds", &s);
  domain->crds = mrc_crds_read(io, s);
  free(s);
}

static void
_mrc_domain_write(struct mrc_obj *obj, struct mrc_io *io)
{
  struct mrc_domain *domain = to_mrc_domain(obj);
  mrc_io_write_attr_string(io, mrc_obj_name(obj), "crds",
			   mrc_crds_name(domain->crds));
  mrc_crds_write(mrc_domain_get_crds(domain), io);
}

struct mrc_patch *
mrc_domain_get_patches(struct mrc_domain *domain, int *nr_patches)
{
  check_is_setup(domain);
  assert(mrc_domain_ops(domain)->get_patches);
  return mrc_domain_ops(domain)->get_patches(domain, nr_patches);
}

struct mrc_crds *
mrc_domain_get_crds(struct mrc_domain *domain)
{
  return domain->crds;
}

int
mrc_domain_get_neighbor_rank(struct mrc_domain *domain, int shift[3])
{
  check_is_setup(domain);
  return mrc_domain_ops(domain)->get_neighbor_rank(domain, shift);
}

void
mrc_domain_get_global_dims(struct mrc_domain *domain, int *dims)
{
  check_is_setup(domain);
  mrc_domain_ops(domain)->get_global_dims(domain, dims);
}

void
mrc_domain_get_local_idx(struct mrc_domain *domain, int *idx)
{
  check_is_setup(domain);
  mrc_domain_ops(domain)->get_local_idx(domain, idx);
}

void
mrc_domain_get_nr_procs(struct mrc_domain *domain, int *nr_procs)
{
  check_is_setup(domain);
  mrc_domain_ops(domain)->get_nr_procs(domain, nr_procs);
}

void
mrc_domain_get_bc(struct mrc_domain *domain, int *bc)
{
  check_is_setup(domain);
  mrc_domain_ops(domain)->get_bc(domain, bc);
}

bool
mrc_domain_is_setup(struct mrc_domain *domain)
{
  return domain->is_setup;
}

// ======================================================================

struct mrc_ddc *
mrc_domain_create_ddc(struct mrc_domain *domain, struct mrc_ddc_params *ddc_par,
		      struct mrc_ddc_ops *ddc_ops)
{
  check_is_setup(domain);
  return mrc_domain_ops(domain)->create_ddc(domain, ddc_par, ddc_ops);
}

// ======================================================================

struct mrc_f3 *
mrc_domain_f3_create(struct mrc_domain *domain, int bnd)
{
  int ldims[3];
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(domain, &nr_patches);
  assert(nr_patches == 1);
  for (int d = 0; d < 3; d++) {
    ldims[d] = patches[0].ldims[d] + 2 * bnd;
  }
  struct mrc_f3 *f3 = mrc_f3_alloc(domain->obj.comm, NULL, ldims);
  f3->domain = domain;
  f3->sw = bnd;
  return f3;
}

// ======================================================================

struct mrc_m3 *
mrc_domain_m3_create(struct mrc_domain *domain, int sw)
{
  struct mrc_m3 *m3 = mrc_m3_create(domain->obj.comm);
  mrc_m3_set_param_int(m3, "sw", sw);
  m3->domain = domain;
  return m3;
}

// ======================================================================
// mrc_domain class

static LIST_HEAD(mrc_domain_subclasses);

void
libmrc_domain_register(struct mrc_domain_ops *ops)
{
  list_add_tail(&ops->list, &mrc_domain_subclasses);
}

// ----------------------------------------------------------------------
// mrc_domain_init

static void
mrc_domain_init()
{
  libmrc_domain_register_simple();
  libmrc_domain_register_multi();
}

struct mrc_class mrc_class_mrc_domain = {
  .name             = "mrc_domain",
  .size             = sizeof(struct mrc_domain),
  .subclasses       = &mrc_domain_subclasses,
  .init             = mrc_domain_init,
  .create           = _mrc_domain_create,
  .destroy          = _mrc_domain_destroy,
  .set_from_options = _mrc_domain_set_from_options,
  .view             = _mrc_domain_view,
  .setup            = _mrc_domain_setup,
  .read             = _mrc_domain_read,
  .write            = _mrc_domain_write,
};

