
#include "psc.h"

#include "psc_particles_as_c.h" // only for checking...

#include <mrc_profile.h>
#include <mrc_params.h>
#include <mrc_io.h>
#include <stdlib.h>
#include <string.h>

// ======================================================================

// FIXME, do we need this at all, and if so, could be a PARAM_OBJ
void
psc_mparticles_set_domain(struct psc_mparticles *mparticles,
			  struct mrc_domain *domain)
{
  mparticles->domain = domain;
  mrc_domain_get_patches(domain, &mparticles->nr_patches);
}

static void
_psc_mparticles_view(struct psc_mparticles *mprts)
{
  MPI_Comm comm = psc_mparticles_comm(mprts);
  mpi_printf(comm, "  n_patches = %d\n", mprts->nr_patches);

  int n_prts_by_patch[mprts->nr_patches];
  psc_mparticles_n_prts_all(mprts, n_prts_by_patch);

  for (int p = 0; p < mprts->nr_patches; p++) {
    mpi_printf(comm, "  p %d : n_prts = %d\n", p, n_prts_by_patch[p]);
  }  
}

static void
_psc_mparticles_setup(struct psc_mparticles *mprts)
{
  mprts->mpatch = calloc(mprts->nr_patches, sizeof(*mprts->mpatch));
}

static void
_psc_mparticles_destroy(struct psc_mparticles *mparticles)
{
  free(mparticles->mpatch);
}

static void
_psc_mparticles_write(struct psc_mparticles *mparticles, struct mrc_io *io)
{
#ifdef USE_CUDA
  if (strcmp(psc_mparticles_type(mparticles), "cuda") == 0) { // FIXME
    extern void psc_mparticles_cuda_reorder(struct psc_mparticles *);
    psc_mparticles_cuda_reorder(mparticles);
  }
#endif

  mrc_io_write_ref(io, mparticles, "domain", mparticles->domain);
  mrc_io_write_int(io, mparticles, "flags", mparticles->flags);
}

static void
_psc_mparticles_read(struct psc_mparticles *mparticles, struct mrc_io *io)
{
  mparticles->domain = mrc_io_read_ref(io, mparticles, "domain", mrc_domain);
  mrc_domain_get_patches(mparticles->domain, &mparticles->nr_patches);
  mrc_io_read_int(io, mparticles, "flags", (int *) &mparticles->flags);

  _psc_mparticles_setup(mparticles);
}

int
psc_mparticles_nr_particles(struct psc_mparticles *mprts)
{
  struct psc_mparticles_ops *ops = psc_mparticles_ops(mprts);

  if (ops->get_nr_particles) {
    return ops->get_nr_particles(mprts);
  }

  int n_prts_by_patch[mprts->nr_patches];
  psc_mparticles_n_prts_all(mprts, n_prts_by_patch);
  int n_part = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    n_part += n_prts_by_patch[p];
  }
  return n_part;
}

void
psc_mparticles_n_prts_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct psc_mparticles_ops *ops = psc_mparticles_ops(mprts);
  assert(ops->get_n_prts_all);

  return ops->get_n_prts_all(mprts, n_prts_by_patch);
}

void
psc_mparticles_resize_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct psc_mparticles_ops *ops = psc_mparticles_ops(mprts);
  assert(ops->resize_all);

  return ops->resize_all(mprts, n_prts_by_patch);
}

void
psc_mparticles_setup_internals(struct psc_mparticles *mprts)
{
  struct psc_mparticles_ops *ops = psc_mparticles_ops(mprts);
  if (ops->setup_internals) {
    ops->setup_internals(mprts);
  }
}

void
psc_mparticles_reserve(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct psc_mparticles_ops *ops = psc_mparticles_ops(mprts);
  assert(ops && ops->reserve);
  ops->reserve(mprts, n_prts_by_patch);
}

struct psc_mparticles *
psc_mparticles_get_as(struct psc_mparticles *mp_base, const char *type,
		      unsigned int flags)
{
  const char *type_base = psc_mparticles_type(mp_base);
  // If we're already the subtype, nothing to be done
  if (strcmp(type_base, type) == 0)
    return mp_base;

  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_get_as", 1., 0, 0);
  }
  prof_start(pr);

  //  mprintf("get_as %s -> %s\n", psc_mparticles_type(mp_base), type);
  int nr_particles_by_patch[mp_base->nr_patches];
  psc_mparticles_n_prts_all(mp_base, nr_particles_by_patch);

  struct psc_mparticles *mp =
    psc_mparticles_create(psc_mparticles_comm(mp_base));
  psc_mparticles_set_type(mp, type);
  psc_mparticles_set_domain(mp, mp_base->domain);
  psc_mparticles_set_param_int(mp, "flags", flags);
  psc_mparticles_setup(mp);
  psc_mparticles_reserve(mp, nr_particles_by_patch);

  if (!(flags & MP_DONT_COPY)) {
    psc_mparticles_copy_to_func_t copy_to, copy_from;
    char s[strlen(type) + 12]; sprintf(s, "copy_to_%s", type);
    copy_to = (psc_mparticles_copy_to_func_t) psc_mparticles_get_method(mp_base, s);
    if (!copy_to) {
      sprintf(s, "copy_from_%s", type_base);
      copy_from = (psc_mparticles_copy_from_func_t) psc_mparticles_get_method(mp, s);
    }
    if (!copy_to && !copy_from) {
      fprintf(stderr, "ERROR: no 'copy_to_%s' in psc_mparticles '%s' and "
	      "no 'copy_from_%s' in '%s'!\n",
	      type, psc_mparticles_type(mp_base), type_base, psc_mparticles_type(mp));
      assert(0);
    }
    
#ifdef USE_CUDA
    if (strcmp(type_base, "cuda") == 0) { // FIXME
      extern void psc_mparticles_cuda_reorder(struct psc_mparticles *);
      psc_mparticles_cuda_reorder(mp_base);
    }
#endif

    if (copy_to) {
      copy_to(mp_base, mp, flags);
    } else {
      copy_from(mp, mp_base, flags);
    }
    psc_mparticles_setup_internals(mp);
  }

  prof_stop(pr);
  return mp;
}

void
psc_mparticles_put_as(struct psc_mparticles *mp, struct psc_mparticles *mp_base,
		      unsigned int flags)
{
  // If we're already the subtype, nothing to be done
  const char *type = psc_mparticles_type(mp);
  const char *type_base = psc_mparticles_type(mp_base);
  if (strcmp(type_base, type) == 0)
    return;

  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_put_as", 1., 0, 0);
  }
  prof_start(pr);

  if (!(flags & MP_DONT_COPY)) {
    psc_mparticles_copy_from_func_t copy_from, copy_to;
    char s[strlen(type) + 12]; sprintf(s, "copy_from_%s", type);
    copy_from = (psc_mparticles_copy_from_func_t) psc_mparticles_get_method(mp_base, s);
    if (!copy_from) {
      sprintf(s, "copy_to_%s", type_base);
      copy_to = (psc_mparticles_copy_from_func_t) psc_mparticles_get_method(mp, s);
    }
    if (!copy_from && !copy_to) {
      fprintf(stderr, "ERROR: no 'copy_from_%s' in psc_mparticles '%s' and "
	      "no 'copy_to_%s' in '%s'!\n",
	      type, psc_mparticles_type(mp_base), type_base, psc_mparticles_type(mp));
      assert(0);
    }

    int n_prts_by_patch[mp->nr_patches];
    psc_mparticles_n_prts_all(mp, n_prts_by_patch);
    psc_mparticles_resize_all(mp_base, n_prts_by_patch);
    
    if (copy_from) {
      copy_from(mp_base, mp, MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS);
    } else {
      copy_to(mp, mp_base, MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS);
    }
    psc_mparticles_setup_internals(mp_base);
  }
  psc_mparticles_destroy(mp);

  prof_stop(pr);
}

void
psc_mparticles_check(struct psc_mparticles *mprts_base)
{
  int fail_cnt = 0;

  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, "c", 0);
  
  psc_foreach_patch(ppsc, p) {
    struct psc_patch *patch = &ppsc->patch[p];
    particle_range_t prts = particle_range_mprts(mprts, p);

    f_real xb[3], xe[3];
    
    // New-style boundary requirements.
    // These will need revisiting when it comes to non-periodic domains.
    
    for (int d = 0; d < 3; d++) {
      xb[d] = patch->xb[d];
      xe[d] = patch->xb[d] + patch->ldims[d] * patch->dx[d];
    }

    PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
      particle_t *part = particle_iter_deref(prt_iter);
      if (part->xi < 0.f || part->xi >= xe[0] - xb[0] || // FIXME xz only!
	  part->zi < 0.f || part->zi >= xe[2] - xb[2]) {
	if (fail_cnt++ < 10) {
	  mprintf("FAIL: xi %g [%g:%g]\n", part->xi, 0., xe[0] - xb[0]);
	  mprintf("      zi %g [%g:%g]\n", part->zi, 0., xe[2] - xb[2]);
	}
      }
    }
  }
  assert(fail_cnt == 0);

  psc_mparticles_put_as(mprts, mprts_base, 0);
}

// ======================================================================

static void
psc_mparticles_init()
{
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_single_by_block_ops);
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_fortran_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_cuda_ops);
#endif
#ifdef USE_CUDA2
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_cuda2_ops);
#endif
#ifdef USE_ACC
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_acc_ops);
#endif
}

#define VAR(x) (void *)offsetof(struct psc_mparticles, x)
static struct param psc_mparticles_descr[] = {
  { "flags"             , VAR(flags)           , PARAM_INT(0)       },
  {},
};
#undef VAR

struct mrc_class_psc_mparticles mrc_class_psc_mparticles = {
  .name             = "psc_mparticles",
  .size             = sizeof(struct psc_mparticles),
  .param_descr      = psc_mparticles_descr,
  .init             = psc_mparticles_init,
  .view             = _psc_mparticles_view,
  .setup            = _psc_mparticles_setup,
  .destroy          = _psc_mparticles_destroy,
  .read             = _psc_mparticles_read,
  .write            = _psc_mparticles_write,
};

