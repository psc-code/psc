
#include "psc.h"

#include "psc_particles_as_double.h" // only for checking...

#include <mrc_profile.h>
#include <mrc_params.h>
#include <mrc_io.h>
#include <stdlib.h>
#include <string.h>

// ======================================================================
// psc_mparticles base class

static void
_psc_mparticles_view(struct psc_mparticles *mprts)
{
  MPI_Comm comm = psc_mparticles_comm(mprts);
  mpi_printf(comm, "  n_patches    = %d\n", mprts->nr_patches);
  mpi_printf(comm, "  n_prts_total = %d\n", psc_mparticles_nr_particles(mprts));

  int n_prts_by_patch[mprts->nr_patches];
  psc_mparticles_get_size_all(mprts, n_prts_by_patch);

  for (int p = 0; p < mprts->nr_patches; p++) {
    mpi_printf(comm, "  p %d: n_prts = %d\n", p, n_prts_by_patch[p]);
  }  
}

int
psc_mparticles_nr_particles(struct psc_mparticles *mprts)
{
  struct psc_mparticles_ops *ops = psc_mparticles_ops(mprts);
  assert(ops->get_nr_particles);

  return ops->get_nr_particles(mprts);
}

void
psc_mparticles_get_size_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct psc_mparticles_ops *ops = psc_mparticles_ops(mprts);
  assert(ops->get_size_all);

  return ops->get_size_all(mprts, n_prts_by_patch);
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
psc_mparticles_reserve_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct psc_mparticles_ops *ops = psc_mparticles_ops(mprts);
  assert(ops && ops->reserve_all);
  ops->reserve_all(mprts, n_prts_by_patch);
}

static void
copy(struct psc_mparticles *mprts_from, struct psc_mparticles *mprts_to,
     const char *type_from, const char *type_to,
     unsigned int flags)
{
  psc_mparticles_copy_func_t copy_to, copy_from;

  char s[strlen(type_to) + 12];
  sprintf(s, "copy_to_%s", type_to);
  copy_to = (psc_mparticles_copy_func_t) psc_mparticles_get_method(mprts_from, s);
  if (!copy_to) {
    sprintf(s, "copy_from_%s", type_from);
    copy_from = (psc_mparticles_copy_func_t) psc_mparticles_get_method(mprts_to, s);
  }
  if (!copy_to && !copy_from) {
    fprintf(stderr, "ERROR: no 'copy_to_%s' in psc_mparticles '%s' and "
	    "no 'copy_from_%s' in '%s'!\n",
	    type_to, psc_mparticles_type(mprts_from), type_from, psc_mparticles_type(mprts_to));
    assert(0);
  }

  if (copy_to) {
    copy_to(mprts_from, mprts_to, flags);
  } else {
    copy_from(mprts_to, mprts_from, flags);
  }
  psc_mparticles_setup_internals(mprts_to);
}

struct psc_mparticles *
psc_mparticles_get_as(struct psc_mparticles *mprts_from, const char *type,
		      unsigned int flags)
{
  const char *type_from = psc_mparticles_type(mprts_from);
  // If we're already the subtype, nothing to be done
  if (strcmp(type_from, type) == 0) {
    return mprts_from;
  }
  
  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_get_as", 1., 0, 0);
  }
  prof_start(pr);

  //  mprintf("get_as %s -> %s from\n", type_from, type);
  //  psc_mparticles_view(mprts_from);

  struct psc_mparticles *mprts = psc_mparticles_create(psc_mparticles_comm(mprts_from));
  psc_mparticles_set_type(mprts, type);
  psc_mparticles_set_param_int(mprts, "nr_patches", mprts_from->nr_patches);
  psc_mparticles_set_param_int(mprts, "flags", flags);
  psc_mparticles_setup(mprts);

  if (!(flags & MP_DONT_RESIZE)) {
    int nr_particles_by_patch[mprts_from->nr_patches];
    psc_mparticles_get_size_all(mprts_from, nr_particles_by_patch);
    psc_mparticles_reserve_all(mprts, nr_particles_by_patch);
    psc_mparticles_resize_all(mprts, nr_particles_by_patch);
  }
  
  if (!(flags & MP_DONT_COPY)) {
    copy(mprts_from, mprts, type_from, type, flags);
  }

  //  mprintf("get_as %s -> %s to\n", type_from, type);
  //  psc_mparticles_view(mprts);

  prof_stop(pr);
  return mprts;
}

void
psc_mparticles_put_as(struct psc_mparticles *mprts, struct psc_mparticles *mprts_to,
		      unsigned int flags)
{
  // If we're already the subtype, nothing to be done
  const char *type = psc_mparticles_type(mprts);
  const char *type_to = psc_mparticles_type(mprts_to);
  if (strcmp(type_to, type) == 0) {
    return;
  }

  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_put_as", 1., 0, 0);
  }
  prof_start(pr);

  //  mprintf("put_as %s -> %s from\n", type, type_to);
  //  psc_mparticles_view(mprts);
  
  int n_prts_by_patch[mprts->nr_patches];
  psc_mparticles_get_size_all(mprts, n_prts_by_patch);
  if (!(flags & MP_DONT_COPY)) {
    psc_mparticles_reserve_all(mprts_to, n_prts_by_patch);
    psc_mparticles_resize_all(mprts_to, n_prts_by_patch);
    copy(mprts, mprts_to, type, type_to, flags);
  } else {
    // let's check that the size of the particle arrays hasn't changed, since
    // it's not obvious what we should do in case it did...
    int n_prts_by_patch_to[mprts_to->nr_patches];
    psc_mparticles_get_size_all(mprts_to, n_prts_by_patch_to);
    assert(mprts_to->nr_patches == mprts->nr_patches);
    for (int p = 0; p < mprts->nr_patches; p++) {
      if (n_prts_by_patch[p] != n_prts_by_patch_to[p]) {
	mprintf("psc_mparticles_put_as: p = %d n_prts %d -- %d\n",
		p, n_prts_by_patch[p], n_prts_by_patch_to[p]);
      }
      assert(n_prts_by_patch[p] == n_prts_by_patch_to[p]);
    }
  }
  psc_mparticles_destroy(mprts);

  //  mprintf("put_as %s -> %s to\n", type, type_to);
  //  psc_mparticles_view(mprts_to);
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
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_vpic_ops);
#endif
}

#define VAR(x) (void *)offsetof(struct psc_mparticles, x)
static struct param psc_mparticles_descr[] = {
  { "nr_patches"        , VAR(nr_patches)      , PARAM_INT(0)       },
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
};

