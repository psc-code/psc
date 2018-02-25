
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
_psc_mparticles_view(struct psc_mparticles *_mprts)
{
  MPI_Comm comm = psc_mparticles_comm(_mprts);
  mparticles_base_t mprts(_mprts);
  mpi_printf(comm, "  n_patches    = %d\n", mprts->n_patches());
  mpi_printf(comm, "  n_prts_total = %d\n", mprts->get_n_prts());

  uint n_prts_by_patch[mprts->n_patches()];
  mprts->get_size_all(n_prts_by_patch);

  for (int p = 0; p < mprts->n_patches(); p++) {
    mpi_printf(comm, "  p %d: n_prts = %d\n", p, n_prts_by_patch[p]);
  }  
}

static void
copy(struct psc_mparticles *mprts_from, struct psc_mparticles *mprts_to,
     const char *type_from, const char *type_to,
     unsigned int flags)
{
  mparticles_base_t mp_from(mprts_from), mp_to(mprts_to);
  psc_mparticles_copy_func_t copy_to, copy_from;

  assert(mp_from->n_patches() == mp_to->n_patches());

  if (flags & MP_DONT_COPY) {
    if (!(flags & MP_DONT_RESIZE)) {
      uint n_prts_by_patch[mp_from->n_patches()];
      mp_from->get_size_all(n_prts_by_patch);
      mp_to->reserve_all(n_prts_by_patch);
      mp_to->resize_all(n_prts_by_patch);
    }
    return;
  }

  assert(!(flags & MP_DONT_RESIZE));

  char s[std::max(strlen(type_from),strlen(type_to)) + 12];
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
}

struct psc_mparticles *
psc_mparticles_get_as(struct psc_mparticles *mprts_from, const char *type,
		      unsigned int flags)
{
  mparticles_base_t mp_from(mprts_from);
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
  mprts->grid = mprts_from->grid;
  psc_mparticles_setup(mprts);

  copy(mprts_from, mprts, type_from, type, flags);

  //  mprintf("get_as %s -> %s to\n", type_from, type);
  //  psc_mparticles_view(mprts);

  prof_stop(pr);
  return mprts;
}

void
psc_mparticles_put_as(struct psc_mparticles *mprts, struct psc_mparticles *mprts_to,
		      unsigned int flags)
{
  mparticles_base_t mp(mprts), mp_to(mprts_to);
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
  
  if (flags & MP_DONT_COPY) {
    // let's check that the size of the particle arrays hasn't changed, since
    // it's not obvious what we should do in case it did...
    uint n_prts_by_patch[mp->n_patches()];
    uint n_prts_by_patch_to[mp_to->n_patches()];

    mp->get_size_all(n_prts_by_patch);
    mp_to->get_size_all(n_prts_by_patch_to);
    assert(mp_to->n_patches() == mp->n_patches());

    for (int p = 0; p < mp->n_patches(); p++) {
      if (n_prts_by_patch[p] != n_prts_by_patch_to[p]) {
	mprintf("psc_mparticles_put_as: p = %d n_prts %d -- %d\n",
		p, n_prts_by_patch[p], n_prts_by_patch_to[p]);
      }
      assert(n_prts_by_patch[p] == n_prts_by_patch_to[p]);
    }

    flags |= MP_DONT_RESIZE;
  }
  
  copy(mprts, mprts_to, type, type_to, flags);
  
  psc_mparticles_destroy(mprts);

  //  mprintf("put_as %s -> %s to\n", type, type_to);
  //  psc_mparticles_view(mprts_to);
  prof_stop(pr);
}

void
psc_mparticles_check(struct psc_mparticles *mprts_base)
{
  int fail_cnt = 0;

  mparticles_t mprts = mprts_base->get_as<PscMparticlesDouble>(0);
  const Grid_t& grid = ppsc->grid();
  
  psc_foreach_patch(ppsc, p) {
    auto& patch = ppsc->grid().patches[p];
    mparticles_t::patch_t& prts = mparticles_t(mprts)[p];

    f_real xb[3], xe[3];
    
    // New-style boundary requirements.
    // These will need revisiting when it comes to non-periodic domains.
    
    for (int d = 0; d < 3; d++) {
      xb[d] = patch.xb[d];
      xe[d] = patch.xb[d] + grid.ldims[d] * grid.dx[d];
    }

    PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
      particle_t *part = &*prt_iter;
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

  mprts.put_as(mprts_base, 0);
}

// ======================================================================

extern struct psc_mparticles_ops psc_mparticles_single_ops;
extern struct psc_mparticles_ops psc_mparticles_double_ops;
extern struct psc_mparticles_ops psc_mparticles_sse2_ops;
extern struct psc_mparticles_ops psc_mparticles_cbe_ops;
extern struct psc_mparticles_ops psc_mparticles_cuda_ops;
extern struct psc_mparticles_ops psc_mparticles_cuda2_ops;
extern struct psc_mparticles_ops psc_mparticles_acc_ops;
extern struct psc_mparticles_ops psc_mparticles_vpic_ops;
extern struct psc_mparticles_ops psc_mparticles_single_by_kind_ops;

static void
psc_mparticles_init()
{
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_single_by_kind_ops);
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

struct mrc_class_psc_mparticles_ : mrc_class_psc_mparticles {
  mrc_class_psc_mparticles_() {
    name             = "psc_mparticles";
    size             = sizeof(struct psc_mparticles);
    init             = psc_mparticles_init;
    view             = _psc_mparticles_view;
  }
} mrc_class_psc_mparticles;

