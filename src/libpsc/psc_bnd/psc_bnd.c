
#include "psc_bnd_private.h"
#include "ddc_particles.h"

#include <mrc_io.h>
#include <mrc_ddc.h>
#include <mrc_profile.h>

// ======================================================================
// psc_bnd

// ----------------------------------------------------------------------
// psc_set_psc

void
psc_bnd_set_psc(struct psc_bnd *bnd, struct psc *psc)
{
  bnd->psc = psc;
}

// ----------------------------------------------------------------------
// psc_bnd_setup

static void
_psc_bnd_setup(struct psc_bnd *bnd)
{
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->create_ddc);
  ops->create_ddc(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_destroy

static void
_psc_bnd_destroy(struct psc_bnd *bnd)
{
  mrc_ddc_destroy(bnd->ddc);
}

// ----------------------------------------------------------------------
// psc_bnd_write

static void
_psc_bnd_write(struct psc_bnd *bnd, struct mrc_io *io)
{
  const char *path = psc_bnd_name(bnd);
  mrc_io_write_obj_ref(io, path, "psc", (struct mrc_obj *) bnd->psc);
}

// ----------------------------------------------------------------------
// psc_bnd_read

static void
_psc_bnd_read(struct psc_bnd *bnd, struct mrc_io *io)
{
  const char *path = psc_bnd_name(bnd);
  bnd->psc = (struct psc *)
    mrc_io_read_obj_ref(io, path, "psc", &mrc_class_psc);

  psc_bnd_setup(bnd);
}

// ----------------------------------------------------------------------
// check_domain
//
// check if the underlying mrc_domain changed since setup(),
// which might happen, e.g., through rebalancing.
// In this case, do setup() over.

static void
check_domain(struct psc_bnd *bnd)
{
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);

  struct mrc_domain *domain = mrc_ddc_get_domain(bnd->ddc);
  if (domain != bnd->psc->mrc_domain) {
    mrc_ddc_destroy(bnd->ddc);

    ops->create_ddc(bnd);
    ops->setup(bnd);
  }
}

// ======================================================================
// forward to subclass

void
psc_bnd_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds, int mb, int me)
{
  static int pr;
  if (!pr) {
    pr = prof_register("add_ghosts", 1., 0, 0);
  }

  check_domain(bnd);

  psc_stats_start(st_time_comm);
  prof_start(pr);

  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->add_ghosts);
  ops->add_ghosts(bnd, flds, mb, me);

  prof_stop(pr);
  psc_stats_stop(st_time_comm);
}

void
psc_bnd_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds, int mb, int me)
{
  static int pr;
  if (!pr) {
    pr = prof_register("fill_ghosts", 1., 0, 0);
  }

  check_domain(bnd);

  psc_stats_start(st_time_comm);
  prof_start(pr);

  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->fill_ghosts);
  ops->fill_ghosts(bnd, flds, mb, me);

  prof_stop(pr);
  psc_stats_stop(st_time_comm);
}

// ======================================================================
// psc_bnd_init

static void
psc_bnd_init()
{
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_auto_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_single_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_mix_ops);
#endif
}

// ======================================================================
// psc_bnd class

struct mrc_class_psc_bnd mrc_class_psc_bnd = {
  .name             = "psc_bnd",
  .size             = sizeof(struct psc_bnd),
  .init             = psc_bnd_init,
  .setup            = _psc_bnd_setup,
  .destroy          = _psc_bnd_destroy,
  .write            = _psc_bnd_write,
  .read             = _psc_bnd_read,
};

