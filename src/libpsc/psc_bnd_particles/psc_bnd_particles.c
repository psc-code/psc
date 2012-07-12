
#include <psc_bnd_particles_private.h>

#include "../psc_bnd/ddc_particles.h"

#include <mrc_io.h>
#include <mrc_ddc.h>
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_bnd_particles_setup

static void
_psc_bnd_particles_setup(struct psc_bnd_particles *bnd)
{
  struct mrc_ddc *ddc = mrc_domain_create_ddc(bnd->psc->mrc_domain);
  //  mrc_ddc_set_funcs(ddc, &ddc_funcs);
  mrc_ddc_set_param_int3(ddc, "ibn", bnd->psc->ibn);
  mrc_ddc_set_param_int(ddc, "max_n_fields", 24);
  mrc_ddc_set_param_int(ddc, "size_of_type", sizeof(float));
  mrc_ddc_setup(ddc);
  bnd->ddc = ddc;
}

// ----------------------------------------------------------------------
// psc_bnd_particles_destroy

static void
_psc_bnd_particles_destroy(struct psc_bnd_particles *bnd)
{
  struct psc_bnd_particles_ops *ops = psc_bnd_particles_ops(bnd);
  assert(ops->unsetup);
  ops->unsetup(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_write

static void
_psc_bnd_particles_write(struct psc_bnd_particles *bnd, struct mrc_io *io)
{
  const char *path = psc_bnd_particles_name(bnd);
  mrc_io_write_obj_ref(io, path, "psc", (struct mrc_obj *) bnd->psc);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_read

static void
_psc_bnd_particles_read(struct psc_bnd_particles *bnd, struct mrc_io *io)
{
  const char *path = psc_bnd_particles_name(bnd);
  bnd->psc = (struct psc *)
    mrc_io_read_obj_ref(io, path, "psc", &mrc_class_psc);

  psc_bnd_particles_setup(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_set_psc

void
psc_bnd_particles_set_psc(struct psc_bnd_particles *bnd, struct psc *psc)
{
  bnd->psc = psc;
}

// ----------------------------------------------------------------------
// check_domain
//
// check if the underlying mrc_domain changed since setup(),
// which might happen, e.g., through rebalancing.
// In this case, do setup() over.

static void
check_domain(struct psc_bnd_particles *bnd)
{
  struct psc_bnd_particles_ops *ops = psc_bnd_particles_ops(bnd);

  struct mrc_domain *domain = mrc_ddc_get_domain(bnd->ddc);
  if (domain != bnd->psc->mrc_domain) {
    ops->unsetup(bnd);
    mrc_ddc_destroy(bnd->ddc);

    ops->setup(bnd);
  }
}

// ----------------------------------------------------------------------
// psc_bnd_particles_exchange

void
psc_bnd_particles_exchange(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts)
{
  static int pr;
  if (!pr) {
    pr = prof_register("xchg_prts", 1., 0, 0);
  }

  check_domain(bnd);

  prof_start(pr);
  psc_stats_start(st_time_comm);
  struct psc_bnd_particles_ops *ops = psc_bnd_particles_ops(bnd);
  assert(ops->exchange_particles);
  ops->exchange_particles(bnd, mprts);
  psc_stats_stop(st_time_comm);
  prof_stop(pr);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_init

static void
psc_bnd_particles_init()
{
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_auto_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_single2_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd_particles, &psc_bnd_particles_mix_ops);
#endif
}

// ======================================================================
// psc_bnd_particles class

struct mrc_class_psc_bnd_particles mrc_class_psc_bnd_particles = {
  .name             = "psc_bnd_particles",
  .size             = sizeof(struct psc_bnd_particles),
  .init             = psc_bnd_particles_init,
  .setup            = _psc_bnd_particles_setup,
  .destroy          = _psc_bnd_particles_destroy,
  .write            = _psc_bnd_particles_write,
  .read             = _psc_bnd_particles_read,
};

