
#include "psc_bnd_private.h"

void
psc_bnd_set_psc(struct psc_bnd *bnd, struct psc *psc)
{
  bnd->psc = psc;
}

// ======================================================================
// forward to subclass

void
psc_bnd_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds, int mb, int me)
{
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->add_ghosts);
  ops->add_ghosts(bnd, flds, mb, me);
}

void
psc_bnd_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds, int mb, int me)
{
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->fill_ghosts);
  ops->fill_ghosts(bnd, flds, mb, me);
}

void
psc_bnd_exchange_particles(struct psc_bnd *bnd, mparticles_base_t *particles)
{
  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->exchange_particles);
  ops->exchange_particles(bnd, particles);
}

void
psc_bnd_exchange_photons(struct psc_bnd *bnd, mphotons_t *mphotons)
{
  int n_total = 0;
  psc_foreach_patch(bnd->psc, p) {
    n_total += mphotons->p[p].nr;
  }
  if (n_total == 0)
    return;

  struct psc_bnd_ops *ops = psc_bnd_ops(bnd);
  assert(ops->exchange_photons);
  ops->exchange_photons(bnd, mphotons);
}

// ======================================================================
// psc_bnd_init

static void
psc_bnd_init()
{
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_fortran_ops);
}

// ======================================================================
// psc_bnd class

struct mrc_class_psc_bnd mrc_class_psc_bnd = {
  .name             = "psc_bnd",
  .size             = sizeof(struct psc_bnd),
  .init             = psc_bnd_init,
};

