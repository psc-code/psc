
#include "psc_bnd_particles_private.h"

#include <mrc_ddc.h>
#include <string.h>

struct sub {
  struct psc_bnd_particles *fwd;
};

#define sub(o) mrc_to_subobj(o, struct sub)

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_setup

static void
psc_bnd_particles_sub_setup(struct psc_bnd_particles *bnd)
{
  struct sub *sub = sub(bnd);

  const char *mprts_type = psc_mparticles_type(ppsc->particles);
  char s[10 + strlen(mprts_type)];
  if (0) {
  } else {
    strcpy(s, mprts_type);
  }

  MPI_Comm comm = psc_bnd_particles_comm(bnd);
  mpi_printf(comm, "INFO: using psc_bnd_particles '%s'\n", s);
  
  sub->fwd = psc_bnd_particles_create(comm);
  psc_bnd_particles_set_type(sub->fwd, s);
  psc_bnd_particles_set_psc(sub->fwd, bnd->psc);
  psc_bnd_particles_setup(sub->fwd);
  psc_bnd_particles_add_child(bnd, (struct mrc_obj *) sub->fwd);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_unsetup

static void
psc_bnd_particles_sub_unsetup(struct psc_bnd_particles *bnd)
{
  // may not need to do anything, since the sub does check_domain() itself
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles

static void
psc_bnd_particles_sub_exchange_particles(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts)
{
  // FIXME, this skips check_domain() in the forwarded class
  struct sub *sub = sub(bnd);
  struct psc_bnd_particles_ops *ops = psc_bnd_particles_ops(sub->fwd);
  assert(ops->exchange_particles);
  ops->exchange_particles(sub->fwd, mprts);
}

// ======================================================================
// psc_bnd_particles: subclass "auto"

struct psc_bnd_particles_ops psc_bnd_particles_auto_ops = {
  .name                    = "auto",
  .size                    = sizeof(struct sub),
  .setup                   = psc_bnd_particles_sub_setup,
  .unsetup                 = psc_bnd_particles_sub_unsetup,
  .exchange_particles      = psc_bnd_particles_sub_exchange_particles,
};
