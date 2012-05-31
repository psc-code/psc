
#include "psc_bnd_fields_private.h"

#include <string.h>

struct sub {
  struct psc_bnd_fields *fwd;
};

#define sub(bnd) mrc_to_subobj(bnd, struct sub)

static void
psc_bnd_fields_auto_setup(struct psc_bnd_fields *bnd)
{
  struct sub *sub = sub(bnd);

  const char *mflds_type = psc_mfields_type(ppsc->flds);
  char s[10 + strlen(mflds_type)];
  sprintf(s, "%s", mflds_type);

  MPI_Comm comm = psc_bnd_fields_comm(bnd);
  mpi_printf(comm, "INFO: using psc_bnd_fields '%s'\n", s);
  
  sub->fwd = psc_bnd_fields_create(comm);
  psc_bnd_fields_set_type(sub->fwd, s);
  psc_bnd_fields_setup(sub->fwd);
  psc_bnd_fields_add_child(bnd, (struct mrc_obj *) sub->fwd);

  psc_bnd_fields_setup_super(bnd);
}

static void
psc_bnd_fields_auto_fill_ghosts_a_E(struct psc_bnd_fields *bnd,
				    struct psc_fields *flds_base)
{
  struct sub *sub = sub(bnd);
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(sub->fwd);
  assert(ops->fill_ghosts_a_E);
  ops->fill_ghosts_a_E(sub->fwd, flds_base);
}

static void
psc_bnd_fields_auto_fill_ghosts_a_H(struct psc_bnd_fields *bnd,
				    struct psc_fields *flds_base)
{
  struct sub *sub = sub(bnd);
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(sub->fwd);
  assert(ops->fill_ghosts_a_H);
  ops->fill_ghosts_a_H(sub->fwd, flds_base);
}

static void
psc_bnd_fields_auto_fill_ghosts_b_E(struct psc_bnd_fields *bnd,
				    struct psc_fields *flds_base)
{
  struct sub *sub = sub(bnd);
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(sub->fwd);
  assert(ops->fill_ghosts_b_E);
  ops->fill_ghosts_b_E(sub->fwd, flds_base);
}

static void
psc_bnd_fields_auto_fill_ghosts_b_H(struct psc_bnd_fields *bnd,
				    struct psc_fields *flds_base)
{
  struct sub *sub = sub(bnd);
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(sub->fwd);
  assert(ops->fill_ghosts_b_H);
  ops->fill_ghosts_b_H(sub->fwd, flds_base);
}

static void
psc_bnd_fields_auto_add_ghosts_J(struct psc_bnd_fields *bnd,
				 struct psc_fields *flds_base)
{
  struct sub *sub = sub(bnd);
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(sub->fwd);
  assert(ops->add_ghosts_J);
  ops->add_ghosts_J(sub->fwd, flds_base);
}

// ======================================================================
// psc_bnd_fields: subclass "auto"

struct psc_bnd_fields_ops psc_bnd_fields_auto_ops = {
  .name                  = "auto",
  .size                  = sizeof(struct sub),
  .setup                 = psc_bnd_fields_auto_setup,
  .fill_ghosts_a_E       = psc_bnd_fields_auto_fill_ghosts_a_E,
  .fill_ghosts_a_H       = psc_bnd_fields_auto_fill_ghosts_a_H,
  .fill_ghosts_b_E       = psc_bnd_fields_auto_fill_ghosts_b_E,
  .fill_ghosts_b_H       = psc_bnd_fields_auto_fill_ghosts_b_H,
  .add_ghosts_J          = psc_bnd_fields_auto_add_ghosts_J,
};

