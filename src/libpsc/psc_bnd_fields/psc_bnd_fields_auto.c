
#include "psc_bnd_fields_private.h"

#include <mrc_io.h>
#include <string.h>

struct psc_bnd_fields_auto {
  struct psc_bnd_fields *fwd;
};

#define psc_bnd_fields_auto(bnd) mrc_to_subobj(bnd, struct psc_bnd_fields_auto)

static void
psc_bnd_fields_auto_setup(struct psc_bnd_fields *bnd)
{
  struct psc_bnd_fields_auto *sub = psc_bnd_fields_auto(bnd);

  const char *mflds_type = psc_mfields_type(ppsc->flds);
  char s[10 + strlen(mflds_type)];
  sprintf(s, "%s", mflds_type);

  MPI_Comm comm = psc_bnd_fields_comm(bnd);
  mpi_printf(comm, "INFO: using psc_bnd_fields '%s'\n", s);

  sub->fwd = psc_bnd_fields_create(psc_bnd_fields_comm(bnd));
  psc_bnd_fields_set_type(sub->fwd, s);
  psc_bnd_fields_setup_member_objs_sub(bnd);
  psc_bnd_fields_setup_super(bnd);
}

static void
psc_bnd_fields_auto_destroy(struct psc_bnd_fields *bnd)
{
  struct psc_bnd_fields_auto *sub = psc_bnd_fields_auto(bnd);

  psc_bnd_fields_destroy(sub->fwd);
}

static void
psc_bnd_fields_auto_write(struct psc_bnd_fields *bnd, struct mrc_io *io)
{
  struct psc_bnd_fields_auto *sub = psc_bnd_fields_auto(bnd);

  psc_bnd_fields_write(bnd, io);
  mrc_io_write_ref(io, bnd, "fwd", sub->fwd);
}

static void
psc_bnd_fields_auto_read(struct psc_bnd_fields *bnd, struct mrc_io *io)
{
  struct psc_bnd_fields_auto *sub = psc_bnd_fields_auto(bnd);

  psc_bnd_fields_read_super(bnd, io);
  sub->fwd = mrc_io_read_ref(io, bnd, "fwd", psc_bnd_fields);
}

static void
psc_bnd_fields_auto_fill_ghosts_a_E(struct psc_bnd_fields *bnd,
				    struct psc_fields *flds_base)
{
  struct psc_bnd_fields_auto *sub = psc_bnd_fields_auto(bnd);
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(sub->fwd);
  assert(ops->fill_ghosts_a_E);
  ops->fill_ghosts_a_E(sub->fwd, flds_base);
}

static void
psc_bnd_fields_auto_fill_ghosts_a_H(struct psc_bnd_fields *bnd,
				    struct psc_fields *flds_base)
{
  struct psc_bnd_fields_auto *sub = psc_bnd_fields_auto(bnd);
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(sub->fwd);
  assert(ops->fill_ghosts_a_H);
  ops->fill_ghosts_a_H(sub->fwd, flds_base);
}

static void
psc_bnd_fields_auto_fill_ghosts_b_E(struct psc_bnd_fields *bnd,
				    struct psc_fields *flds_base)
{
  struct psc_bnd_fields_auto *sub = psc_bnd_fields_auto(bnd);
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(sub->fwd);
  assert(ops->fill_ghosts_b_E);
  ops->fill_ghosts_b_E(sub->fwd, flds_base);
}

static void
psc_bnd_fields_auto_fill_ghosts_b_H(struct psc_bnd_fields *bnd,
				    struct psc_fields *flds_base)
{
  struct psc_bnd_fields_auto *sub = psc_bnd_fields_auto(bnd);
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(sub->fwd);
  assert(ops->fill_ghosts_b_H);
  ops->fill_ghosts_b_H(sub->fwd, flds_base);
}

static void
psc_bnd_fields_auto_add_ghosts_J(struct psc_bnd_fields *bnd,
				 struct psc_fields *flds_base)
{
  struct psc_bnd_fields_auto *sub = psc_bnd_fields_auto(bnd);
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(sub->fwd);
  assert(ops->add_ghosts_J);
  ops->add_ghosts_J(sub->fwd, flds_base);
}

// ======================================================================
// psc_bnd_fields: subclass "auto"

struct psc_bnd_fields_ops psc_bnd_fields_auto_ops = {
  .name                  = "auto",
  .size                  = sizeof(struct psc_bnd_fields_auto),
  .setup                 = psc_bnd_fields_auto_setup,
  .destroy               = psc_bnd_fields_auto_destroy,
  .write                 = psc_bnd_fields_auto_write,
  .read                  = psc_bnd_fields_auto_read,
  .fill_ghosts_a_E       = psc_bnd_fields_auto_fill_ghosts_a_E,
  .fill_ghosts_a_H       = psc_bnd_fields_auto_fill_ghosts_a_H,
  .fill_ghosts_b_E       = psc_bnd_fields_auto_fill_ghosts_b_E,
  .fill_ghosts_b_H       = psc_bnd_fields_auto_fill_ghosts_b_H,
  .add_ghosts_J          = psc_bnd_fields_auto_add_ghosts_J,
};

