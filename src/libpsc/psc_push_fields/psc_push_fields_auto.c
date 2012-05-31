
#include "psc_push_fields_private.h"

#include <string.h>

struct sub {
  struct psc_push_fields *fwd;
};

#define sub(push) mrc_to_subobj(push, struct sub)

static void
psc_push_fields_auto_setup(struct psc_push_fields *push)
{
  struct sub *sub = sub(push);

  const char *mflds_type = psc_mfields_type(ppsc->flds);
  char s[10 + strlen(mflds_type)];
  sprintf(s, "%s", mflds_type);

  MPI_Comm comm = psc_push_fields_comm(push);
  mpi_printf(comm, "INFO: using psc_push_fields '%s'\n", s);
  
  sub->fwd = psc_push_fields_create(comm);
  psc_push_fields_set_type(sub->fwd, s);
  psc_push_fields_setup(sub->fwd);
  psc_push_fields_add_child(push, (struct mrc_obj *) sub->fwd);

  psc_push_fields_setup_super(push);
}

static void
psc_push_fields_auto_push_a_E(struct psc_push_fields *push,
			      struct psc_fields *flds_base)
{
  struct sub *sub = sub(push);
  struct psc_push_fields_ops *ops = psc_push_fields_ops(sub->fwd);
  assert(ops->push_a_E);
  ops->push_a_E(sub->fwd, flds_base);
}

static void
psc_push_fields_auto_push_a_H(struct psc_push_fields *push,
			      struct psc_fields *flds_base)
{
  struct sub *sub = sub(push);
  struct psc_push_fields_ops *ops = psc_push_fields_ops(sub->fwd);
  assert(ops->push_a_H);
  ops->push_a_H(sub->fwd, flds_base);
}

static void
psc_push_fields_auto_push_b_H(struct psc_push_fields *push,
			      struct psc_fields *flds_base)
{
  struct sub *sub = sub(push);
  struct psc_push_fields_ops *ops = psc_push_fields_ops(sub->fwd);
  assert(ops->push_b_H);
  ops->push_b_H(sub->fwd, flds_base);
}

static void
psc_push_fields_auto_push_b_E(struct psc_push_fields *push,
			      struct psc_fields *flds_base)
{
  struct sub *sub = sub(push);
  struct psc_push_fields_ops *ops = psc_push_fields_ops(sub->fwd);
  assert(ops->push_b_E);
  ops->push_b_E(sub->fwd, flds_base);
}

// ======================================================================
// psc_push_fields: subclass "auto"

struct psc_push_fields_ops psc_push_fields_auto_ops = {
  .name                  = "auto",
  .size                  = sizeof(struct sub),
  .setup                 = psc_push_fields_auto_setup,
  .push_a_E              = psc_push_fields_auto_push_a_E,
  .push_a_H              = psc_push_fields_auto_push_a_H,
  .push_b_H              = psc_push_fields_auto_push_b_H,
  .push_b_E              = psc_push_fields_auto_push_b_E,
};

