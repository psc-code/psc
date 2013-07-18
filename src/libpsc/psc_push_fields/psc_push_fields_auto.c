
#include "psc_push_fields_private.h"
#include "psc_bnd_fields_private.h"

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

  psc_bnd_fields_destroy(push->bnd_fields);
  push->bnd_fields = psc_bnd_fields_get(sub->fwd->bnd_fields);
}

static void
psc_push_fields_auto_push_E(struct psc_push_fields *push,
			    struct psc_fields *flds)
{
  struct sub *sub = sub(push);
  struct psc_push_fields_ops *ops = psc_push_fields_ops(sub->fwd);
  assert(ops->push_E);
  ops->push_E(sub->fwd, flds);
}

static void
psc_push_fields_auto_push_H(struct psc_push_fields *push,
			    struct psc_fields *flds)
{
  struct sub *sub = sub(push);
  struct psc_push_fields_ops *ops = psc_push_fields_ops(sub->fwd);
  assert(ops->push_H);
  ops->push_H(sub->fwd, flds);
}

// ======================================================================
// psc_push_fields: subclass "auto"

struct psc_push_fields_ops psc_push_fields_auto_ops = {
  .name                  = "auto",
  .size                  = sizeof(struct sub),
  .setup                 = psc_push_fields_auto_setup,
  .push_E                = psc_push_fields_auto_push_E,
  .push_H                = psc_push_fields_auto_push_H,
};

