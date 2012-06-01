
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
psc_push_fields_auto_step_a(struct psc_push_fields *push,
			    struct psc_mfields *mflds)
{
  struct sub *sub = sub(push);
  struct psc_push_fields_ops *ops = psc_push_fields_ops(sub->fwd);
  if (ops->step_a) {
    ops->step_a(sub->fwd, mflds);
  } else {
    psc_push_fields_step_a(sub->fwd, mflds);
  }
}

static void
psc_push_fields_auto_step_b_H(struct psc_push_fields *push,
			    struct psc_mfields *mflds)
{
  struct sub *sub = sub(push);
  psc_push_fields_step_b_H(sub->fwd, mflds);
}

static void
psc_push_fields_auto_step_b_E(struct psc_push_fields *push,
			      struct psc_mfields *mflds)
{
  struct sub *sub = sub(push);
  psc_push_fields_step_b_E(sub->fwd, mflds);
}

// ======================================================================
// psc_push_fields: subclass "auto"

struct psc_push_fields_ops psc_push_fields_auto_ops = {
  .name                  = "auto",
  .size                  = sizeof(struct sub),
  .setup                 = psc_push_fields_auto_setup,
  .step_a                = psc_push_fields_auto_step_a,
  .step_b_H              = psc_push_fields_auto_step_b_H,
  .step_b_E              = psc_push_fields_auto_step_b_E,
};

