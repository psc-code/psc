
#include "psc_bnd_private.h"
#include "psc_bnd_fld.h"

#include <mrc_ddc.h>
#include <string.h>

struct sub {
  struct psc_bnd *fwd;
};

#define sub(o) mrc_to_subobj(o, struct sub)

// ----------------------------------------------------------------------
// psc_bnd_sub_setup

static void
psc_bnd_sub_setup(struct psc_bnd *bnd)
{
  struct sub *sub = sub(bnd);

  const char *mprts_type = psc_mparticles_type(ppsc->particles);
  char s[10 + strlen(mprts_type)];
  if (0) {
  } else {
    strcpy(s, mprts_type);
  }

  MPI_Comm comm = psc_bnd_comm(bnd);
  mpi_printf(comm, "INFO: using psc_bnd '%s'\n", s);
  
  sub->fwd = psc_bnd_create(comm);
  psc_bnd_set_type(sub->fwd, s);
  psc_bnd_set_psc(sub->fwd, bnd->psc);
  psc_bnd_setup(sub->fwd);
  psc_bnd_add_child(bnd, (struct mrc_obj *) sub->fwd);

  psc_bnd_setup_super(bnd);
}

// ----------------------------------------------------------------------
// psc_bnd_sub_fields_create

static void
psc_bnd_sub_fields_create(struct psc_bnd *bnd)
{
  struct sub *sub = sub(bnd);
  bnd->ddc = mrc_ddc_get(sub->fwd->ddc);
  assert(bnd->ddc);
}

// ----------------------------------------------------------------------
// psc_bnd_sub_fields_add_ghosts

static void
psc_bnd_sub_fields_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds,
			      int mb, int me)
{
  struct sub *sub = sub(bnd);
  struct psc_bnd_ops *ops = psc_bnd_ops(sub->fwd);
  assert(ops->add_ghosts);
  ops->add_ghosts(sub->fwd, mflds, mb, me);
}

// ----------------------------------------------------------------------
// psc_bnd_sub_fields_fill_ghosts

static void
psc_bnd_sub_fields_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds,
			      int mb, int me)
{
  struct sub *sub = sub(bnd);
  struct psc_bnd_ops *ops = psc_bnd_ops(sub->fwd);
  assert(ops->fill_ghosts);
  ops->fill_ghosts(sub->fwd, mflds, mb, me);
}

// ======================================================================
// psc_bnd: subclass "auto"

struct psc_bnd_ops psc_bnd_auto_ops = {
  .name                    = "auto",
  .size                    = sizeof(struct sub),
  .setup                   = psc_bnd_sub_setup,
  .create_ddc              = psc_bnd_sub_fields_create,
  .add_ghosts              = psc_bnd_sub_fields_add_ghosts,
  .fill_ghosts             = psc_bnd_sub_fields_fill_ghosts,
};
