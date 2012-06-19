
#include "psc_output_fields_item_private.h"

#include <string.h>

struct sub {
  struct psc_output_fields_item *fwd;
};

#define sub(bnd) mrc_to_subobj(bnd, struct sub)

// ======================================================================
// generic fwd

static const char *
fwd_get_component_name(struct psc_output_fields_item *item, int m)
{
  struct sub *sub = sub(item);
  struct psc_output_fields_item_ops *ops = psc_output_fields_item_ops(sub->fwd);
  if (ops->get_component_name) {
    return ops->get_component_name(sub->fwd, m);
  }

  return ops->fld_names[m];
}

static int
fwd_get_nr_components(struct psc_output_fields_item *item)
{
  struct sub *sub = sub(item);
  struct psc_output_fields_item_ops *ops = psc_output_fields_item_ops(sub->fwd);
  if (ops->get_nr_components) {
    return ops->get_nr_components(sub->fwd);
  }

  return ops->nr_comp;
}

static void
fwd_run(struct psc_output_fields_item *item, struct psc_fields *flds,
	struct psc_particles *prts, struct psc_fields *res)
{
  struct sub *sub = sub(item);
  struct psc_output_fields_item_ops *ops = psc_output_fields_item_ops(sub->fwd);
  assert(ops->run);
  ops->run(item, flds, prts, res);
}

// ======================================================================

#define MAKE_FWD(what, flgs)						\
static void								\
what##_setup(struct psc_output_fields_item *item)			\
{									\
  struct sub *sub = sub(item);						\
									\
  const char *mprts_type = psc_mparticles_type(ppsc->particles);	\
  char s[10 + strlen(mprts_type)];					\
  if (strcmp(mprts_type, "cuda") == 0) {				\
    sprintf(s, #what "_single");					\
  } else {								\
    sprintf(s, #what "_%s", mprts_type);				\
  }									\
  									\
  MPI_Comm comm = psc_output_fields_item_comm(item);			\
  mpi_printf(comm, "INFO: using psc_output_fields_item '%s'\n", s);	\
  									\
  sub->fwd = psc_output_fields_item_create(comm);			\
  psc_output_fields_item_set_type(sub->fwd, s);				\
  psc_output_fields_item_set_psc_bnd(sub->fwd, item->bnd);		\
  psc_output_fields_item_setup(sub->fwd);				\
  psc_output_fields_item_add_child(item, (struct mrc_obj *) sub->fwd);	\
}									\
									\
struct psc_output_fields_item_ops psc_output_fields_item_##what##_ops = { \
  .name               = #what,						\
  .size               = sizeof(struct sub),				\
  .setup              = what##_setup,					\
  .get_nr_components  = fwd_get_nr_components,				\
  .get_component_name = fwd_get_component_name,				\
  .run                = fwd_run,					\
  .flags              = flgs,						\
};

MAKE_FWD(n_1st, POFI_ADD_GHOSTS | POFI_BY_KIND) // FIXME, flags duplicated
MAKE_FWD(v_1st, POFI_ADD_GHOSTS | POFI_BY_KIND)
MAKE_FWD(vv_1st, POFI_ADD_GHOSTS | POFI_BY_KIND)

