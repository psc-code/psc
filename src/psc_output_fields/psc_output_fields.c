
#include "psc_output_fields_private.h"

// ======================================================================
// forward to subclass

void
psc_output_fields_run(struct psc_output_fields *output_fields,
		      mfields_base_t *flds, mparticles_base_t *particles)
{
  struct psc_output_fields_ops *ops = psc_output_fields_ops(output_fields);
  assert(ops->run);
  ops->run(output_fields, flds, particles);
}

// ======================================================================
// psc_output_fields_init

static void
psc_output_fields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_output_fields, &psc_output_fields_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields, &psc_output_fields_fortran_ops);
}

// ======================================================================
// psc_output_fields class

struct mrc_class_psc_output_fields mrc_class_psc_output_fields = {
  .name             = "psc_output_fields",
  .size             = sizeof(struct psc_output_fields),
  .init             = psc_output_fields_init,
};

