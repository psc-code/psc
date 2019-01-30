
#include "psc_diag_item_private.h"

int
psc_diag_item_nr_values(struct psc_diag_item *item)
{
  struct psc_diag_item_ops *ops = psc_diag_item_ops(item);

  return ops->nr_values;
}

const char *
psc_diag_item_title(struct psc_diag_item *item, int i)
{
  struct psc_diag_item_ops *ops = psc_diag_item_ops(item);

  return ops->title[i];
}

void
psc_diag_item_run(struct psc_diag_item *item,
		  MparticlesBase& mprts, MfieldsStateBase& mflds,
		  double *result)
{
  struct psc_diag_item_ops *ops = psc_diag_item_ops(item);

  ops->run(item, mprts, mflds, result);
}

// ----------------------------------------------------------------------
// psc_diag_item_init

extern struct psc_diag_item_ops psc_diag_item_field_energy_ops;
extern struct psc_diag_item_ops psc_diag_item_particle_energy_ops;

static void
psc_diag_item_init()
{
  mrc_class_register_subclass(&mrc_class_psc_diag_item,
			      &psc_diag_item_field_energy_ops);
  mrc_class_register_subclass(&mrc_class_psc_diag_item,
			      &psc_diag_item_particle_energy_ops);
}

// ----------------------------------------------------------------------
// psc_diag_item class

struct mrc_class_psc_diag_item_ : mrc_class_psc_diag_item {
  mrc_class_psc_diag_item_() {
    name             = "psc_diag_item";
    size             = sizeof(struct psc_diag_item);
    init             = psc_diag_item_init;
  }
} mrc_class_psc_diag_item;

