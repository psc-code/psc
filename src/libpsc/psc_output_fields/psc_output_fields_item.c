
#include "psc_output_fields_item_private.h"

// ======================================================================
// psc_output_fields_item_init

static void
psc_output_fields_item_init()
{
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_j_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_j_nc_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_j_ec_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_e_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_e_nc_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_e_ec_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_h_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_h_nc_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_h_fc_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_jdote_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_poyn_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_e2_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_h2_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_densities_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_v_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_vv_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_photon_n_ops);
}

// ======================================================================
// psc_output_fields_item class

struct mrc_class_psc_output_fields_item mrc_class_psc_output_fields_item = {
  .name             = "psc_output_fields_item",
  .size             = sizeof(struct psc_output_fields_item),
  .init             = psc_output_fields_item_init,
};

