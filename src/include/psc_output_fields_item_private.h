
#ifndef PSC_OUTPUT_FIELDS_ITEM_PRIVATE_H
#define PSC_OUTPUT_FIELDS_ITEM_PRIVATE_H

#include <psc_output_fields_item.h>

struct psc_output_fields_item {
  struct mrc_obj obj;
};

struct psc_output_fields_item_ops {
  MRC_SUBCLASS_OPS(struct psc_output_fields_item);
  void (*run)(struct psc_output_fields_item *item,
	      mfields_base_t *flds, mparticles_base_t *particles,
	      mfields_c_t *res);
  int nr_comp;
  char *fld_names[6];
};

#define psc_output_fields_item_ops(item)			\
  ((struct psc_output_fields_item_ops *)((item)->obj.ops))

// ======================================================================

extern struct psc_output_fields_item_ops psc_output_fields_item_j_nc_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_j_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_j_ec_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_e_nc_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_e_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_e_ec_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_h_nc_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_h_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_h_fc_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_jdote_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_poyn_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_e2_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_h2_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_densities_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_v_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_vv_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_photon_n_ops;

#endif
