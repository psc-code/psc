
#ifndef PSC_OUTPUT_FIELDS_ITEM_PRIVATE_H
#define PSC_OUTPUT_FIELDS_ITEM_PRIVATE_H

#include <psc_output_fields_item.h>

struct psc_output_fields_item {
  struct mrc_obj obj;
  struct psc_bnd *bnd;
};

enum {
  POFI_ADD_GHOSTS = 1, // this item needs to have ghost points added to interior points
};

struct psc_output_fields_item_ops {
  MRC_SUBCLASS_OPS(struct psc_output_fields_item);
  int  (*get_nr_components)(struct psc_output_fields_item *item);
  const char *(*get_component_name)(struct psc_output_fields_item *item, int m);
  void (*run)(struct psc_output_fields_item *item,
	      struct psc_fields *flds, struct psc_particles *prts,
	      struct psc_fields *res);
  // when the following are known at compile time, they can just be part
  // of the subclass definition.
  // otherwise, one should overload ::get_nr_components() and 
  // ::get_component_name()
  int nr_comp;
  char *fld_names[6];
  unsigned int flags;
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
extern struct psc_output_fields_item_ops psc_output_fields_item_n_1st_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_v_1st_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_vv_1st_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_n_1st_c_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_v_1st_c_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_vv_1st_c_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_n_1st_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_v_1st_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_vv_1st_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_n_1st_nc_c_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_n_1st_nc_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_n_2nd_nc_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_v_2nd_nc_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_vv_2nd_nc_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_n_photon_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_coll_stats_ops;

#endif
