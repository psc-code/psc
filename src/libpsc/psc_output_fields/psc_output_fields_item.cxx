
#include "psc_output_fields_item_private.h"
#include "psc_bnd.h"
#include "bnd.hxx"

#include <string.h>

// ----------------------------------------------------------------------
// psc_output_fields_item_set_psc_bnd

void
psc_output_fields_item_set_psc_bnd(struct psc_output_fields_item *item,
				   struct psc_bnd *bnd)
{
  item->bnd = bnd; // FIXME, ref counting?
}

// ======================================================================
// psc_output_fields_item_init

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
extern struct psc_output_fields_item_ops psc_output_fields_item_divb_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_divj_ops;

extern struct psc_output_fields_item_ops psc_output_fields_item_n_1st_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_p_1st_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_T_1st_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_v_1st_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_vv_1st_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_Tvv_1st_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_nvt_1st_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_nvp_1st_single_ops;

extern struct psc_output_fields_item_ops psc_output_fields_item_n_1st_double_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_p_1st_double_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_T_1st_double_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_v_1st_double_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_vv_1st_double_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_Tvv_1st_double_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_nvt_1st_double_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_nvp_1st_double_ops;


extern struct psc_output_fields_item_ops psc_output_fields_item_n_1st_nc_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_n_1st_nc_double_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_n_2nd_nc_double_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_rho_1st_nc_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_rho_1st_nc_double_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_rho_2nd_nc_double_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_v_1st_nc_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_v_1st_nc_double_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_n_photon_ops;

extern struct psc_output_fields_item_ops psc_output_fields_item_coll_stats_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_coll_rei_single_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_coll_stats_double_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_coll_rei_double_ops;

extern struct psc_output_fields_item_ops psc_output_fields_item_dive_c_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_dive_single_ops;

extern struct psc_output_fields_item_ops psc_output_fields_item_dive_cuda_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_rho_1st_nc_cuda_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_n_1st_cuda_ops;

extern struct psc_output_fields_item_ops psc_output_fields_item_vpic_fields_ops;
extern struct psc_output_fields_item_ops psc_output_fields_item_vpic_hydro_ops;

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
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_divb_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_divj_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_p_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_T_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_v_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_vv_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_Tvv_1st_single_ops);
  // mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_nvt_1st_single_ops);
  // mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_nvp_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_p_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_T_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_v_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_vv_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_Tvv_1st_double_ops);
  // mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_nvt_1st_double_ops);
  // mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_nvp_1st_double_ops);

  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_1st_nc_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_1st_nc_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_2nd_nc_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_rho_1st_nc_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_rho_1st_nc_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_rho_2nd_nc_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_v_1st_nc_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_v_1st_nc_double_ops);

  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_coll_stats_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_coll_rei_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_coll_stats_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_coll_rei_double_ops);

  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_dive_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_dive_single_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_dive_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_rho_1st_nc_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_1st_cuda_ops);
#endif
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_vpic_fields_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_vpic_hydro_ops);
#endif
}

// ======================================================================
// psc_output_fields_item class

struct mrc_class_psc_output_fields_item_ : mrc_class_psc_output_fields_item {
  mrc_class_psc_output_fields_item_() {
    name             = "psc_output_fields_item";
    size             = sizeof(struct psc_output_fields_item);
    init             = psc_output_fields_item_init;
  }
} mrc_class_psc_output_fields_item;

