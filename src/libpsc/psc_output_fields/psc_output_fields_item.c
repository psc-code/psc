
#include "psc_output_fields_item_private.h"
#include "psc_bnd.h"

#include <string.h>

// ----------------------------------------------------------------------
// psc_output_fields_item_set_psc_bnd

void
psc_output_fields_item_set_psc_bnd(struct psc_output_fields_item *item,
				   struct psc_bnd *bnd)
{
  item->bnd = bnd; // FIXME, ref counting?
}

// ----------------------------------------------------------------------
// psc_output_fields_item_create_mfields

struct psc_mfields *
psc_output_fields_item_create_mfields(struct psc_output_fields_item *item)
{
  struct psc_output_fields_item_ops *ops = psc_output_fields_item_ops(item);
  struct psc_mfields *flds = psc_mfields_create(psc_output_fields_item_comm(item));
  if (strcmp(psc_output_fields_item_type(item), "n_1st_cuda") == 0) { // FIXME
    psc_mfields_set_type(flds, "cuda");
  } else {
    psc_mfields_set_type(flds, "c");
  }
  psc_mfields_set_param_obj(flds, "domain", ppsc->mrc_domain);
  int nr_comp;
  if (ops->flags & POFI_BY_KIND) {
    nr_comp = ops->nr_comp * ppsc->nr_kinds;
  } else {
    nr_comp = ops->nr_comp;
  }
  psc_mfields_set_param_int(flds, "nr_fields", nr_comp);
  psc_mfields_set_param_int3(flds, "ibn", ppsc->ibn);
  psc_mfields_setup(flds);
  assert(ops->nr_comp <= POFI_MAX_COMPS);
  for (int m = 0; m < nr_comp; m++) {
    if (ops->flags & POFI_BY_KIND) {
      int mm = m % ops->nr_comp;
      int k = m / ops->nr_comp;
      char s[strlen(ops->fld_names[mm]) + strlen(ppsc->kinds[k].name) + 2];
      sprintf(s, "%s_%s", ops->fld_names[mm], ppsc->kinds[k].name);
      psc_mfields_set_comp_name(flds, m, s);
    } else {
      psc_mfields_set_comp_name(flds, m, ops->fld_names[m]);
    }
  }

  return flds;
}

// ----------------------------------------------------------------------
// psc_output_fields_item_run

void
psc_output_fields_item_run(struct psc_output_fields_item *item,
			   struct psc_mfields *flds, struct psc_mparticles *particles,
			   struct psc_mfields *res)
{
  struct psc_output_fields_item_ops *ops = psc_output_fields_item_ops(item);
  assert(ops->run_all);

  ops->run_all(item, flds, particles, res);
  if (ops->flags & POFI_ADD_GHOSTS) {
    assert(item->bnd);
    psc_bnd_add_ghosts(item->bnd, res, 0, res->nr_fields);
  }
}

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
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_divb_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_divj_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_1st_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_p_1st_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_T_1st_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_v_1st_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_vv_1st_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_p_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_T_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_v_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_vv_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_Tvv_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_nvt_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_nvp_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_p_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_T_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_v_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_vv_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_Tvv_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_nvt_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_nvp_1st_double_ops);


  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_1st_nc_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_1st_nc_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_1st_nc_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_2nd_nc_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_rho_1st_nc_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_rho_1st_nc_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_rho_1st_nc_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_rho_2nd_nc_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_2nd_nc_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_v_2nd_nc_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_vv_2nd_nc_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_coll_stats_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_coll_stats_single_ops);

  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_dive_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_dive_single_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_dive_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_rho_1st_nc_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_n_1st_cuda_ops);
#endif
}

// ======================================================================
// psc_output_fields_item class

struct mrc_class_psc_output_fields_item mrc_class_psc_output_fields_item = {
  .name             = "psc_output_fields_item",
  .size             = sizeof(struct psc_output_fields_item),
  .init             = psc_output_fields_item_init,
};

