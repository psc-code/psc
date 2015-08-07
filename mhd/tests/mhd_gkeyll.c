
#include <ggcm_mhd_private.h>
#include <ggcm_mhd_step.h>
#include <ggcm_mhd_ic_private.h>
#include <ggcm_mhd_crds_private.h>
#include <ggcm_mhd_crds_gen.h>
#include <ggcm_mhd_bnd.h>
#include <ggcm_mhd_bndsw.h>
#include <ggcm_mhd_diag.h>

#include <mrc_fld_as_double_aos.h>
#include <mrc_domain.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

#include <ggcm_mhd_diag_item_private.h>


extern struct ggcm_mhd_step_ops ggcm_mhd_step_gkeyll_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_gkeyll_ops;
extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_item_ops_gkeyll_e;
extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_item_ops_gkeyll_i;
extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_item_ops_gkeyll_em;

// ======================================================================
// ggcm_mhd subclass "gkeyll"

// ----------------------------------------------------------------------
// ggcm_mhd_gkeyll_create

static void
ggcm_mhd_gkeyll_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_double3(crds, "l", (double[3]) {  0.0, 0.0, 0.0 });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {  1.0, 1.0, 0.1 });
}

// ----------------------------------------------------------------------
// ggcm_mhd_gkeyll_ops

struct ggcm_mhd_gkeyll {
};

#define VAR(x) (void *)offsetof(struct ggcm_mhd_gkeyll, x)
static struct param ggcm_mhd_gkeyll_descr[] = {
  {},
};
#undef VAR

static struct ggcm_mhd_ops ggcm_mhd_gkeyll_ops = {
  .name             = "gkeyll",
  .size             = sizeof(struct ggcm_mhd_gkeyll),
  .param_descr      = ggcm_mhd_gkeyll_descr,
  .create           = ggcm_mhd_gkeyll_create,
};

// ======================================================================
// main

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_step, &ggcm_mhd_step_gkeyll_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_gkeyll_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_gkeyll_e);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_gkeyll_i);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_gkeyll_em);  

  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_gkeyll_ops);  
 
  return ggcm_mhd_main(&argc, &argv);
}

