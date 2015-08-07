
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
// ggcm_mhd_ic_ot_run

static void
ggcm_mhd_ic_ot_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct mrc_patch_info pinfo;
  double dx0[3], dx[3];
  mrc_crds_get_dx_base(crds, dx0);

  int gdims[3], p1x, p1y, p1z;
  mrc_domain_get_global_dims(mhd->domain, gdims);
  p1x = (gdims[0] > 1);
  p1y = (gdims[1] > 1);
  p1z = (gdims[2] > 1);

  struct mrc_fld *Az = mrc_domain_fld_create(mhd->domain, 2, "Az");
  mrc_fld_set_type(Az, FLD_TYPE);
  mrc_fld_setup(Az);
  mrc_fld_view(Az);

  mrc_fld_data_t B0 = 1.;
  mrc_fld_data_t rr0 = 25. / 9.;
  mrc_fld_data_t v0 = 1.;
  mrc_fld_data_t pp0 = 5. / 3.;

  /* Initialize vector potential */

  double max_divb = 0.;
  
  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    /* get dx for this patch */
    mrc_domain_get_local_patch_info(mhd->domain, p, &pinfo);
    for (int d = 0; d < 3; d++){
      float refine = 1.0;
      if (pinfo.ldims[d] > 1) {
        refine = 1.0 / (1 << pinfo.level);
      }
      dx[d] = dx0[d] * refine;
    }
    
    /* Initialize vector potential */
    mrc_fld_foreach(fld, ix,iy,iz, 1, 2) {
      mrc_fld_data_t xx = MRC_MCRDX(crds, ix, p) - .5 * dx[0];
      mrc_fld_data_t yy = MRC_MCRDY(crds, iy, p) - .5 * dx[1];
      M3(Az, 0, ix,iy,iz, p) = B0 / (4.*M_PI) * cos(4.*M_PI * xx) + B0 / (2.*M_PI) * cos(2.*M_PI * yy);
    } mrc_fld_foreach_end;
    
    /* Initialize face-centered fields */
    mrc_fld_foreach(fld, ix,iy,iz, 1, 1) {
      BX_(fld, ix,iy,iz, p) =  (M3(Az, 0, ix    , iy+p1y, iz, p) - M3(Az, 0, ix,iy,iz, p)) / dx[1];
      BY_(fld, ix,iy,iz, p) = -(M3(Az, 0, ix+p1x, iy    , iz, p) - M3(Az, 0, ix,iy,iz, p)) / dx[0];
    } mrc_fld_foreach_end;

    /* Initialize density, momentum, total energy */
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      mrc_fld_data_t xx = MRC_MCRDX(crds, ix, p), yy = MRC_MCRDY(crds, iy, p);
      
      RR_(fld, ix,iy,iz, p) = rr0;
      VX_(fld, ix,iy,iz, p) = -v0*sin(2.*M_PI * yy);
      VY_(fld, ix,iy,iz, p) =  v0*sin(2.*M_PI * xx);
      PP_(fld, ix,iy,iz, p) = pp0;
    } mrc_fld_foreach_end;    

    /* calc max divb */    
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      double val =
        (BX(fld, ix+p1x, iy    , iz) - BX(fld, ix,iy,iz)) / dx[0] +
        (BY(fld, ix    , iy+p1y, iz) - BY(fld, ix,iy,iz)) / dx[1];

      max_divb = fmax(max_divb, fabs(val));
    } mrc_fld_foreach_end;
  }

  mprintf("max divb = %g\n", max_divb);

  mrc_fld_destroy(Az);
  mrc_fld_put_as(fld, mhd->fld);

  ggcm_mhd_convert_from_primitive(mhd, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_ot_ops = {
  .name        = "ot",
  .run         = ggcm_mhd_ic_ot_run,
};


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
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_ot_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_gkeyll_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_gkeyll_e);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_gkeyll_i);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag_item, &ggcm_mhd_diag_item_ops_gkeyll_em);  

  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_gkeyll_ops);  
 
  return ggcm_mhd_main(&argc, &argv);
}

