
#include <ggcm_mhd_private.h>
#include <ggcm_mhd_ic_private.h>
#include <ggcm_mhd_bnd.h>
#include <ggcm_mhd_diag.h>
#include <ggcm_mhd_ic.h>
#include <ggcm_mhd_step.h>
#include <ggcm_mhd_defs.h>

#include <mrc_domain.h>
#include <mrc_io.h>
#include <mrc_fld_as_double.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

// ======================================================================
// ggcm_mhd_ic subclass "test"

// ----------------------------------------------------------------------
// ggcm_mhd_ic_test_run

static void
ggcm_mhd_ic_test_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
#if 0
    /* Initialize face-centered fields */
    mrc_fld_foreach(fld, ix,iy,iz, 1, 1) {
      BX_(fld, ix,iy,iz, p) = 0.;
      BY_(fld, ix,iy,iz, p) = 0.;
    } mrc_fld_foreach_end;
#endif

    /* Initialize cell-centered fields */
    mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
      mrc_fld_data_t xx[3] =
       	{ MRC_MCRDX(crds, ix, p), MRC_MCRDY(crds, iy, p), MRC_MCRDZ(crds, iz, p), };
      mrc_fld_data_t rr = sqrt(sqr(xx[0]) + sqr(xx[1]) + sqr(xx[2]));
      
      RR_(fld, ix,iy,iz, p) = rr;
      PP_(fld, ix,iy,iz, p) = rr;
    } mrc_fld_foreach_end;    
  }

  mrc_fld_put_as(fld, mhd->fld);

  ggcm_mhd_convert_from_primitive(mhd, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_test_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_test_ops = {
  .name        = "test",
  .run         = ggcm_mhd_ic_test_run,
};


// ======================================================================
// ggcm_mhd subclass "test"

// ----------------------------------------------------------------------
// ggcm_mhd_test_create

static void
ggcm_mhd_test_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_bnd_set_type(mhd->bnd, "sphere_sc_double");
  ggcm_mhd_step_set_type(mhd->step , "c3_double"); // FIXME, if not double, the conversion mess up the "view" field

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_double3(crds, "l", (double[3]) {  -3., -3., -3. });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {   3.,  3.,  3. });
}

// ----------------------------------------------------------------------
// ggcm_mhd_test_ops

static struct ggcm_mhd_ops ggcm_mhd_test_ops = {
  .name             = "test",
  .create           = ggcm_mhd_test_create,
};

// ======================================================================
// main

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_test_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_test_ops);  
 
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);
  ggcm_mhd_register();

  struct ggcm_mhd *mhd = ggcm_mhd_create(MPI_COMM_WORLD);
  ggcm_mhd_set_from_options(mhd);
  ggcm_mhd_setup(mhd);
  ggcm_mhd_view(mhd);

  ggcm_mhd_ic_run(mhd->ic);
  ggcm_mhd_bnd_fill_ghosts(mhd->bnd, mhd->fld, 0.);
  ggcm_mhd_diag_run_now(mhd->diag, mhd->fld, DIAG_TYPE_3D, 0);

  ggcm_mhd_destroy(mhd);

  MPI_Finalize();
  return 0;
}

