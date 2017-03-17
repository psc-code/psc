
#include <ggcm_mhd_private.h>
#include <ggcm_mhd_bnd.h>
#include <ggcm_mhd_diag.h>
#include <ggcm_mhd_ic.h>
#include <ggcm_mhd_step.h>

#include <mrc_domain.h>
#include <mrc_physics.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

// ======================================================================
// ggcm_mhd subclass "mirdip"

// ----------------------------------------------------------------------
// ggcm_mhd_mirdip_create

static void
ggcm_mhd_mirdip_create(struct ggcm_mhd *mhd)
{
  // normalization based on RE, equitorial field, ionospheric density
  mhd->par.norm_length = C_RE; // earth radius in m
  mhd->par.norm_B = C_DIPOLESTRENGTH; // equitorial dipole field in T
  mhd->par.norm_density = 1e10 * C_AMU; // approximate ionospheric density in kg/m^3

  // select output / input unit prefixes as in OpenGGCM
  mhd->par.bbnorm0  = 1e-9;  // uT
  mhd->par.vvnorm0  = 1e3;   // km/s
  mhd->par.rrnorm0  = 1e6 * C_AMU; // amu/cm^3
  mhd->par.ppnorm0  = 1e-12; // pPa
  mhd->par.ccnorm0  = 1e-6;  // uA/m^2
  mhd->par.eenorm0  = 1e-3;  // mV/m
  mhd->par.resnorm0 = 1.f;   // ...
  mhd->par.tnorm0   = 1.f;   // s

  ggcm_mhd_ic_set_type(mhd->ic, "mirdip_double");
  ggcm_mhd_bnd_set_type(mhd->bnd, "inoutflow_sc_double");
  ggcm_mhd_step_set_type(mhd->step , "c3_double"); // FIXME, if not double, the conversion mess up the "view" field

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_double3(crds, "l", (double[3]) {  -20., -20., -20. });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {   20.,  20.,  20. });
}

// ----------------------------------------------------------------------
// ggcm_mhd_mirdip_ops

static struct ggcm_mhd_ops ggcm_mhd_mirdip_ops = {
  .name             = "mirdip",
  .create           = ggcm_mhd_mirdip_create,
};

// ======================================================================
// main

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_mirdip_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
 
  return ggcm_mhd_main(&argc, &argv);
}

