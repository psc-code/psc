
#include <ggcm_mhd_private.h>
#include <ggcm_mhd_step.h>
#include <ggcm_mhd_ic_private.h>
#include <ggcm_mhd_crds_private.h>
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

// ======================================================================
// ggcm_mhd_ic subclass "bowshock3d"

struct ggcm_mhd_ic_bowshock3d {
  float rho_obstacle; // initial density in obstacle
  float r_obstacle;
  float x_obstacle;
  float rho0; // initial density outside
  float p0; // initial pressure
  float v0; // initial velocity
  float Bx0;
  float By0;
  float Bz0;
};

// ----------------------------------------------------------------------
// ggcm_mhd_ic_bowshock3d_run
//
// assuming domain [-1.5,-1.5]x[-1.5,1.5]x[-0.75,0.75]

static mrc_fld_data_t linear_change(mrc_fld_data_t x,
    mrc_fld_data_t x1, mrc_fld_data_t val1,
    mrc_fld_data_t x2, mrc_fld_data_t val2)
{
  mrc_fld_data_t val = val1;
  if (x > x2) {
    val = val2;
  } else if (x > x1) {
    val = val1 + (val2 - val1) * (x - x1) / (x2 - x1);
  }
  return val;
}

static void
ggcm_mhd_ic_bowshock3d_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_bowshock3d *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_bowshock3d);
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

  struct mrc_fld *A = mrc_domain_fld_create(mhd->domain, SW_2, NULL);
  mrc_fld_set_type(A, FLD_TYPE);
  mrc_fld_set_param_int(A, "nr_comps", 3);
  mrc_fld_setup(A);
  mrc_fld_view(A);

  double l[3];
  mrc_crds_get_param_double3(mrc_domain_get_crds(mhd->domain), "l", l);

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
      M3(A, 0, ix,iy,iz, p) = 
        linear_change(xx, l[0], - sub->Bz0 * yy,
            sub->x_obstacle - sub->r_obstacle, 0.);
      M3(A, 2, ix,iy,iz, p) = 
        linear_change(xx, l[0], + sub->Bx0 * yy - sub->By0 * xx,
            sub->x_obstacle - sub->r_obstacle, 0.);
    } mrc_fld_foreach_end;

    /* Initialize face-centered fields B = curl(A) */
    mrc_fld_foreach(fld, ix,iy,iz, 1, 1) {
      BX_(fld, ix,iy,iz, p) =
        + (M3(A, 2, ix    , iy+p1y, iz, p) - M3(A, 2, ix,iy,iz, p)) / dx[1]
        - (M3(A, 1, ix    , iy, iz+p1z, p) - M3(A, 1, ix,iy,iz, p)) / dx[2];
      BY_(fld, ix,iy,iz, p) =
        + (M3(A, 0, ix, iy    , iz+p1z, p) - M3(A, 0, ix,iy,iz, p)) / dx[2]
        - (M3(A, 2, ix+p1x, iy    , iz, p) - M3(A, 2, ix,iy,iz, p)) / dx[0];
      BZ_(fld, ix,iy,iz, p) =
        + (M3(A, 1, ix+p1x, iy    , iz, p) - M3(A, 1, ix,iy,iz, p)) / dx[0]
        - (M3(A, 0, ix, iy    , iz+p1y, p) - M3(A, 0, ix,iy,iz, p)) / dx[1];
    } mrc_fld_foreach_end;

    /* Initialize density, momentum, total energy */
    mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {
      mrc_fld_data_t xx = MRC_CRD(crds, 0, ix) - sub->x_obstacle;
      mrc_fld_data_t yy = MRC_CRD(crds, 1, iy);
      mrc_fld_data_t zz = MRC_CRD(crds, 2, iz);

      PP_(fld, ix, iy, iz, p) = sub->p0;
      if( sqrt((xx*xx) + (yy*yy) + (zz*zz)) <= sub->r_obstacle ){
        RR_(fld, ix, iy, iz, p) = sub->rho_obstacle;
        VX_(fld, ix, iy, iz, p) = 0.;
      } else{
        RR_(fld, ix, iy, iz, p) = sub->rho0;
        VX_(fld, ix, iy, iz, p) = sub->v0;
        if( xx > 0. && sqrt((yy*yy) + (zz*zz)) <= sub->r_obstacle ){
          VX_(fld, ix, iy, iz, p) = 0.;
        }
      }
    } mrc_fld_foreach_end;
  }

  mrc_fld_destroy(A);
  mrc_fld_put_as(fld, mhd->fld);

  ggcm_mhd_convert_from_primitive(mhd, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_bowshock3d_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_bowshock3d, x)
static struct param ggcm_mhd_ic_bowshock3d_descr[] = {
  { "rho_obstacle", VAR(rho_obstacle), PARAM_FLOAT(100.f)   },
  { "r_obstacle"  , VAR(r_obstacle)  , PARAM_FLOAT(.0625f)  },
  { "x_obstacle"  , VAR(x_obstacle)  , PARAM_FLOAT(-.75)    },
  { "rho0"        , VAR(rho0)        , PARAM_FLOAT(.01f)    },
  { "p0"          , VAR(p0)          , PARAM_FLOAT(.0015f)  },
  { "v0"          , VAR(v0)          , PARAM_FLOAT(1.f)     },
  { "Bx0"         , VAR(Bx0)         , PARAM_FLOAT(0.f)     },
  { "By0"         , VAR(By0)         , PARAM_FLOAT(0.001f)  },
  { "Bz0"         , VAR(Bz0)         , PARAM_FLOAT(0.f)     },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_bowshock3d_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_bowshock3d_ops = {
  .name        = "bowshock3d",
  .size        = sizeof(struct ggcm_mhd_ic_bowshock3d),
  .param_descr = ggcm_mhd_ic_bowshock3d_descr,
  .run         = ggcm_mhd_ic_bowshock3d_run,
};

// ======================================================================
// ggcm_mhd_ic subclass "ot"

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_run
//
// assuming 2d domain [-1,1]x[-1,1]

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
        (BX_(fld, ix+p1x, iy    , iz, p) - BX_(fld, ix,iy,iz, p)) / dx[0] +
        (BY_(fld, ix    , iy+p1y, iz, p) - BY_(fld, ix,iy,iz, p)) / dx[1];

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

// ======================================================================
// ggcm_mhd subclass "gkeyll"

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
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_bowshock3d_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_ot_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_gkeyll_ops);  

  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_gkeyll_ops);  
 
  return ggcm_mhd_main(&argc, &argv);
}

