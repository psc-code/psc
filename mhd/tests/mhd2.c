
//#define BOUNDS_CHECK

#include "mhd2.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_ic_private.h"
#include "ggcm_mhd_step_private.h"
#include "ggcm_mhd_flds.h"
#include "ggcm_mhd_crds_private.h"
#include "ggcm_mhd_crds_gen.h"
#include "ggcm_mhd_bnd.h"
#include "ggcm_mhd_diag.h"

#include <mrc_ts.h>
#include <mrc_ts_monitor.h>
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_params.h>
#include <mrc_ddc.h>
#include <mrctest.h>
#include <mrc_io.h> 
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

static void
ggcm_mhd_cweno_create(struct ggcm_mhd *mhd)
{
  mhd->par.rrnorm = 1.f;
  mhd->par.ppnorm = 1.f;
  mhd->par.vvnorm = 1.f;
  mhd->par.bbnorm = 1.f;
  mhd->par.ccnorm = 1.f;
  mhd->par.eenorm = 1.f;
  mhd->par.resnorm = 1.f;
  mhd->par.diffco = 0.f;

  ggcm_mhd_bnd_set_type(mhd->bnd, "none");

  mrc_domain_set_param_int(mhd->domain, "bcx", BC_PERIODIC);
  mrc_domain_set_param_int(mhd->domain, "bcy", BC_PERIODIC);
  mrc_domain_set_param_int(mhd->domain, "bcz", BC_PERIODIC);

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "gaussian_2D");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_float3(crds, "l", (float[3]) {  0.0, 0.0, -1.0 });
  mrc_crds_set_param_float3(crds, "h", (float[3]) {  2.*M_PI, 2.*M_PI,  1.0 });

  /* set defaults for the ddc, this does the communication */
  struct mrc_ddc *ddc = mrc_domain_get_ddc(mhd->domain);
  mrc_ddc_set_param_int(ddc, "max_n_fields", 8);
  mrc_ddc_set_param_int3(ddc, "ibn", (int[3]) { SW_2, SW_2, SW_2 });

  // generate MHD solver grid from mrc_crds
  ggcm_mhd_crds_gen_set_type(mhd->crds->crds_gen, "mrc");
  ggcm_mhd_set_param_float(mhd, "isphere", 0.);
  ggcm_mhd_set_param_float(mhd, "diffsphere", 0.);
  ggcm_mhd_set_param_float(mhd, "speedlimit", 1e9);
}

static struct ggcm_mhd_ops ggcm_mhd_cweno_ops = {
  .name             = "cweno",
  .create           = ggcm_mhd_cweno_create,
};


// ======================================================================

#if 0

/* central difference, ith component of gradient */
/* i is determined by dind, so [1,0,0] is d/dx   */
#define DDXI(fld, m, crds, i, ind, dind)                               \
 ((MRC_F3(fld, (m), ind[0]+dind[0], ind[1]+dind[1], ind[2]+dind[2]) -  \
   MRC_F3(fld, (m), ind[0]-dind[0], ind[1]-dind[1], ind[2]-dind[2])) / \
  (MRC_CRD(crds, i, ind[i]+1) - MRC_CRD(crds, i, ind[i]-1)))

/* central difference 2nd derivative */
#define D2DX2I(fld, m, crds, i, ind, dind)                             \
 ((MRC_F3(fld, (m), ind[0]+dind[0], ind[1]+dind[1], ind[2]+dind[2]) -  \
 2.0*MRC_F3(fld, (m), ind[0]        , ind[1]        , ind[2]        ) +  \
   MRC_F3(fld, (m), ind[0]-dind[0], ind[1]-dind[1], ind[2]-dind[2])) / \
  ((MRC_CRD(crds, i, ind[i]+1) - MRC_CRD(crds, i, ind[i]-1)) *           \
   (MRC_CRD(crds, i, ind[i]+1) - MRC_CRD(crds, i, ind[i]-1)) ))

#define ZIP_AVG2(flda, a, fldb, b, ind, dind)                               \
  (0.5*((MRC_F3(flda, a, ind[0],         ind[1],         ind[2]        ) *  \
         MRC_F3(fldb, b, ind[0]+dind[0], ind[1]+dind[1], ind[2]+dind[2])) + \
        (MRC_F3(flda, a, ind[0]+dind[0], ind[1]+dind[1], ind[2]+dind[2]) *  \
         MRC_F3(fldb, b, ind[0],         ind[1],         ind[2]        ))))

#define ZIP_AVG3(flda, a, fldb, b, fldc, c, ind, dind)                       \
  (0.25*((MRC_F3(flda, a, ind[0],         ind[1],         ind[2]        ) *  \
          MRC_F3(fldb, b, ind[0],         ind[1],         ind[2]        ) *  \
          MRC_F3(fldc, c, ind[0]+dind[0], ind[1]+dind[1], ind[2]+dind[2])) + \
         (MRC_F3(flda, a, ind[0]+dind[0], ind[1]+dind[1], ind[2]+dind[2]) *  \
          MRC_F3(fldb, b, ind[0]+dind[0], ind[1]+dind[1], ind[2]+dind[2]) *  \
          MRC_F3(fldc, c, ind[0],         ind[1],         ind[2]        )) + \
         (MRC_F3(flda, a, ind[0]+dind[0], ind[1]+dind[1], ind[2]+dind[2]) *  \
          MRC_F3(fldb, b, ind[0],         ind[1],         ind[2]        ) *  \
          MRC_F3(fldc, c, ind[0]+dind[0], ind[1]+dind[1], ind[2]+dind[2])) + \
         (MRC_F3(flda, a, ind[0],         ind[1],         ind[2]        ) *  \
          MRC_F3(fldb, b, ind[0]+dind[0], ind[1]+dind[1], ind[2]+dind[2]) *  \
          MRC_F3(fldc, c, ind[0],         ind[1],         ind[2]        ))))

static void
calc_zip_flux_rho(struct mrc_f3 *flux, struct ggcm_mhd *mhd, struct mrc_f3 *fld,
		  struct mrc_f3 *v)
{
  /* d/dt (rho) = -div(rho*v) */
  /* face centered fluxes calculated using zip average */
  mrc_f3_foreach(fld, ix, iy, iz, 0, 1) {
    int ind[3] = { ix, iy, iz };

    for(int i=0; i<3; i++) {
      int dind[3] = {0, 0, 0};
      dind[i] = -1;

      MRC_F3(flux, i, ix,iy,iz) = ZIP_AVG2(fld, _RR1, v, i, ind, dind);
    }
  } mrc_f3_foreach_end;
}

static void
calc_zip_flux_pxyz(struct mrc_f3 *flux, int j, struct ggcm_mhd *mhd, struct mrc_f3 *fld, struct mrc_f3 *v)
{
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  float nu = 1e-2f;
  float k = 1.f;
  float T = 1.f;
  
  /* d/dt (rho*v_j) = -div_i(pi*vj) + ... */
  mrc_f3_foreach(fld, ix, iy, iz, 0, 1) {
    int ind[3] = { ix, iy, iz };
    
    for(int i=0; i<3; i++){
      int dind[3] = {0, 0, 0};
      dind[i] = -1;
      // pi * vj
      MRC_F3(flux, i, ix,iy,iz) = ZIP_AVG2(fld, _RV1X+i, v, j, ind, dind);
      // - Trho delta_ij
      if (i == j) {
	MRC_F3(flux, i, ix,iy,iz) +=
	  k * T *.5 * (MRC_F3(fld, _RR1, ix,iy,iz) +
		       MRC_F3(fld, _RR1, ix+dind[0],iy+dind[1],iz+dind[2]));
      }
      // nu \partial_i pj
      MRC_F3(flux, i, ix,iy,iz) +=
	- nu * (MRC_F3(fld, _RV1X+j, ix,iy,iz) -
		MRC_F3(fld, _RV1X+j, ix+dind[0],iy+dind[1],iz+dind[2])) / 
	(MRC_CRD(crds, i, ind[i]) - MRC_CRD(crds, i, ind[i]+dind[i]));
    }
  } mrc_f3_foreach_end;
}

static void __unused
calc_zip_fluxes(struct mrc_f3 **flux, struct ggcm_mhd *mhd, struct mrc_f3 *fld)
{
  struct mrc_ddc *ddc = mrc_domain_get_ddc(mhd->domain);
  struct mrc_f3 *v = ggcm_mhd_get_fields(mhd, "v", 3);

  mrc_ddc_fill_ghosts(ddc, 0, _RV1Z + 1, fld);

  /* get v_i = p_i / rho  (p here is momentum) */
  mrc_f3_foreach(fld, ix, iy, iz, 1, 1) {
    for(int i=0; i<3; i++) {
      MRC_F3(v, i, ix, iy, iz) =
	MRC_F3(fld, _RV1X+i, ix, iy, iz) / MRC_F3(fld, _RR1, ix, iy, iz);
    }
  } mrc_f3_foreach_end;

  calc_zip_flux_rho(flux[0], mhd, fld, v);
  for(int j=0; j<3; j++){
    calc_zip_flux_pxyz(flux[j+1], j, mhd, fld, v);
  }

  mrc_f3_destroy(v);
}

static void
calc_fluxes(struct mrc_f3 **flux, struct ggcm_mhd *mhd, struct mrc_f3 *fld)
{
  for (int i = 0; i < 3; i++) {
    calc_fluxes_per_face(flux, mhd, fld, i);
  }
}

static void __unused
calc_avg_fluxes(struct mrc_f3 **flux, struct ggcm_mhd *mhd, struct mrc_f3 *fld)
{
  struct mrc_ddc *ddc = mrc_domain_get_ddc(mhd->domain);

  struct mrc_f3 *flux_cc[_RV1Z + 1];
  for (int m = 0; m < _RV1Z + 1; m++) {
    flux_cc[m] = ggcm_mhd_get_fields(mhd, "flux_cc", 3);
  }
  
  calc_fluxes(flux_cc, mhd, fld);

  for (int m = 0; m < _RV1Z + 1; m++) {
    mrc_ddc_fill_ghosts(ddc, 0, 3, flux_cc[m]);

    mrc_f3_foreach(fld, ix,iy,iz, 0, 1) {
      MRC_F3(flux[m], 0, ix,iy,iz) = .5f * (MRC_F3(flux_cc[m], 0,ix-1,iy,iz) +
					    MRC_F3(flux_cc[m], 0,ix  ,iy,iz));
      MRC_F3(flux[m], 1, ix,iy,iz) = .5f * (MRC_F3(flux_cc[m], 1,ix,iy-1,iz) +
					    MRC_F3(flux_cc[m], 1,ix,iy  ,iz));
      MRC_F3(flux[m], 2, ix,iy,iz) = .5f * (MRC_F3(flux_cc[m], 2,ix,iy,iz-1) +
					    MRC_F3(flux_cc[m], 2,ix,iy,iz  ));
    } mrc_f3_foreach_end;
  }

  for (int m = 0; m < _RV1Z + 1; m++) {
    mrc_f3_destroy(flux_cc[m]);
  }
}

static void __unused 
zero_flux(struct mrc_f3 **flux, struct ggcm_mhd *mhd, struct mrc_f3 *u)
{  
  /*
  mrc_f3_foreach(u, ix,iy,iz, 2, 2) {    
    for (int ff = 0; ff < _UU1+1; ff++) {      
      MRC_F3(flux[ff],0,ix,iy,iz) = 0.0 ; 
      //MRC_F3(flux[ff],0,ix,iy,iz) = 0.0 ; 
    } mrc_f3_foreach_end;  
  */
  struct mrc_f3 *f3 = ggcm_mhd_flds_get_mrc_f3(mhd->flds_base);
  const int *dims = mrc_f3_dims(f3);
  int nx = dims[0], ny = dims[1], nz = dims[2];
  int sw = SW_2;
  int gdims[3];
  int bc[3];
  struct mrc_patch_info info;
  //  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_domain_get_global_dims(mhd->domain, gdims);
  mrc_domain_get_local_patch_info(mhd->domain, 0, &info);
  //mrc_domain_get_param_int(crds,'sw', &sw[0] ); 
  mrc_domain_get_param_int(mhd->domain, "bcx", &bc[0]); // FIXME in libmrc
  mrc_domain_get_param_int(mhd->domain, "bcy", &bc[1]);
  mrc_domain_get_param_int(mhd->domain, "bcz", &bc[2]);


  if (bc[1] != BC_PERIODIC && info.off[1] == 0) { // x lo
    for (int iz = -sw; iz < nz + sw; iz++) {
      for (int iy = -sw; iy < 1; iy++) {	 
        for (int ix = -sw; ix < nx + sw; ix++) {	 
	  for (int ff = 0; ff < _JZ+1; ff++) { 
	    MRC_F3(flux[ff],0, ix,iy,iz) = 0.0 ; 
	    MRC_F3(flux[ff],1, ix,iy,iz) = 0.0 ; 	 	 
	    MRC_F3(flux[ff],2, ix,iy,iz) = 0.0 ; 	 
	  //MRC_F3(flux[ff],2, ix,iy,iz) = 0.0 ; 	 	 	 
	  }
        }
      }
    }
  }
  if (bc[1] != BC_PERIODIC && info.off[1] + info.ldims[1] == gdims[1]) { // x hi
    for (int iz = -sw; iz < nz + sw; iz++) { 
      for (int iy = ny; iy < ny+sw; iy++) {	 
	for (int ix = -sw; ix < nx + sw; ix++) {	 
	  for (int ff = 0; ff < _JZ+1; ff++) { 
	    MRC_F3(flux[ff],0, ix,iy,iz) = 0.0 ; 
	    MRC_F3(flux[ff],1, ix,iy,iz) = 0.0 ; 	 	 
	    MRC_F3(flux[ff],2, ix,iy,iz) = 0.0 ; 	 
	    //MRC_F3(flux[ff],2, ix,iy,iz) = 0.0 ; 	 	 	 
	  }
	}
      }
    }
  }
}
#endif
  
static void
calc_mhd_rhs(void *ctx, struct mrc_obj *_rhs, float time, struct mrc_obj *_fld)
{
  struct ggcm_mhd *mhd = ctx;
  struct mrc_f3 *rhs = (struct mrc_f3 *) _rhs;
  struct mrc_f3 *fld = (struct mrc_f3 *) _fld;
  
  mhd->time = time;
  ggcm_mhd_step_calc_rhs(mhd->step, rhs, fld);
}


// ======================================================================

static void
diag_write(void *_mhd, float time, struct mrc_obj *_x, FILE *file)
{
  struct mrc_f3 *x = (struct mrc_f3 *) _x;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0)
    return;

  fprintf(file, "%g", time);
  for (int m = _RR1; m <= _B1Z; m++) {
    fprintf(file, " %g", MRC_F3(x, m, 0,0,0));
  }
  fprintf(file, "\n");
}

// ======================================================================

extern struct mrc_crds_ops mrc_crds_two_gaussian_ops;
extern struct mrc_crds_ops mrc_crds_gaussian_ops;
extern struct mrc_crds_ops mrc_crds_gaussian_2D_ops;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_conducting_ops;
extern struct ggcm_mhd_step_ops ggcm_mhd_step_cweno_ops;

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_cweno_ops);  

  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_two_gaussian_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_gaussian_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_gaussian_2D_ops);

  mrc_class_register_subclass(&mrc_class_mrc_ts_monitor, &mrc_ts_monitor_ggcm_ops);

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_fadeev_ops);  
#if 0
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_whistler_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_bw_ops); 
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_ot_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_otzi_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_hydroblast_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_mhdblast_ops);    
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_harris_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_kh_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_ici_ops); 
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_wave_sound_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_wave_alfven_ops);
 
  mrc_class_register_subclass(&mrc_class_mrc_ts, &mrc_ts_ggcm_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_harris_ops);
#endif
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_conducting_ops);

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_step, &ggcm_mhd_step_cweno_ops);

  struct ggcm_mhd *mhd = ggcm_mhd_create(MPI_COMM_WORLD);
  ggcm_mhd_set_type(mhd, "cweno");
  ggcm_mhd_step_set_type(mhd->step, "cweno");
  ggcm_mhd_flds_set_type(mhd->flds_base, "c");
  ggcm_mhd_set_from_options(mhd);
  ggcm_mhd_setup(mhd);
  ggcm_mhd_view(mhd);

  // set up initial condition
  ggcm_mhd_ic_run(mhd->ic);

  // run time integration
  struct mrc_ts *ts = mrc_ts_create(mrc_domain_comm(mhd->domain));
  mrc_ts_set_type(ts, "rk2");
  mrc_ts_set_context(ts, ggcm_mhd_to_mrc_obj(mhd));

  struct mrc_ts_monitor *mon_output =
    mrc_ts_monitor_create(mrc_ts_comm(ts));
  mrc_ts_monitor_set_type(mon_output, "ggcm");
  mrc_ts_monitor_set_name(mon_output, "mrc_ts_output");
  mrc_ts_add_monitor(ts, mon_output);

  struct mrc_ts_monitor *mon_diag =
    mrc_ts_monitor_create(mrc_ts_comm(ts));
  mrc_ts_monitor_set_type(mon_diag, "diag");
  mrc_ts_monitor_set_name(mon_diag, "mrc_ts_diag");
  mrc_ts_monitor_diag_set_function(mon_diag, diag_write, mhd);
  mrc_ts_add_monitor(ts, mon_diag);

  mrc_ts_set_dt(ts, 1e-6);
  mrc_ts_set_solution(ts, mrc_f3_to_mrc_obj(ggcm_mhd_flds_get_mrc_f3(mhd->flds_base)));
  mrc_ts_set_rhs_function(ts, calc_mhd_rhs, mhd);
  mrc_ts_set_from_options(ts);
  mrc_ts_view(ts);
  mrc_ts_setup(ts);
  mrc_ts_solve(ts);
  mrc_ts_view(ts);
  mrc_ts_destroy(ts);  
  ggcm_mhd_destroy(mhd);

  MPI_Finalize();
  return 0;
}

