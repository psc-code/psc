
//#define BOUNDS_CHECK

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_step.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_ic.h"

#include <mrc_ts.h>
#include <mrc_ts_monitor.h>
#include <mrc_fld.h>
#include <mrc_domain.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

static void
ggcm_mhd_cweno_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "gaussian_2D");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_float3(crds, "l", (float[3]) {  0.0, 0.0, -1.0 });
  mrc_crds_set_param_float3(crds, "h", (float[3]) {  2.*M_PI, 2.*M_PI,  1.0 });
}

static struct ggcm_mhd_ops ggcm_mhd_cweno_ops = {
  .name             = "cweno",
  .create           = ggcm_mhd_cweno_create,
};

// ======================================================================

static void
diag_write(void *_mhd, float time, struct mrc_obj *_x, FILE *file)
{
  struct mrc_fld *x = (struct mrc_fld *) _x;

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

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_mirdip_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_mirdip2_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_mirdip3_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_whistler_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_otzi_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_ot_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_harris_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_fadeev_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_bw_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_hydroblast_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_mhdblast_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_ici_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_harris;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_kh_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_wave_sound_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_wave_alfven_ops;

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);
  ggcm_mhd_register();

  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_cweno_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_whistler_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_bw_ops); 
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_ot_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_otzi_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_hydroblast_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_mhdblast_ops);    
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_kh_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_ici_ops); 
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_wave_sound_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_wave_alfven_ops);
 
  struct ggcm_mhd *mhd = ggcm_mhd_create(MPI_COMM_WORLD);
  ggcm_mhd_set_type(mhd, "cweno");
  ggcm_mhd_step_set_type(mhd->step, "cweno");
  mrc_fld_set_type(mhd->fld, "float");
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
  mrc_ts_set_solution(ts, mrc_fld_to_mrc_obj(mhd->fld));
  mrc_ts_set_rhs_function(ts, ts_ggcm_mhd_step_calc_rhs, mhd);
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

