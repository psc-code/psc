
#include "ggcm_mhd.h"
#include "ggcm_mhd_bnd.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_step.h"

#include <mrc_ts_monitor.h>

extern struct mrc_ts_ops mrc_ts_ggcm_ops;

extern struct mrc_ts_monitor_ops mrc_ts_monitor_ggcm_ops;
extern struct mrc_ts_monitor_ops mrc_ts_monitor_conservation_ops;

extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_conducting_ops;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_conducting_x_ops;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_conducting_y_ops;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_open_x_ops;


extern struct ggcm_mhd_step_ops ggcm_mhd_step_cweno_ops;


void
ggcm_mhd_register()
{
  mrc_class_register_subclass(&mrc_class_mrc_ts_monitor, &mrc_ts_monitor_ggcm_ops);
  mrc_class_register_subclass(&mrc_class_mrc_ts_monitor, &mrc_ts_monitor_conservation_ops);

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_conducting_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_conducting_x_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_conducting_y_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_open_x_ops);

  // force regular subclasses to be registered first
  struct mrc_fld *fld = mrc_fld_create(MPI_COMM_NULL);
  mrc_fld_destroy(fld);
#if 0
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_harris_ops);
#endif
}
