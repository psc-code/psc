
#include "ggcm_mhd.h"
#include "ggcm_mhd_bnd.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_step.h"

#include <mrc_ts_monitor.h>

extern struct mrc_crds_ops mrc_crds_two_gaussian_ops;
extern struct mrc_crds_ops mrc_crds_gaussian_ops;
extern struct mrc_crds_ops mrc_crds_gaussian_2D_ops;

extern struct mrc_ts_ops mrc_ts_ggcm_ops;

extern struct mrc_ts_monitor_ops mrc_ts_monitor_ggcm_ops;

extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_conducting_ops;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_conducting_x_ops;
extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_open_x_ops;


extern struct ggcm_mhd_step_ops ggcm_mhd_step_cweno_ops;


void
ggcm_mhd_register()
{
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_two_gaussian_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_gaussian_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_gaussian_2D_ops);

  mrc_class_register_subclass(&mrc_class_mrc_ts_monitor, &mrc_ts_monitor_ggcm_ops);

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_conducting_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_conducting_x_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_open_x_ops);

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_step, &ggcm_mhd_step_cweno_ops);

#if 0
  mrc_class_register_subclass(&mrc_class_mrc_ts, &mrc_ts_ggcm_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_bnd, &ggcm_mhd_bnd_harris_ops);
#endif
}
