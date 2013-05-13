
#include "ggcm_mhd.h"
#include "ggcm_mhd_bnd.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_ic.h"

#include <mrc_ts_monitor.h>

extern struct mrc_crds_ops mrc_crds_two_gaussian_ops;
extern struct mrc_crds_ops mrc_crds_gaussian_ops;
extern struct mrc_crds_ops mrc_crds_gaussian_2D_ops;

extern struct mrc_ts_ops mrc_ts_ggcm_ops;

extern struct mrc_ts_monitor_ops mrc_ts_monitor_ggcm_ops;

extern struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_harris_ops;

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


void
ggcm_mhd_register()
{
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_two_gaussian_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_gaussian_ops);
  mrc_class_register_subclass(&mrc_class_mrc_crds, &mrc_crds_gaussian_2D_ops);

  mrc_class_register_subclass(&mrc_class_mrc_ts_monitor, &mrc_ts_monitor_ggcm_ops);

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
}
