
#ifndef GGCM_MHD_PRIVATE_H
#define GGCM_MHD_PRIVATE_H

#include "ggcm_mhd.h"

enum {
  MAGDIFFU_NL1,
  MAGDIFFU_RES1,
  MAGDIFFU_CONST,
};

struct ggcm_mhd_params {
  float gamm;
  float rrmin;
  float bbnorm, vvnorm, rrnorm, ppnorm;
  float ccnorm, eenorm, resnorm;
  float diffco, diffth;
  float diffsphere;
  float speedlimit, thx;
  float isphere, timelo;
  float d_i;
  float dtmin;
  double dbasetime;
  int modnewstep;
  int magdiffu;
};

struct ggcm_mhd {
  struct mrc_obj obj;
  struct ggcm_mhd_params par;
  struct mrc_domain *domain;
  struct mrc_fld *fld;
  struct ggcm_mhd_crds *crds;
  struct ggcm_mhd_step *step;
  struct ggcm_mhd_commu *commu;
  struct ggcm_mhd_bnd *bnd;
  struct ggcm_mhd_diag *diag;
  struct ggcm_mhd_ic *ic;

  // mhd state
  float time; // current time
  float dt;   // current timestep (parameter to pred/corr, so can be .5 dt)
  int istep;
  float timla;
  double dacttime;

  float bndt; // .5 * current timestep in sec, not alfven times

  // for easy access, cached from ::domain
  int im[3];  // local domain excl ghost points
  int img[3]; // local domain incl ghost points
};

struct ggcm_mhd_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd);
  void (*newstep)(struct ggcm_mhd *, float *dtn);
  void (*get_state)(struct ggcm_mhd *mhd);
  void (*set_state)(struct ggcm_mhd *mhd);
};

extern struct ggcm_mhd_ops ggcm_mhd_ops_box;

// ----------------------------------------------------------------------

void primvar1_c(struct ggcm_mhd *mhd);
void primbb_c(struct ggcm_mhd *mhd);
void zmaskn_c(struct ggcm_mhd *mhd);

#endif
