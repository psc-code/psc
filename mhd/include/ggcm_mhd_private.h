
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
  float ccnorm, eenorm, resnorm, tnorm;
  float diffco, diffth;
  float diffsphere;
  float speedlimit, thx;
  float isphere, timelo;
  float diff_timelo;
  float diff_swbnd;
  int diff_obnd;
  float d_i;
  float dtmin;
  double dbasetime;
  int modnewstep;
  int magdiffu;

  bool monitor_conservation;
};

struct ggcm_mhd {
  struct mrc_obj obj;
  struct ggcm_mhd_params par;
  int amr; //< turn on if > 0, value selects initial domain refinement
  struct mrc_ddc *ddc_amr_cc;
  struct mrc_ddc *ddc_amr_flux_x;
  struct mrc_ddc *ddc_amr_flux_y;
  struct mrc_ddc *ddc_amr_flux_z;
  struct mrc_ddc *ddc_amr_E;

  struct mrc_domain *domain;
  struct mrc_fld *fld;
  struct mrc_fld *ymask;
  struct ggcm_mhd_crds *crds;
  struct ggcm_mhd_step *step;
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
  void (*get_state)(struct ggcm_mhd *mhd);
  void (*set_state)(struct ggcm_mhd *mhd);
};

extern struct ggcm_mhd_ops ggcm_mhd_ops_box;

// ----------------------------------------------------------------------

// helpers for subclasses to use

struct mrc_fld *ggcm_mhd_get_3d_fld(struct ggcm_mhd *mhd, int nr_comps);
void ggcm_mhd_put_3d_fld(struct ggcm_mhd *mhd, struct mrc_fld *f);


// direct access to coords for a given cell
// (ideally avoided for performance critical parts, because it's slower)
//
// FIXME, this should maybe become a feature of mrc_crds, or go away entirely
// because it duplicates already existing functionality to access coordinates via
// MRC_MCRD macros, though the latter only support cell-centered coords at this
// time

void ggcm_mhd_get_crds_cc(struct ggcm_mhd *mhd, int ix, int iy, int iz, int p,
			  float crd[3]);
void ggcm_mhd_get_crds_nc(struct ggcm_mhd *mhd, int ix, int iy, int iz, int p,
			  float crd[3]);
void ggcm_mhd_get_crds_fc(struct ggcm_mhd *mhd, int ix, int iy, int iz, int p,
			  int d, float crd[3]);
void ggcm_mhd_get_crds_ec(struct ggcm_mhd *mhd, int ix, int iy, int iz, int p,
			  int d, float crd[3]);

void primvar_c(struct ggcm_mhd *mhd, int m_curr);
void primvar_float(struct ggcm_mhd *mhd, int m_curr);
void primvar_double(struct ggcm_mhd *mhd, int m_curr);
void primvar1_c(struct ggcm_mhd *mhd);
void primbb_c(struct ggcm_mhd *mhd, int m_curr);
void primbb_float(struct ggcm_mhd *mhd, int m_curr);
void primbb_double(struct ggcm_mhd *mhd, int m_curr);
void ymaskn_c(struct ggcm_mhd *mhd, struct mrc_fld *ymask);
void zmaskn_c(struct ggcm_mhd *mhd);
void zmaskn_float(struct ggcm_mhd *mhd);
void zmaskn_double(struct ggcm_mhd *mhd);
void newstep(struct ggcm_mhd *mhd, float *dtn);

#endif
