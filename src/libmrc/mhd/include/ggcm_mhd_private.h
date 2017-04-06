
#ifndef GGCM_MHD_PRIVATE_H
#define GGCM_MHD_PRIVATE_H

#include "ggcm_mhd.h"

#include "mrc_crds.h"

enum {
  MAGDIFFU_NL1,
  MAGDIFFU_RES1,
  MAGDIFFU_CONST,
};

#define GK_NR_FLUIDS_MAX (5)

struct ggcm_mhd_params {
  float gamm;
  float rrmin;

  double xxnorm0;
  double bbnorm0, vvnorm0, rrnorm0, ppnorm0;
  double ccnorm0, eenorm0, resnorm0, tnorm0;
  double qqnorm0;
  double norm_length; // normalizing length (in m)
  double norm_B; // normalizing magnetic field (in T)
  double norm_density; // normalizing density (in 1/m^3)
  // norm_mu0, sort of like the other quantities, is the ratio between
  // code unit and external unit. In particular, for Alfven-normalized
  // units, norm_mu0 = mu0. In the case of SI-SI, or normalized-normalized,
  // it's equal to 1
  double norm_mu0;
  // This is the mu0 we're using in the actual equations we're solving
  // (traditionally, we're using normalized units, so mu0_code = 1)
  double mu0_code;

  float diffco, diffth;
  float diffsphere;
  float speedlimit, thx;
  float r_db_dt, isphere, timelo;
  float diff_timelo;
  float diff_swbnd;
  int diff_obnd;
  float d_i;
  float dtmin;
  int modnewstep;
  int magdiffu;
  bool do_badval_checks; // check for NaN or negative density / pressure

  bool do_limit2;
  bool do_limit3;
  bool limit_aspect_low;
  bool calce_aspect_low;

  // params for multi-fluid moment runs
  // to be obtained from gkeyll instead
  // of from command line options
  double gk_speed_of_light;
  int gk_nr_fluids;
  int gk_nr_moments;
  // charge, mass per species
  // (will define things like d_i)
  struct mrc_param_float_array gk_charge;
  struct mrc_param_float_array gk_mass;
  // pressure ratios are only used in the initial condition to distribute the MHD
  // total pressure onto species
  struct mrc_param_float_array gk_pressure_ratios;

  // derived quantities reusable to handle gkeyll data
  // first index of all species
  // holding five species at max
  int gk_idx[5];
  // relative mass ratios of all species
  float gk_mass_ratios[5];
  // q/m ratios of all species
  float gk_q_m[5];

  bool gk_norm;
  double gk_norm_speed_of_light;
  double gk_norm_mi_over_me;
  double gk_norm_ppi_over_ppe;
  double gk_norm_rr;

  bool monitor_conservation;
};

struct ggcm_mhd {
  struct mrc_obj obj;
  struct ggcm_mhd_params par;
  int amr; //< turn on if > 0, value selects initial domain refinement
  char *amr_grid_file;  // used if mhd->amr == 999  
  struct mrc_ddc *ddc_amr_cc;
  struct mrc_ddc *ddc_amr_flux[3];
  struct mrc_ddc *ddc_amr_E;

  struct mrc_domain *domain;
  struct mrc_fld *fld;
  struct mrc_fld *ymask;
  struct mrc_fld *bnd_mask;
  // background B field
  // b0 = NULL means there is none
  struct mrc_fld *b0;
  struct ggcm_mhd_crds *crds;
  struct ggcm_mhd_step *step;
  struct ggcm_mhd_bnd *bnd;
  struct ggcm_mhd_bnd *bnd1;
  struct ggcm_mhd_diag *diag;
  struct ggcm_mhd_ic *ic;

  // mhd state
  // normalization parameters
  // multiplying the internal normalized quantities by these will produce
  // physical values in SI units, but with a prefix given by the corresponding
  // XXnorm0 parameter
  double xxnorm;
  double bbnorm, vvnorm, rrnorm, ppnorm;
  double ccnorm, eenorm, resnorm, tnorm;
  double qqnorm;

  float time_code; // current time in code (normalized) units
  float dt_code; // current timestep in code (normalized) units
  int istep;
  float timla;
  double dacttime;

  // for easy access, cached from ::domain
  int im[3];  // local domain excl ghost points
  int img[3]; // local domain incl ghost points
};

struct ggcm_mhd_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd);
  void (*set_state)(struct ggcm_mhd *mhd);
  void (*pre_step)(struct ggcm_mhd *mhd, struct mrc_ts *ts, struct mrc_fld *fld);
  void (*post_step)(struct ggcm_mhd *mhd, struct mrc_ts *ts, struct mrc_fld *fld);
};

extern struct ggcm_mhd_ops ggcm_mhd_ops_box;

// ----------------------------------------------------------------------

void ggcm_mhd_calc_currcc_bgrid_fc_ggcm(struct ggcm_mhd *mhd, struct mrc_fld *f, int m,
					struct mrc_fld *c);
void ggcm_mhd_calc_currcc_bgrid_fc(struct ggcm_mhd *mhd, struct mrc_fld *f, int m,
				   struct mrc_fld *c);
void ggcm_mhd_calc_currcc_bgrid_cc(struct ggcm_mhd *mhd, struct mrc_fld *f, int m,
				   struct mrc_fld *c);
void ggcm_mhd_calc_divb_bgrid_fc_ggcm(struct ggcm_mhd *mhd, struct mrc_fld *f, struct mrc_fld *d);
void ggcm_mhd_calc_divb_bgrid_fc(struct ggcm_mhd *mhd, struct mrc_fld *f, struct mrc_fld *d);
void ggcm_mhd_calc_divb_bgrid_cc(struct ggcm_mhd *mhd, struct mrc_fld *f, struct mrc_fld *d);

void ggcm_mhd_calc_rr_scons(struct ggcm_mhd *mhd, struct mrc_fld *rr, struct mrc_fld *fld);
void ggcm_mhd_calc_rr_fcons_fc(struct ggcm_mhd *mhd, struct mrc_fld *rr, struct mrc_fld *fld);
void ggcm_mhd_calc_rr_fcons_cc(struct ggcm_mhd *mhd, struct mrc_fld *rr, struct mrc_fld *fld);
void ggcm_mhd_calc_rr_gkeyll(struct ggcm_mhd *mhd, struct mrc_fld *rr, struct mrc_fld *fld);

void ggcm_mhd_calc_v_scons(struct ggcm_mhd *mhd, struct mrc_fld *v, struct mrc_fld *fld);
void ggcm_mhd_calc_v_fcons_fc(struct ggcm_mhd *mhd, struct mrc_fld *v, struct mrc_fld *fld);
void ggcm_mhd_calc_v_fcons_cc(struct ggcm_mhd *mhd, struct mrc_fld *v, struct mrc_fld *fld);
void ggcm_mhd_calc_v_gkeyll(struct ggcm_mhd *mhd, struct mrc_fld *v, struct mrc_fld *fld);

void ggcm_mhd_calc_pp_scons(struct ggcm_mhd *mhd, struct mrc_fld *pp, struct mrc_fld *fld);
void ggcm_mhd_calc_pp_fcons_fc(struct ggcm_mhd *mhd, struct mrc_fld *pp, struct mrc_fld *fld);
void ggcm_mhd_calc_pp_fcons_cc(struct ggcm_mhd *mhd, struct mrc_fld *pp, struct mrc_fld *fld);
void ggcm_mhd_calc_pp_gkeyll(struct ggcm_mhd *mhd, struct mrc_fld *pp, struct mrc_fld *fld);

// helpers for subclasses to use

struct mrc_fld *ggcm_mhd_get_3d_fld(struct ggcm_mhd *mhd, int nr_comps);
void ggcm_mhd_put_3d_fld(struct ggcm_mhd *mhd, struct mrc_fld *f);

struct mrc_ddc *ggcm_mhd_create_amr_ddc(struct ggcm_mhd *mhd);
struct mrc_ddc *ggcm_mhd_create_amr_ddc_flux(struct ggcm_mhd *mhd, int d);
struct mrc_ddc *ggcm_mhd_create_amr_ddc_E(struct ggcm_mhd *mhd);
void ggcm_mhd_setup_amr_domain(struct ggcm_mhd *mhd);

// reference implementation only
void ggcm_mhd_amr_fill_ghosts_b(struct ggcm_mhd *mhd, struct mrc_fld *fld);

#endif
