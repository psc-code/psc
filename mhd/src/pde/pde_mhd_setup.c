
#ifndef PDE_MHD_SETUP_C
#define PDE_MHD_SETUP_C

#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_gkeyll.h"

#include "pde/pde_setup.c"

// ======================================================================
// MHD parameters, we keep these around statically

static mrc_fld_data_t s_gamma;  // adiabatic exponent
static mrc_fld_data_t s_gamma_m1;  // adiabatic exponent - 1
static mrc_fld_data_t s_gamma_m1_inv;  // 1 / (adiabatic exponent - 1)
static mrc_fld_data_t s_eta;    // (constant) resistivity 
static mrc_fld_data_t s_d_i;    // ion skin depth
static mrc_fld_data_t s_cfl;    // CFL number
static mrc_fld_data_t s_mu0;    // mu0 (in code units)
static mrc_fld_data_t s_mu0_inv;// 1 /mu0 (in code units)

static int s_gk_nr_fluids;
static int *s_gk_idx;
static int s_gk_idx_em;
// FIXME, this should be changed to mrc_fld_data_t
static float *s_gk_mass_ratios;
static float *s_gk_pressure_ratios;

// s_n_state is mostly the same as s_n_comps, except for legacy ggcm which
// has the one field with so many components
static int s_n_state; // number of components in the state vector

static int s_magdiffu;
static mrc_fld_data_t s_diffco; // same as s_eta, but not normalized
static mrc_fld_data_t s_diff_swbnd;
static int s_diff_obnd;
static mrc_fld_data_t s_diff_timelo;
static mrc_fld_data_t s_diffsphere;
static mrc_fld_data_t s_diffth;
static mrc_fld_data_t s_timelo;
static mrc_fld_data_t s_speedlimit_code;
static mrc_fld_data_t s_isphere;
static bool s_do_badval_checks;

static bool s_do_limit2;
static bool s_do_limit3;
static bool s_limit_aspect_low;
static bool s_calce_aspect_low;

// FIXME, these could/should be s_opt_*
static mrc_fld_data_t s_divb_glm_alpha; // ratio of hyperbolic / parabolic divb timescales
static mrc_fld_data_t s_divb_glm_ch_fac; // multiply ch by this factor
static mrc_fld_data_t s_limiter_mc_beta; // beta for MC limiter

// ----------------------------------------------------------------------

static mrc_fld_data_t s_divb_glm_ch _mrc_unused; // hyperbolic c_h from Dedner et al for divb cleaning

// ======================================================================
// options

struct mhd_options {
  int eqn;
  int limiter;
  int riemann;
  int divb;
  int resistivity;
  int hall;
  int time_integrator;
  int get_dt;
  bool background;
  bool bc_reconstruct;

  int mhd_primvar;
  int mhd_primbb;
  int mhd_zmaskn;
  int mhd_rmaskn;
  int mhd_newstep;
  int mhd_pushpred;
  int mhd_pushcorr;
  int mhd_pushfluid1;
  int mhd_pushfluid2;
  int mhd_pushfield1;
  int mhd_pushfield2;
  int mhd_push_ej;
  int mhd_pfie3;
  int mhd_bpush1;
  int mhd_calce;
  int mhd_calc_resis;

  double divb_glm_alpha;
  double divb_glm_ch_fac;
  double limiter_mc_beta;
};

// ----------------------------------------------------------------------
// eqn

static struct mrc_param_select opt_eqn_descr[] _mrc_unused = {
  { .val = OPT_EQN_MHD_FCONS    , .str = "mhd_fcons"   },
  { .val = OPT_EQN_MHD_SCONS    , .str = "mhd_scons"   },
  { .val = OPT_EQN_HD           , .str = "hd"          },
  {},
};

#ifdef OPT_EQN
static const int s_opt_eqn = OPT_EQN;
#else
static int s_opt_eqn _mrc_unused;
#endif

// ----------------------------------------------------------------------
// limiter

static struct mrc_param_select opt_limiter_descr[] _mrc_unused = {
  { .val = OPT_LIMITER_FLAT     , .str = "flat"        },
  { .val = OPT_LIMITER_MINMOD   , .str = "minmod"      },
  { .val = OPT_LIMITER_MC       , .str = "mc"          },
  { .val = OPT_LIMITER_GMINMOD  , .str = "gminmod"     },
  {},
};

#ifdef OPT_LIMITER
static const int s_opt_limiter = OPT_LIMITER;
#else
static int s_opt_limiter _mrc_unused;
#endif

// ----------------------------------------------------------------------
// resistivity

static struct mrc_param_select opt_resistivity_descr[] _mrc_unused = {
  { .val = OPT_RESISTIVITY_NONE     , .str = "none"        },
  { .val = OPT_RESISTIVITY_CONST    , .str = "const"        },
  {},
};

#ifdef OPT_RESISTIVITY
static const int s_opt_resistivity = OPT_RESISTIVITY;
#else
static int s_opt_resistivity _mrc_unused;
#endif

// ----------------------------------------------------------------------
// hall

static struct mrc_param_select opt_hall_descr[] _mrc_unused = {
  { .val = OPT_HALL_NONE     , .str = "none"        },
  { .val = OPT_HALL_CONST    , .str = "const"       },
  { .val = OPT_HALL_YES      , .str = "yes"         },
  {},
};

#ifdef OPT_HALL
static const int s_opt_hall = OPT_HALL;
#else
static int s_opt_hall _mrc_unused;
#endif

// ----------------------------------------------------------------------
// riemann

static struct mrc_param_select opt_riemann_descr[] _mrc_unused = {
  { .val = OPT_RIEMANN_RUSANOV  , .str = "rusanov"     },
  { .val = OPT_RIEMANN_HLL      , .str = "hll"         },
  { .val = OPT_RIEMANN_HLLC     , .str = "hllc"        },
  { .val = OPT_RIEMANN_HLLD     , .str = "hlld"        },
  {},
};

#ifdef OPT_RIEMANN
static const int s_opt_riemann = OPT_RIEMANN;
#else
static int s_opt_riemann _mrc_unused;
#endif

// ----------------------------------------------------------------------
// divb

static struct mrc_param_select opt_divb_descr[] _mrc_unused = {
  { .val = OPT_DIVB_NONE     , .str = "none"        },
  { .val = OPT_DIVB_GLM      , .str = "glm"         },
  {},
};

#ifdef OPT_DIVB
static const int s_opt_divb = OPT_DIVB;
#else
static int s_opt_divb _mrc_unused;
#endif

// ----------------------------------------------------------------------
// time_integrator

static struct mrc_param_select opt_time_integrator_descr[] _mrc_unused = {
  { .val = OPT_TIME_INTEGRATOR_EULER   , .str = "euler"     },
  { .val = OPT_TIME_INTEGRATOR_PREDCORR, .str = "predcorr"  },
  { .val = OPT_TIME_INTEGRATOR_TVD_RK2 , .str = "tvd_rk2"   },
  {},
};

#ifdef OPT_TIME_INTEGRATOR
static const int s_opt_time_integrator = OPT_TIME_INTEGRATOR;
#else
static int s_opt_time_integrator _mrc_unused;
#endif

// ----------------------------------------------------------------------
// get_dt

static struct mrc_param_select opt_get_dt_descr[] _mrc_unused = {
  { .val = OPT_GET_DT_MHD_GGCM, .str = "mhd_ggcm"     },
  { .val = OPT_GET_DT_MHD     , .str = "mhd"          },
  { .val = OPT_GET_DT_MHD_CT  , .str = "mhd_ct"       },
  { .val = OPT_GET_DT_HD      , .str = "hd"           },
  {},
};

#ifdef OPT_GET_DT
static const int s_opt_get_dt = OPT_GET_DT;
#else
static int s_opt_get_dt _mrc_unused;
#endif

// ----------------------------------------------------------------------
// background

#ifdef OPT_BACKGROUND
static const bool s_opt_background = OPT_BACKGROUND;
#else
static bool s_opt_background _mrc_unused;
#endif

// ----------------------------------------------------------------------
// bc_reconstruct

#ifdef OPT_BC_RECONSTRUCT
static const bool s_opt_bc_reconstruct = OPT_BC_RECONSTRUCT;
#else
static bool s_opt_bc_reconstruct _mrc_unused;
#endif

// ======================================================================
// (legacy) mhd options

static struct mrc_param_select opt_mhd_descr[] _mrc_unused = {
  { .val = OPT_MHD_C      , .str = "c"           },
  { .val = OPT_MHD_FORTRAN, .str = "fortran"     },
  { .val = OPT_MHD_C_V2   , .str = "c_v2"        },
  {},
};

// ----------------------------------------------------------------------
// mhd_primvar

#ifdef OPT_MHD_PRIMVAR
static const int s_opt_mhd_primvar = OPT_MHD_PRIMVAR;
#else
static int s_opt_mhd_primvar _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_primbb

#ifdef OPT_MHD_PRIMBB
static const int s_opt_mhd_primbb = OPT_MHD_PRIMBB;
#else
static int s_opt_mhd_primbb _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_zmaskn

#ifdef OPT_MHD_ZMASKN
static const int s_opt_mhd_zmaskn = OPT_MHD_ZMASKN;
#else
static int s_opt_mhd_zmaskn _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_rmaskn

#ifdef OPT_MHD_RMASKN
static const int s_opt_mhd_rmaskn = OPT_MHD_RMASKN;
#else
static int s_opt_mhd_rmaskn _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_newstep

#ifdef OPT_MHD_NEWSTEP
static const int s_opt_mhd_newstep = OPT_MHD_NEWSTEP;
#else
static int s_opt_mhd_newstep _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_pushpred

#ifdef OPT_MHD_PUSHPRED
static const int s_opt_mhd_pushpred = OPT_MHD_PUSHPRED;
#else
static int s_opt_mhd_pushpred _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_pushcorr

#ifdef OPT_MHD_PUSHCORR
static const int s_opt_mhd_pushcorr = OPT_MHD_PUSHCORR;
#else
static int s_opt_mhd_pushcorr _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_pushfluid1

#ifdef OPT_MHD_PUSHFLUID1
static const int s_opt_mhd_pushfluid1 = OPT_MHD_PUSHFLUID1;
#else
static int s_opt_mhd_pushfluid1 _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_pushfluid2

#ifdef OPT_MHD_PUSHFLUID2
static const int s_opt_mhd_pushfluid2 = OPT_MHD_PUSHFLUID2;
#else
static int s_opt_mhd_pushfluid2 _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_pushfield1

#ifdef OPT_MHD_PUSHFIELD1
static const int s_opt_mhd_pushfield1 = OPT_MHD_PUSHFIELD1;
#else
static int s_opt_mhd_pushfield1 _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_pushfield2

#ifdef OPT_MHD_PUSHFIELD2
static const int s_opt_mhd_pushfield2 = OPT_MHD_PUSHFIELD2;
#else
static int s_opt_mhd_pushfield2 _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_push_ej

#ifdef OPT_MHD_PUSH_EJ
static const int s_opt_mhd_push_ej = OPT_MHD_PUSH_EJ;
#else
static int s_opt_mhd_push_ej _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_pfie3

#ifdef OPT_MHD_PFIE3
static const int s_opt_mhd_pfie3 = OPT_MHD_PFIE3;
#else
static int s_opt_mhd_pfie3 _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_bpush1

#ifdef OPT_MHD_BPUSH1
static const int s_opt_mhd_bpush1 = OPT_MHD_BPUSH1;
#else
static int s_opt_mhd_bpush1 _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_calce

#ifdef OPT_MHD_CALCE
static const int s_opt_mhd_calce = OPT_MHD_CALCE;
#else
static int s_opt_mhd_calce _mrc_unused;
#endif

// ----------------------------------------------------------------------
// mhd_calc_resis

#ifdef OPT_MHD_CALC_RESIS
static const int s_opt_mhd_calc_resis = OPT_MHD_CALC_RESIS;
#else
static int s_opt_mhd_calc_resis _mrc_unused;
#endif

// ======================================================================
// calculated options follow

// ----------------------------------------------------------------------
// need_current (calculated)

// FIXME? these aren't really all cases, e.g., OPT_RESISTIVITY != NONE by itself would be enough for yes
#if defined(OPT_RESISTIVITY) && defined(OPT_HALL)
#if OPT_RESISTIVITY != OPT_RESISTIVITY_NONE || OPT_HALL != OPT_HALL_NONE
#define OPT_NEED_CURRENT 1
#else
#define OPT_NEED_CURRENT 0
#endif
#endif

#ifdef OPT_NEED_CURRENT
static const bool s_opt_need_current = OPT_NEED_CURRENT;
#else
static bool s_opt_need_current _mrc_unused;
#endif

// ----------------------------------------------------------------------
// pde_mhd_set_options

static void _mrc_unused
pde_mhd_set_options(struct ggcm_mhd *mhd, struct mhd_options *opt)
{
  // limiter
#ifdef OPT_LIMITER
  assert(OPT_LIMITER == opt->limiter);
#else
  s_opt_limiter = opt->limiter;
#endif

  // resistivity
#ifdef OPT_RESISTIVITY
  assert(OPT_RESISTIVITY == opt->resistivity);
#else
  s_opt_resistivity = opt->resistivity;
#endif

  // hall
#ifdef OPT_HALL
  assert(OPT_HALL == opt->hall);
#else
  s_opt_hall = opt->hall;
#endif

  // riemann
#ifdef OPT_RIEMANN
  assert(OPT_RIEMANN == opt->riemann);
#else
  s_opt_riemann = opt->riemann;
#endif

  // divb
#ifdef OPT_DIVB
  assert(OPT_DIVB == opt->divb);
#else
  s_opt_divb = opt->divb;
#endif

  // time_integrator
#ifdef OPT_TIME_INTEGRATOR
  assert(OPT_TIME_INTEGRATOR == opt->time_integrator);
#else
  s_opt_time_integrator = opt->time_integrator;
#endif

  // get_dt
#ifdef OPT_GET_DT
  assert(OPT_GET_DT == opt->get_dt);
#else
  s_opt_get_dt = opt->get_dt;
#endif

  // background
#ifdef OPT_BACKGROUND
  assert(OPT_BACKGROUND == opt->background);
#else
  s_opt_background = opt->background;
#endif

  // bc_reconstruct
#ifdef OPT_BC_RECONSTRUCT
  assert(OPT_BC_RECONSTRUCT == opt->bc_reconstruct);
#else
  s_opt_bc_reconstruct = opt->bc_reconstruct;
#endif

  // mhd_primvar
#ifdef OPT_MHD_PRIMVAR
  assert(OPT_MHD_PRIMVAR == opt->mhd_primvar);
#else
  s_opt_mhd_primvar = opt->mhd_primvar;
#endif

  // mhd_primbb
#ifdef OPT_MHD_PRIMBB
  assert(OPT_MHD_PRIMBB == opt->mhd_primbb);
#else
  s_opt_mhd_primbb = opt->mhd_primbb;
#endif

  // mhd_zmaskn
#ifdef OPT_MHD_ZMASKN
  assert(OPT_MHD_ZMASKN == opt->mhd_zmaskn);
#else
  s_opt_mhd_zmaskn = opt->mhd_zmaskn;
#endif

  // mhd_rmaskn
#ifdef OPT_MHD_RMASKN
  assert(OPT_MHD_RMASKN == opt->mhd_rmaskn);
#else
  s_opt_mhd_rmaskn = opt->mhd_rmaskn;
#endif

  // mhd_newstep
#ifdef OPT_MHD_NEWSTEP
  assert(OPT_MHD_NEWSTEP == opt->mhd_newstep);
#else
  s_opt_mhd_newstep = opt->mhd_newstep;
#endif

  // mhd_pushpred
#ifdef OPT_MHD_PUSHPRED
  assert(OPT_MHD_PUSHPRED == opt->mhd_pushpred);
#else
  s_opt_mhd_pushpred = opt->mhd_pushpred;
#endif

  // mhd_pushcorr
#ifdef OPT_MHD_PUSHCORR
  assert(OPT_MHD_PUSHCORR == opt->mhd_pushcorr);
#else
  s_opt_mhd_pushcorr = opt->mhd_pushcorr;
#endif

  // mhd_pushfluid1
#ifdef OPT_MHD_PUSHFLUID1
  assert(OPT_MHD_PUSHFLUID1 == opt->mhd_pushfluid1);
#else
  s_opt_mhd_pushfluid1 = opt->mhd_pushfluid1;
#endif

  // mhd_pushfluid1
#ifdef OPT_MHD_PUSHFLUID2
  assert(OPT_MHD_PUSHFLUID2 == opt->mhd_pushfluid2);
#else
  s_opt_mhd_pushfluid2 = opt->mhd_pushfluid2;
#endif

  // mhd_pushfield1
#ifdef OPT_MHD_PUSHFIELD1
  assert(OPT_MHD_PUSHFIELD1 == opt->mhd_pushfield1);
#else
  s_opt_mhd_pushfield1 = opt->mhd_pushfield1;
#endif

  // mhd_pushfield2
#ifdef OPT_MHD_PUSHFIELD2
  assert(OPT_MHD_PUSHFIELD2 == opt->mhd_pushfield2);
#else
  s_opt_mhd_pushfield2 = opt->mhd_pushfield2;
#endif

  // mhd_push_ej
#ifdef OPT_MHD_PUSH_EJ
  assert(OPT_MHD_PUSH_EJ == opt->mhd_push_ej);
#else
  s_opt_mhd_push_ej = opt->mhd_push_ej;
#endif

  // mhd_pfie3
#ifdef OPT_MHD_PFIE3
  assert(OPT_MHD_PFIE3 == opt->mhd_pfie3);
#else
  s_opt_mhd_pfie3 = opt->mhd_pfie3;
#endif

  // mhd_bpush1
#ifdef OPT_MHD_BPUSH1
  assert(OPT_MHD_BPUSH1 == opt->mhd_bpush1);
#else
  s_opt_mhd_bpush1 = opt->mhd_bpush1;
#endif

  // mhd_calce
#ifdef OPT_MHD_CALCE
  assert(OPT_MHD_CALCE == opt->mhd_calce);
#else
  s_opt_mhd_calce = opt->mhd_calce;
#endif

  // mhd_calc_resis
#ifdef OPT_MHD_CALC_RESIS
  assert(OPT_MHD_CALC_RESIS == opt->mhd_calc_resis);
#else
  s_opt_mhd_calc_resis = opt->mhd_calc_resis;
#endif

  // ----------------------------------------------------------------------
  // now set "calculated options"
  // this is replicated runtime code from compile time above

  // need current
#ifndef OPT_NEED_CURRENT
  s_opt_need_current = (s_opt_resistivity != OPT_RESISTIVITY_NONE ||
			s_opt_hall != OPT_HALL_NONE);
#endif

  // ----------------------------------------------------------------------
  // parameters

  s_divb_glm_alpha = opt->divb_glm_alpha;
  s_divb_glm_ch_fac = opt->divb_glm_ch_fac;
  s_limiter_mc_beta = opt->limiter_mc_beta;
}

// ----------------------------------------------------------------------
// pde_mhd_setup

static void
pde_mhd_setup(struct ggcm_mhd *mhd, int n_comps)
{
  pde_setup(mhd->fld, n_comps);
  
  // general (x)mhd params
  s_gamma = mhd->par.gamm;
  s_gamma_m1 = s_gamma - 1.f;
  s_gamma_m1_inv = 1.f / s_gamma_m1;
  // FIXME, this isn't really a good place to normalize
  s_eta   = mhd->par.diffco / mhd->resnorm;
  s_d_i   = mhd->par.d_i;
  s_cfl   = mhd->par.thx;
  s_mu0   = mhd->par.mu0_code;
  s_mu0_inv = 1.f / s_mu0;

  s_n_state = mrc_fld_nr_comps(mhd->fld);
  // FIXME, hacky as usual, to deal with the legacy all-in-one big array
  if (s_n_state == _NR_FLDS) {
    s_n_state = 8;
  }

  // gkeyll parameters
  assert(ggcm_mhd_gkeyll_nr_moments(mhd) == 5);
  s_gk_nr_fluids = mhd->par.gk_nr_fluids;
  s_gk_idx = mhd->par.gk_idx;
  s_gk_idx_em = ggcm_mhd_gkeyll_em_fields_index(mhd);
  s_gk_mass_ratios = mhd->par.gk_mass_ratios;
  s_gk_pressure_ratios = mhd->par.gk_pressure_ratios.vals;

  // openggcm specific params
  s_magdiffu = mhd->par.magdiffu;
  s_diffco = mhd->par.diffco;
  s_diff_swbnd = mhd->par.diff_swbnd;
  s_diff_obnd = mhd->par.diff_obnd;
  s_diff_timelo = mhd->par.diff_timelo;
  s_diffsphere = mhd->par.diffsphere;
  s_diffth = mhd->par.diffth;
  s_timelo = mhd->par.timelo;
  s_speedlimit_code = mhd->par.speedlimit / mhd->vvnorm;
  s_isphere = mhd->par.isphere;
  s_do_badval_checks = mhd->par.do_badval_checks;

  s_do_limit2 = mhd->par.do_limit2;
  s_do_limit3 = mhd->par.do_limit3;
  s_limit_aspect_low = mhd->par.limit_aspect_low;
  s_calce_aspect_low = mhd->par.calce_aspect_low;
}

// ======================================================================
// FIXME this is not a good place

#define BT_(p_B, d, i,j,k)  (F3S(p_B, d, i,j,k) + (s_opt_background ? F3S(s_p_aux.b0, d, i,j,k) : 0))

#define BTX_(U, i,j,k, p)  (s_opt_background ? (BX_(U, i,j,k, p) + B0(b0, 0, i,j,k, p)) : BX_(U, i,j,k, p))
#define BTY_(U, i,j,k, p)  (s_opt_background ? (BY_(U, i,j,k, p) + B0(b0, 1, i,j,k, p)) : BY_(U, i,j,k, p))
#define BTZ_(U, i,j,k, p)  (s_opt_background ? (BZ_(U, i,j,k, p) + B0(b0, 2, i,j,k, p)) : BZ_(U, i,j,k, p))

// ======================================================================
// mhd auxiliary fields
//
// kept around statically so we don't have to pass all this crap

struct mhd_aux {
  // background B field
  fld1d_vec_t b0;
  // bnd_mask for specifying ghost point (1), next-to-ghost point (2)
  fld1d_t bnd_mask;
  // current density for Hall term, resistivity
  fld1d_vec_t j;
};

static struct mhd_aux s_aux;

static void _mrc_unused
pde_mhd_aux_setup()
{
  fld1d_setup(&s_aux.bnd_mask);
  fld1d_vec_setup(&s_aux.j);
  fld1d_vec_setup(&s_aux.b0);
}

// ----------------------------------------------------------------------
// s_p_aux

struct mhd_p_aux {
  fld3d_t b0;
  fld3d_t bnd_mask;
  fld3d_t Jcc; // needed for Hall and constant resistivity
};

static struct mhd_p_aux s_p_aux;

static void _mrc_unused
pde_mhd_p_aux_setup_b0(struct mrc_fld *b0)
{
  if (b0) {
    fld3d_setup(&s_p_aux.b0, b0);
  }
}

static void _mrc_unused
pde_mhd_p_aux_setup_bnd_mask(struct mrc_fld *bnd_mask)
{
  if (bnd_mask) {
    fld3d_setup(&s_p_aux.bnd_mask, bnd_mask);
  }
}

static void _mrc_unused
pde_mhd_p_aux_get(int p)
{
  if (fld3d_is_setup(s_p_aux.b0)) {
    fld3d_get(&s_p_aux.b0, p);
  }
  if (fld3d_is_setup(s_p_aux.bnd_mask)) {
    fld3d_get(&s_p_aux.bnd_mask, p);
  }
}

static void _mrc_unused
pde_mhd_p_aux_put(int p)
{
  if (fld3d_is_setup(s_p_aux.b0)) {
    fld3d_put(&s_p_aux.b0, p);
  }
  if (fld3d_is_setup(s_p_aux.bnd_mask)) {
    fld3d_put(&s_p_aux.bnd_mask, p);
  }
}


// ======================================================================

// ----------------------------------------------------------------------
// EC_TO_CC
//
// edge-center to cell-center averaging

#if OPT_STAGGER == OPT_STAGGER_GGCM

#define EC_TO_CC(f, m, i,j,k)						\
  ({									\
    (.25f * (F3S(f, m, i-(m!=0?di:0),j            ,k            ) +	\
	     F3S(f, m, i            ,j-(m!=1?dj:0),k            ) +	\
	     F3S(f, m, i            ,j            ,k-(m!=2?dk:0)) +	\
	     F3S(f, m, i-(m!=0?di:0),j-(m!=1?dj:0),k-(m!=2?dk:0))));	\
  })

#else

#define EC_TO_CC(f, m, i,j,k)						\
  ({									\
    (.25f * (F3S(f, m, i+(m!=0?di:0),j+(m!=1?dj:0),k+(m!=2?dk:0)) +	\
	     F3S(f, m, i+(m!=0?di:0),j            ,k            ) +	\
	     F3S(f, m, i            ,j+(m!=1?dj:0),k            ) +	\
	     F3S(f, m, i            ,j            ,k+(m!=2?dk:0))));	\
  })

#endif

// ----------------------------------------------------------------------
// CC_TO_EC
//
// average (p_f, m) from cell center to edge center in direction M

#if OPT_STAGGER == OPT_STAGGER_GGCM

#define CC_TO_EC(p_f, m, i,j,k, M)					\
  (.25f * (F3S(p_f, m, i            ,j            ,k            ) +	\
	   F3S(p_f, m, i            ,j+(M!=1?dj:0),k+(M!=2?dk:0)) +	\
	   F3S(p_f, m, i+(M!=0?di:0),j            ,k+(M!=2?dk:0)) +	\
	   F3S(p_f, m, i+(M!=0?di:0),j+(M!=1?dj:0),k           )))
  
#else
  
#define CC_TO_EC(p_f, m, i,j,k, M)					\
  (.25f * (F3S(p_f, m, i-(M!=0?di:0),j-(M!=1?dj:0),k-(M!=2?dk:0)) +	\
	   F3S(p_f, m, i-(M!=0?di:0),j            ,k            ) +	\
	   F3S(p_f, m, i            ,j-(M!=1?dj:0),k            ) +	\
	   F3S(p_f, m, i            ,j            ,k-(M!=2?dk:0))))

#endif

// ----------------------------------------------------------------------
// CURL_FC
// 
// this curl gives values on face centers, assumes the input is on cell edges

#if OPT_STAGGER == OPT_STAGGER_GGCM

#define CURLX_FC(p_E, i,j,k) (PDE_INV_DY(j) * (F3S(p_E, 2, i,j,k) - F3S(p_E, 2, i,j-dj,k)) - \
			      PDE_INV_DZ(k) * (F3S(p_E, 1, i,j,k) - F3S(p_E, 1, i,j,k-dk)))
#define CURLY_FC(p_E, i,j,k) (PDE_INV_DZ(k) * (F3S(p_E, 0, i,j,k) - F3S(p_E, 0, i,j,k-dk)) - \
			      PDE_INV_DX(i) * (F3S(p_E, 2, i,j,k) - F3S(p_E, 2, i-di,j,k)))
#define CURLZ_FC(p_E, i,j,k) (PDE_INV_DX(i) * (F3S(p_E, 1, i,j,k) - F3S(p_E, 1, i-di,j,k)) - \
			      PDE_INV_DY(j) * (F3S(p_E, 0, i,j,k) - F3S(p_E, 0, i,j-dj,k)))

#else

#define CURLX_FC(p_E, i,j,k) (PDE_INV_DY(j) * (F3S(p_E, 2, i,j+dj,k) - F3S(p_E, 2, i,j,k)) - \
			      PDE_INV_DZ(k) * (F3S(p_E, 1, i,j,k+dk) - F3S(p_E, 1, i,j,k)))
#define CURLY_FC(p_E, i,j,k) (PDE_INV_DZ(k) * (F3S(p_E, 0, i,j,k+dk) - F3S(p_E, 0, i,j,k)) - \
			      PDE_INV_DX(i) * (F3S(p_E, 2, i+di,j,k) - F3S(p_E, 2, i,j,k)))
#define CURLZ_FC(p_E, i,j,k) (PDE_INV_DX(i) * (F3S(p_E, 1, i+di,j,k) - F3S(p_E, 1, i,j,k)) - \
			      PDE_INV_DY(j) * (F3S(p_E, 0, i,j+dj,k) - F3S(p_E, 0, i,j,k)))

#endif

#if OPT_STAGGER == OPT_STAGGER_GGCM
#define fld3d_foreach_stagger(i,j,k, l, r) fld3d_foreach(i,j,k, (l)+1, (r)-1)
#else
#define fld3d_foreach_stagger(i,j,k, l, r) fld3d_foreach(i,j,k, l, r)
#endif

#endif
