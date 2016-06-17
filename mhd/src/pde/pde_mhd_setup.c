
// ======================================================================
// MHD parameters, we keep these around statically

static mrc_fld_data_t s_gamma;  // adiabatic exponent
static mrc_fld_data_t s_eta;    // (constant) resistivity 
static mrc_fld_data_t s_d_i;    // ion skin depth
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
pde_mhd_setup(struct ggcm_mhd *mhd)
{
  s_gamma = mhd->par.gamm;
  s_eta   = mhd->par.diffco;
  s_d_i   = mhd->par.d_i;
}

// ======================================================================
// FIXME this is not a good place

#define BTX_(U, i,j,k, p)  (s_opt_background ? (BX_(U, i,j,k, p) + B0(b0, 0, i,j,k, p)) : BX_(U, i,j,k, p))
#define BTY_(U, i,j,k, p)  (s_opt_background ? (BY_(U, i,j,k, p) + B0(b0, 1, i,j,k, p)) : BY_(U, i,j,k, p))
#define BTZ_(U, i,j,k, p)  (s_opt_background ? (BZ_(U, i,j,k, p) + B0(b0, 2, i,j,k, p)) : BZ_(U, i,j,k, p))

// ======================================================================
// mhd auxiliary fields
//
// kept around statically so we don't have to pass all this crap

struct mhd_aux {
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
}
