
// ======================================================================
// options

struct mhd_options {
  int eqn;
  int limiter;
  int riemann;
  int time_integrator;
  int get_dt;
  bool background;
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
  { .val = OPT_LIMITER_GMINMOD  , .str = "gminmod"     },
  {},
};

#ifdef OPT_LIMITER
static const int s_opt_limiter = OPT_LIMITER;
#else
static int s_opt_limiter _mrc_unused;
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

  // riemann
#ifdef OPT_RIEMANN
  assert(OPT_RIEMANN == opt->riemann);
#else
  s_opt_riemann = opt->riemann;
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
}

// ======================================================================
// MHD parameters, we keep these around statically

static double s_gamma;  // adiabatic exponent

// ----------------------------------------------------------------------
// pde_mhd_setup

static void
pde_mhd_setup(struct ggcm_mhd *mhd)
{
  s_gamma = mhd->par.gamm;
}

// ======================================================================
// FIXME this is not a good place

#define BTX_(U, i,j,k, p)  (s_opt_background ? (BX_(U, i,j,k, p) + B0(b0, 0, i,j,k, p)) : BX_(U, i,j,k, p))
#define BTY_(U, i,j,k, p)  (s_opt_background ? (BY_(U, i,j,k, p) + B0(b0, 1, i,j,k, p)) : BY_(U, i,j,k, p))
#define BTZ_(U, i,j,k, p)  (s_opt_background ? (BZ_(U, i,j,k, p) + B0(b0, 2, i,j,k, p)) : BZ_(U, i,j,k, p))

