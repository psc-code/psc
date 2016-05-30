
// ======================================================================
// options

struct mhd_options {
  int riemann;
};

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
// pde_mhd_set_options

static void _mrc_unused
pde_mhd_set_options(struct ggcm_mhd *mhd, struct mhd_options *opt)
{
  // riemann
#ifdef OPT_RIEMANN
  assert(OPT_RIEMANN == opt->riemann);
#else
  s_opt_riemann = opt->riemann;
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

