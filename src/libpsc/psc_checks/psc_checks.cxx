
#include "checks_impl.hxx"

#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "psc_fields_single.h"
#include "psc_fields_c.h"

// ----------------------------------------------------------------------
// psc_checks_init

struct psc_checks_1st_single_ops : psc_checks_ops {
  using Wrapper_t = ChecksWrapper<Checks_<MparticlesSingle, MfieldsSingle, checks_order_1st>>;
  psc_checks_1st_single_ops() {
    name                            = "1st_single";
    size                            = Wrapper_t::size;
    setup                           = Wrapper_t::setup;
    destroy                         = Wrapper_t::destroy;
    //read                            = psc_checks_sub_read;
  }
} psc_checks_1st_single_ops;

struct psc_checks_1st_double_ops : psc_checks_ops {
  using Wrapper_t = ChecksWrapper<Checks_<MparticlesDouble, MfieldsC, checks_order_1st>>;
  psc_checks_1st_double_ops() {
    name                            = "1st_double";
    size                            = Wrapper_t::size;
    setup                           = Wrapper_t::setup;
    destroy                         = Wrapper_t::destroy;
    //read                            = psc_checks_sub_read;
  }
} psc_checks_1st_double_ops;

struct psc_checks_2nd_single_ops : psc_checks_ops {
  using Wrapper_t = ChecksWrapper<Checks_<MparticlesSingle, MfieldsSingle, checks_order_2nd>>;
  psc_checks_2nd_single_ops() {
    name                            = "2nd_single";
    size                            = Wrapper_t::size;
    setup                           = Wrapper_t::setup;
    destroy                         = Wrapper_t::destroy;
    //read                            = psc_checks_sub_read;
  }
} psc_checks_2nd_single_ops;

struct psc_checks_2nd_double_ops : psc_checks_ops {
  using Wrapper_t = ChecksWrapper<Checks_<MparticlesDouble, MfieldsC, checks_order_2nd>>;
  psc_checks_2nd_double_ops() {
    name                            = "2nd_double";
    size                            = Wrapper_t::size;
    setup                           = Wrapper_t::setup;
    destroy                         = Wrapper_t::destroy;
    //read                            = psc_checks_sub_read;
  }
} psc_checks_2nd_double_ops;

static void
psc_checks_init()
{
  mrc_class_register_subclass(&mrc_class_psc_checks, &psc_checks_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_checks, &psc_checks_1st_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_checks, &psc_checks_2nd_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_checks, &psc_checks_2nd_single_ops);
}

// ----------------------------------------------------------------------
// psc_checks_descr

#define VAR(x) (void *)offsetof(struct psc_checks, params.x)

static struct param psc_checks_descr[] = {
  { "continuity_every_step" , VAR(continuity_every_step) , PARAM_INT(-1)       },
  { "continuity_threshold"  , VAR(continuity_threshold)  , PARAM_DOUBLE(1e-14) },
  { "continuity_verbose"    , VAR(continuity_verbose)    , PARAM_BOOL(false)   },
  { "continuity_dump_always", VAR(continuity_dump_always), PARAM_BOOL(false)   },

  { "gauss_every_step"      , VAR(gauss_every_step)      , PARAM_INT(-1)       },
  { "gauss_threshold"       , VAR(gauss_threshold)       , PARAM_DOUBLE(1e-14) },
  { "gauss_verbose"         , VAR(gauss_verbose)         , PARAM_BOOL(false)   },
  { "gauss_dump_always"     , VAR(gauss_dump_always)     , PARAM_BOOL(false)   },

  // { "rho_m"                 , VAR(rho_m)                 , MRC_VAR_OBJ(psc_mfields) },
  // { "rho_p"                 , VAR(rho_p)                 , MRC_VAR_OBJ(psc_mfields) },
  {},
};

#undef VAR

// ----------------------------------------------------------------------
// psc_checks class

struct mrc_class_psc_checks_ : mrc_class_psc_checks {
  mrc_class_psc_checks_() {
    name             = "psc_checks";
    size             = sizeof(struct psc_checks);
    param_descr      = psc_checks_descr;
    init             = psc_checks_init;
  }
} mrc_class_psc_checks;


