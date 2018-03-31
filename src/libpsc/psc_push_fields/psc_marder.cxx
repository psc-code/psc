
#include "psc_marder_private.h"
#include "psc_bnd.h"
#include "psc_output_fields_item.h"
#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "psc_fields_single.h"
#include "fields3d.hxx"
#include "bnd.hxx"
#include "fields_item.hxx"

#include <mrc_io.h>
#include <mrc_profile.h>

#include <math.h>

// ----------------------------------------------------------------------
// marder_calc_aid_fields

static void
marder_calc_aid_fields(struct psc_marder *marder, 
		       PscMfieldsBase mflds_base, PscMparticlesBase mprts_base)
{
  PscFieldsItemBase item_div_e(marder->item_div_e);
  PscFieldsItemBase item_rho(marder->item_rho);
  item_div_e(mflds_base, mprts_base); // FIXME, should accept NULL for particles
  
  if (marder->dump) {
    static int cnt;
    mrc_io_open(marder->io, "w", cnt, cnt);//ppsc->timestep, ppsc->timestep * ppsc->dt);
    cnt++;
    psc_mfields_write_as_mrc_fld(item_rho->mres().mflds(), marder->io);
    psc_mfields_write_as_mrc_fld(item_div_e->mres().mflds(), marder->io);
    mrc_io_close(marder->io);
  }

  item_div_e->mres()->axpy_comp(0, -1., *item_rho->mres().sub(), 0);
  // FIXME, why is this necessary?
  auto bnd = PscBndBase(marder->bnd);
  bnd.fill_ghosts(item_div_e->mres(), 0, 1);
}

// ----------------------------------------------------------------------
// psc_marder_init

#include "marder_impl.hxx"

using marder_ops_c = marder_ops<MparticlesDouble, MfieldsC>;

static struct psc_marder_ops_c : psc_marder_ops {
  psc_marder_ops_c() {
    name                  = "c";
    setup                 = marder_ops_c::setup;
    destroy               = marder_ops_c::destroy;
    run                   = marder_ops_c::run;
  }
} psc_marder_c_ops;

using marder_ops_single = marder_ops<MparticlesSingle, MfieldsSingle>;

static struct psc_marder_ops_single : psc_marder_ops {
  psc_marder_ops_single() {
    name                  = "single";
    setup                 = marder_ops_single::setup;
    destroy               = marder_ops_single::destroy;
    run                   = marder_ops_single::run;
  }
} psc_marder_single_ops;

extern struct psc_marder_ops psc_marder_cuda_ops;
extern struct psc_marder_ops psc_marder_vpic_ops;

static void
psc_marder_init()
{
  mrc_class_register_subclass(&mrc_class_psc_marder, &psc_marder_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_marder, &psc_marder_single_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_marder, &psc_marder_cuda_ops);
#endif
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_marder, &psc_marder_vpic_ops);
#endif
}

// ----------------------------------------------------------------------
// psc_marder description

#define VAR(x) (void *)offsetof(struct psc_marder, x)
static struct param psc_marder_descr[] = {
  { "every_step"       , VAR(every_step)       , PARAM_INT(-1),     },
  { "diffusion"        , VAR(diffusion)        , PARAM_DOUBLE(0.9), },
  { "loop"             , VAR(loop)             , PARAM_INT(1),      },
  { "dump"             , VAR(dump)             , PARAM_BOOL(false), },

  { "clean_div_e_interval", VAR(clean_div_e_interval), PARAM_INT(0),     },
  { "clean_div_b_interval", VAR(clean_div_b_interval), PARAM_INT(0),     },
  { "sync_shared_interval", VAR(sync_shared_interval), PARAM_INT(0),     },
  { "num_div_e_round"     , VAR(num_div_e_round)     , PARAM_INT(0),     },
  { "num_div_b_round"     , VAR(num_div_b_round)     , PARAM_INT(0),     },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_marder class description

struct mrc_class_psc_marder_ : mrc_class_psc_marder {
  mrc_class_psc_marder_() {
    name             = "psc_marder";
    size             = sizeof(struct psc_marder);
    param_descr      = psc_marder_descr;
    init             = psc_marder_init;
  }
} mrc_class_psc_marder;

