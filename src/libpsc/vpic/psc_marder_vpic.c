
#include "psc_marder_vpic.h"

#include "psc_fields_vpic.h"
#include "psc_particles_vpic.h"
#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_marder_vpic_setup

static void
psc_marder_vpic_setup(struct psc_marder *marder)
{
  struct psc_marder_vpic *sub = psc_marder_vpic(marder);

  struct vpic_marder_info info = {
    .clean_div_e_interval = marder->clean_div_e_interval,
    .clean_div_b_interval = marder->clean_div_b_interval,
    .sync_shared_interval = marder->sync_shared_interval,
    .num_div_e_round = marder->num_div_e_round,
    .num_div_b_round = marder->num_div_b_round,
  };
  sub->vmarder = vpic_marder_create();
  vpic_marder_ctor(sub->vmarder, &info);
}

// ----------------------------------------------------------------------
// psc_marder_vpic_run

static void
psc_marder_vpic_run(struct psc_marder *marder,
		    struct psc_mfields *mflds_base,
		    struct psc_mparticles *mprts_base)
{
  struct vpic_marder *vmarder = psc_marder_vpic(marder)->vmarder;
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds_base)->vmflds;
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(mprts_base)->vmprts;
    
  struct psc *psc = ppsc; // FIXME

  vpic_marder_run(vmarder, vmflds, vmprts, psc->timestep);
}

// ----------------------------------------------------------------------
// psc_marder: subclass "vpic"

struct psc_marder_ops psc_marder_vpic_ops = {
  .name                  = "vpic",
  .size                  = sizeof(struct psc_marder_vpic),
  .setup                 = psc_marder_vpic_setup,
  .run                   = psc_marder_vpic_run,
};

