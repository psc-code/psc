
#include "psc_particles_vpic.h"

#include "vpic_iface.h"

// ======================================================================
// conversion to "single"

static void
psc_mparticles_vpic_copy_from_single(struct psc_mparticles *mprts,
				    struct psc_mparticles *mprts_single, unsigned int flags)
{
}

static void
psc_mparticles_vpic_copy_to_single(struct psc_mparticles *mprts,
				  struct psc_mparticles *mprts_single, unsigned int flags)
{
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_methods

static struct mrc_obj_method psc_mparticles_vpic_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_mparticles_vpic_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_mparticles_vpic_copy_from_single),
  {}
};

// ----------------------------------------------------------------------
// psc_mparticles_vpic_setup

static void
psc_mparticles_vpic_setup(struct psc_mparticles *mprts)
{
  struct psc_mparticles_vpic *sub = psc_mparticles_vpic(mprts);

  sub->vmprts = vpic_mparticles_create();
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_reserve_all

static void
psc_mparticles_vpic_reserve_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_get_size_all

static void
psc_mparticles_vpic_get_size_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_resize_all

static void
psc_mparticles_vpic_resize_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
}

// ----------------------------------------------------------------------
// psc_mparticles_vpic_get_nr_particles

static unsigned int
psc_mparticles_vpic_get_nr_particles(struct psc_mparticles *mprts)
{
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(mprts)->vmprts;
  
  return vpic_mparticles_get_nr_particles(vmprts);
}

// ----------------------------------------------------------------------
// psc_mparticles: subclass "vpic"
  
struct psc_mparticles_ops psc_mparticles_vpic_ops = {
  .name                    = "vpic",
  .size                    = sizeof(struct psc_mparticles_vpic),
  .methods                 = psc_mparticles_vpic_methods,
  .setup                   = psc_mparticles_vpic_setup,
  .reserve_all             = psc_mparticles_vpic_reserve_all,
  .get_size_all            = psc_mparticles_vpic_get_size_all,
  .resize_all              = psc_mparticles_vpic_resize_all,
  .get_nr_particles        = psc_mparticles_vpic_get_nr_particles,
};

