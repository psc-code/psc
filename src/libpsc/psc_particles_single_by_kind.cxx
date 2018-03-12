
#include "bk_mparticles_iface.h"

#include "psc.h"
#include "psc_particles_single_by_kind.h"

#define PARTICLE_TYPE "single_by_kind"
#define PFX(x) psc_mparticles_single_by_kind_ ## x
#define psc_mparticles_sub psc_mparticles_single_by_kind

// ======================================================================
// psc_mparticles: subclass "single_by_kind"

static struct mrc_obj_method psc_mparticles_single_by_kind_methods[] = {
  {}
};

// ----------------------------------------------------------------------
// psc_mparticles_sub_setup

static void
PFX(setup)(struct psc_mparticles *_mprts)
{
  PscMparticlesSingleByKind mprts(_mprts);

  new(mprts.sub()) MparticlesSingleByKind{ppsc->grid()};
}

static void
PFX(destroy)(struct psc_mparticles *_mprts)
{
  auto mprts = PscMparticlesSingleByKind{_mprts};

  bk_mparticles_delete(mprts->bkmprts);
}

static unsigned int
PFX(get_nr_particles)(struct psc_mparticles *_mprts)
{
  auto mprts = PscMparticlesSingleByKind{_mprts};

  return bk_mparticles_n_prts(mprts->bkmprts);
}

// ----------------------------------------------------------------------
// psc_mparticles_ops

struct PFX(ops) : psc_mparticles_ops {
  using Wrapper_t = MparticlesWrapper<MparticlesSingleByKind>;
  PFX(ops)() {
    name                    = Wrapper_t::name;
    size                    = Wrapper_t::size;
    methods                 = PFX(methods);
    setup                   = Wrapper_t::setup;
    destroy                 = Wrapper_t::destroy;
#if 0
    write                   = PFX(write);
    read                    = PFX(read);
#endif
  }
} PFX(ops);

