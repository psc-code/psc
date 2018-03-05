
#pragma once

#include "psc_bnd_particles_private.h"

// ======================================================================
// BndParticlesBase

struct BndParticlesBase
{
  virtual void reset(struct mrc_domain*, const Grid_t& grid) = 0;
};

// ======================================================================
// PscBndParticles

template<typename S>
struct PscBndParticles
{
  using sub_t = S;
  
  static_assert(std::is_convertible<sub_t*, BndParticlesBase*>::value,
  		"sub classes used in PscBndParticles must derive from BndParticlesBase");
  
  explicit PscBndParticles(psc_bnd_particles *bndp)
    : bndp_(bndp)
  {}

  void operator()(PscMparticlesBase mprts)
  {
    psc_bnd_particles_exchange(bndp_, mprts.mprts());
  }
  
  sub_t* sub() { return mrc_to_subobj(bndp_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_bnd_particles* bndp_;
};

using PscBndParticlesBase = PscBndParticles<BndParticlesBase>;

// ======================================================================
// PscBndParticlesWrapper

template<typename BndParticles>
class PscBndParticlesWrapper
{
public:
  const static size_t size = sizeof(BndParticles);
  
  static void setup(psc_bnd_particles* _bndp)
  {
    PscBndParticles<BndParticles> bndp(_bndp);
    new(bndp.sub()) BndParticles(_bndp->psc->mrc_domain, _bndp->psc->grid());
  }

  static void destroy(psc_bnd_particles* _bndp)
  {
    PscBndParticles<BndParticles> bndp(_bndp);
    bndp->~BndParticles();
  }

  static void reset(psc_bnd_particles* _bndp)
  {
    PscBndParticles<BndParticles> bndp(_bndp);
    if (!bndp->ddcp) return; // FIXME, hack around being called before constructed
    
    bndp->reset(_bndp->psc->mrc_domain, _bndp->psc->grid());
  }

  static void exchange_particles(struct psc_bnd_particles *_bndp,
				 struct psc_mparticles *mprts_base)
  {
    PscBndParticles<BndParticles> bndp(_bndp);
    
    bndp->exchange_particles(mprts_base);
  }
  
};

