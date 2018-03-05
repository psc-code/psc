
#pragma once

#include "psc_bnd_particles_private.h"

#include <mrc_profile.h>

// ======================================================================
// BndParticlesBase

struct BndParticlesBase
{
  virtual void reset(struct mrc_domain*, const Grid_t& grid) = 0;
  virtual void exchange_particles(PscMparticlesBase mprts_base) = 0;
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
    static int pr;
    if (!pr) {
      pr = prof_register("xchg_prts", 1., 0, 0);
    }
    
    prof_start(pr);
    psc_stats_start(st_time_comm);
    sub()->exchange_particles(mprts);
    psc_stats_stop(st_time_comm);
    prof_stop(pr);
  }

  void reset()
  {
    if (!bndp_->obj.is_setup) return; // FIXME, hack around being called before constructed
    sub()->reset(bndp_->psc->mrc_domain, bndp_->psc->grid());
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
};

