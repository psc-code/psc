
#pragma once

#include "psc_collision_private.h"
#include "particles.hxx"

#include <mrc_profile.h>

// ======================================================================
// PscCollision

template<typename S>
struct PscCollision
{
  using sub_t = S;
  
  explicit PscCollision(psc_collision *collision)
    : collision_(collision)
  {}

  void operator()(PscMparticlesBase mprts)
  {
    static int st_time_collision;
    if (!st_time_collision) {
      st_time_collision = psc_stats_register("time collision");
    }
    
    psc_stats_start(st_time_collision);
    (*this)->run(mprts);
    psc_stats_stop(st_time_collision);
  }
  
  sub_t* sub() { return mrc_to_subobj(collision_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_collision *collision_;
};

// ======================================================================
// CollisionBase

class CollisionBase
{
public:
  virtual void run(PscMparticlesBase mprts_base) = 0;
};

using PscCollisionBase = PscCollision<CollisionBase>;

template<typename Derived, typename MP>
struct CollisionCRTP : CollisionBase
{
  using Mparticles = MP;
  
  CollisionCRTP(int interval)
    : interval_(interval)
  {}

  // ----------------------------------------------------------------------
  // run

  void run(PscMparticlesBase mprts_base) override
  {
    static int pr;
    if (!pr) {
      pr = prof_register("collision", 1., 0, 0);
    }

    if (interval_ > 0 && ppsc->timestep % interval_ == 0) {
      auto mprts = mprts_base.get_as<PscMparticles<Mparticles>>();
      
      prof_start(pr);
      auto& derived = *static_cast<Derived*>(this);
      derived(mprts);
      prof_stop(pr);
      
      mprts.put_as(mprts_base);
    }
  }

protected:
  int interval_;
};

// ======================================================================
// CollisionWrapper

template<typename Collision>
class CollisionWrapper
{
public:
  const static size_t size = sizeof(Collision);

  constexpr static char const* const name = Collision::name;
  
  static void setup(struct psc_collision* _collision)
  {
    PscCollision<Collision> collision(_collision);
    
    new(collision.sub()) Collision(psc_collision_comm(_collision),
				   _collision->every, _collision->nu);
  }

  static void destroy(struct psc_collision* _collision)
  {
    PscCollision<Collision> collision(_collision);
    
    collision->~Collision();
  }
};

template<typename Collision>
struct psc_collision_ops_ : psc_collision_ops {
  using CollisionWrapper = CollisionWrapper<Collision>;
  psc_collision_ops_() {
    name                  = CollisionWrapper::name;
    size                  = CollisionWrapper::size;
    setup                 = CollisionWrapper::setup;
    destroy               = CollisionWrapper::destroy;
  }
};

