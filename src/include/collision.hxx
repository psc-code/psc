
#pragma once

#include "psc_collision_private.h"
#include "particles.hxx"

#include <mrc_profile.h>

// ======================================================================
// CollisionBase

class CollisionBase
{
public:
  virtual void run(MparticlesBase& mprts_base) = 0;
};

// ======================================================================
// PscCollision

template<typename S>
struct PscCollision
{
  using sub_t = S;
  
  static_assert(std::is_convertible<sub_t*, CollisionBase*>::value,
  		"sub classes used in PscCollision must derive from CollisionBase");
  
  explicit PscCollision(psc_collision *collision)
    : collision_(collision)
  {}

  void operator()(MparticlesBase& mprts)
  {
    static int st_time_collision;
    if (!st_time_collision) {
      st_time_collision = psc_stats_register("time collision");
    }
    
    psc_stats_start(st_time_collision);
    if (collision_->every > 0 && ppsc->timestep % collision_->every == 0) {
      sub()->run(mprts);
    }
    psc_stats_stop(st_time_collision);
  }
  
  sub_t* sub() { return mrc_to_subobj(collision_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_collision *collision_;
};

using PscCollisionBase = PscCollision<CollisionBase>;

// ======================================================================
// CollisionConvert

template<typename Collision_t>
struct CollisionConvert : Collision_t, CollisionBase
{
  using Base = Collision_t;
  using Mparticles = typename Collision_t::Mparticles;
  using Base::Base;

  void run(MparticlesBase& mprts_base) override
  {
    auto& mprts = mprts_base.get_as<Mparticles>();
    (*this)(mprts);
    mprts_base.put_as(mprts);
  }
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
  using Wrapper = CollisionWrapper<Collision>;
  psc_collision_ops_() {
    name                  = Wrapper::name;
    size                  = Wrapper::size;
    setup                 = Wrapper::setup;
    destroy               = Wrapper::destroy;
  }
};

