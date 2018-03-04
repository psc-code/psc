
#pragma once

#include "psc_collision_private.h"
#include "particles.hxx"

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

    struct psc_collision_ops *ops = psc_collision_ops(collision_);
    assert(ops->run);
    ops->run(collision_, mprts.mprts());
    
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
  virtual void run(struct psc_mparticles *mprts_base) = 0;
};

using PscCollisionBase = PscCollision<CollisionBase>;

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

  static void run(struct psc_collision* _collision, struct psc_mparticles* mprts)
  {
    PscCollision<Collision> collision(_collision);

    collision->run(mprts);
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
    run                   = CollisionWrapper::run;
  }
};

