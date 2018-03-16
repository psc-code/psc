
#pragma once

#include "psc_collision_private.h"
#include "particles.hxx"

#include <mrc_profile.h>

// ======================================================================
// CollisionBase

class CollisionBase
{
public:
  virtual void run(PscMparticlesBase mprts_base) = 0;
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

  void operator()(PscMparticlesBase mprts)
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

template<typename Derived, typename MP>
struct CollisionCRTP : CollisionBase
{
  using Mparticles = MP;
  
  CollisionCRTP(int interval)
    : interval_(interval)
  {}

  // ----------------------------------------------------------------------
  // operator()
  
  void operator()(Mparticles& mprts)
  {
    static int pr;
    if (!pr) {
      pr = prof_register("collision", 1., 0, 0);
    }

    MPI_Comm comm = MPI_COMM_WORLD; // FIXME
    mpi_printf(comm, "***** Performing collisions...\n");
    prof_start(pr);
    static_cast<Derived*>(this)->collide(mprts);
    prof_stop(pr);
  }

  // ----------------------------------------------------------------------
  // run

  void run(PscMparticlesBase mprts_base) override
  {
    static int pr;
    if (!pr) {
      pr = prof_register("collision", 1., 0, 0);
    }

    auto mprts = mprts_base.get_as<PscMparticles<Mparticles>>();
    
    MPI_Comm comm = MPI_COMM_WORLD; // FIXME
    mpi_printf(comm, "***** Performing collisions...\n");
    prof_start(pr);
    static_cast<Derived*>(this)->collide(*mprts.sub());
    prof_stop(pr);
    
    mprts.put_as(mprts_base);
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

