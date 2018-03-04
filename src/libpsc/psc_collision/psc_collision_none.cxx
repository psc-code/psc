
#include "psc_collision_private.h"
#include "collision.hxx"

#include <mrc_profile.h>

// ======================================================================
// psc_collision_none

struct PscCollisionNone
{
  PscCollisionNone(MPI_Comm comm, int every, double nu)
  {}
  
  void run(struct psc_mparticles *mprts_base)
  {}
};

// ======================================================================
// psc_collision: subclass "none"

struct psc_collision_none_ops : psc_collision_ops {
  using Collision = CollisionWrapper<PscCollisionNone>;
  psc_collision_none_ops() {
    name                  = "none";
    size                  = Collision::size;
    setup                 = Collision::setup;
    destroy               = Collision::destroy;
    run                   = Collision::run;
  }
} psc_collision_none_ops;

