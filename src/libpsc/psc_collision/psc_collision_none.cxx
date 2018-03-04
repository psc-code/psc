
#include "psc_collision_private.h"
#include "collision.hxx"

#include <mrc_profile.h>

// ======================================================================
// psc_collision_none

struct PscCollisionNone
{
  constexpr static char const* const name = "none";
  
  PscCollisionNone(MPI_Comm comm, int every, double nu)
  {}
  
  void run(struct psc_mparticles *mprts_base)
  {}
};

// ======================================================================
// psc_collision: subclass "none"

psc_collision_ops_<PscCollisionNone> psc_collision_none_ops;
