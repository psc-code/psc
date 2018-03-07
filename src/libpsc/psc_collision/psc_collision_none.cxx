
#include "collision.hxx"

struct CollisionNone : CollisionBase
{
  constexpr static const char* name = "none";

  CollisionNone(MPI_Comm comm, int interval, double nu) {}

  virtual void run(PscMparticlesBase mprts_base) {}
};

psc_collision_ops_<CollisionNone> psc_collision_none_ops;
