
#pragma once

#include "collision.hxx"
#include "vpic_iface.h"
#include "psc_method.h"

// ======================================================================
// PscCollisionVpic

class PscCollisionVpic : public CollisionBase
{
public:
  constexpr static char const* const name = "vpic";

  PscCollisionVpic(MPI_Comm comm, int interval, double nu)
  {
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim_);
  }

  void operator()(MparticlesBase& mprts_base) override
  {
    Simulation_collision_run(sim_);
  }

private:
  Simulation *sim_;
};

