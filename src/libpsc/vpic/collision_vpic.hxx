
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
    : interval_{interval}
  {
    
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim_);
  }

  void operator()(MparticlesBase& mprts_base) override
  {
    sim_->collision_run();
  }

  int interval() const { return interval_; }
  
private:
  int interval_;
  Simulation *sim_;
};

