
#pragma once

#include "collision.hxx"
#include "vpic_iface.h"

// ======================================================================
// PscCollisionVpic

class PscCollisionVpic : public CollisionBase
{
public:
  constexpr static char const* const name = "vpic";

  PscCollisionVpic(const Grid_t& grid, int interval, double nu)
    : interval_{interval}
  {}

  void operator()(MparticlesBase& mprts_base) override
  {
#if 0
    // Note: Particles should not have moved since the last performance sort
    // when calling collision operators.
    // FIXME: Technically, this placement of the collision operators only
    // yields a first order accurate Trotter factorization (not a second
    // order accurate factorization).
    
    if (collision_op_list) {
      // FIXME: originally, vpic_clear_accumulator_array() was called before this.
      // It's now called later, though. I'm not sure why that would be necessary here,
      // but it needs to be checked.
      // The assert() below doesn't unfortunately catch all cases where this might go wrong
      // (ie., it's missing the user_particle_collisions())
      
      assert(0);
      TIC ::apply_collision_op_list(collision_op_list); TOC(collision_model, 1);
    }
    TIC user_particle_collisions(); TOC(user_particle_collisions, 1);
#endif
  }

  int interval() const { return interval_; }
  
private:
  int interval_;
};

