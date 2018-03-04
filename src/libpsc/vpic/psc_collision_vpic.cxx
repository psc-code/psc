
#include "psc_collision_private.h"
#include "collision.hxx"

#include "psc_particles_vpic.h"
#include "psc_method.h"

#include "vpic_iface.h"

// ----------------------------------------------------------------------

class PscCollisionVpic
{
public:
  constexpr static char const* const name = "vpic";

  PscCollisionVpic(MPI_Comm comm, int interval, double nu)
  {
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim_);
  }

  void run(psc_mparticles* mprts_base)
  {
    Simulation_collision_run(sim_);
  }

private:
  Simulation *sim_;
};

// ----------------------------------------------------------------------
// psc_collision: subclass "vpic"

psc_collision_ops_<PscCollisionVpic> psc_collision_vpic_ops;
