
#pragma once

#include "collision.hxx"
#include "psc_particles_cuda.h"
#include "psc_particles_single.h"
#include "psc_fields_single.h"

#include "../libpsc/psc_sort/psc_sort_impl.hxx"
#include "../libpsc/psc_collision/psc_collision_impl.hxx"

// ----------------------------------------------------------------------
// CollisionCudaHost

template<typename BS>
class CollisionCudaHost : CollisionBase
{
public:
  using Mparticles = MparticlesCuda<BS>;
  
  CollisionCudaHost(MPI_Comm comm, int interval, double nu);
  
  virtual void run(MparticlesBase& mprts_base) { assert(0); }

  void operator()(Mparticles& mprts)
  {
    auto& mprts = _mprts.template get_as<MparticlesSingle>();
    sort_(mprts);
    coll_(mprts);
    _mprts.put_as(mprts);
  }

private:
  SortCountsort2<MparticlesSingle> sort_;
  Collision_<MparticlesSingle, MfieldsSingle> coll_;
};

