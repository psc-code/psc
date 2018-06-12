
#pragma once

#include "collision.hxx"
#include "psc_particles_cuda.h"
#include "psc_particles_single.h"
#include "psc_fields_single.h"

#include "../libpsc/psc_sort/psc_sort_impl.hxx"
#include "../libpsc/psc_collision/psc_collision_impl.hxx"

template<typename cuda_mparticles>
struct cuda_collision;

// ----------------------------------------------------------------------
// CollisionCuda

template<typename BS>
class CollisionCuda : CollisionBase
{
public:
  using Mparticles = MparticlesCuda<BS>;
  
  CollisionCuda(MPI_Comm comm, int interval, double nu);
  
  virtual void run(PscMparticlesBase mprts_base) { assert(0); }

  void operator()(Mparticles& _mprts);

private:
  SortCountsort2<MparticlesSingle> sort_;
  Collision_<MparticlesSingle, MfieldsSingle> coll_;
  cuda_collision<cuda_mparticles<BS>> *fwd_;
};

