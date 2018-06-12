
#pragma once

#include "collision.hxx"
#include "psc_particles_cuda.h"

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
  cuda_collision<cuda_mparticles<BS>> *fwd_;
};

