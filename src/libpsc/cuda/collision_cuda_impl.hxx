
#pragma once

#include "collision.hxx"
#include "psc_particles_cuda.h"

template<typename cuda_mparticles, typename RngState>
struct CudaCollision;

struct RngStateCuda;
struct RngStateFake;

// ----------------------------------------------------------------------
// CollisionCuda

template<typename BS, typename RngState = RngStateCuda>
class CollisionCuda : CollisionBase
{
public:
  using Mparticles = MparticlesCuda<BS>;
  
  CollisionCuda(MPI_Comm comm, int interval, double nu);
  
  virtual void run(PscMparticlesBase mprts_base) { assert(0); }

  void operator()(Mparticles& _mprts);

private:
  CudaCollision<cuda_mparticles<BS>, RngState> *fwd_;
};

