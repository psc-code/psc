
#pragma once

#include "collision.hxx"

template<typename cuda_mparticles, typename RngState>
struct CudaCollision;

struct RngStateCuda;
struct RngStateFake;

// ----------------------------------------------------------------------
// CollisionCuda

template<typename MP, typename RngState = RngStateCuda>
struct CollisionCuda : CollisionBase
{
  using Mparticles = MP;
  
  CollisionCuda(const Grid_t& grid, int interval, double nu);
  
  void operator()(Mparticles& _mprts);
  int interval() const;
  
private:
  CudaCollision<typename Mparticles::CudaMparticles, RngState> *fwd_;
};

