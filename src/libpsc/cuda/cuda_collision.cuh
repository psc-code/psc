
#pragma once

#include "cuda_mparticles.h"
#include "cuda_mparticles_sort.cuh"

#define THREADS_PER_BLOCK 128

__global__ static void
k_collide(uint *d_off)
{
  uint end = d_off[blockIdx.x + 1];
  for (uint n = d_off[blockIdx.x] + threadIdx.x; n < end; n += THREADS_PER_BLOCK) {
    printf("%d/%d: n = %d\n", blockIdx.x, threadIdx.x, n);
  }
}

// ======================================================================
// cuda_collision

template<typename cuda_mparticles>
struct cuda_collision
{
  cuda_collision(int interval, double nu, int nicell, double dt)
    : interval_{interval}, nu_{nu}, nicell_(nicell), dt_(dt)
  {}
  
  void operator()(cuda_mparticles& cmprts)
  {
    auto sort_by_cell = cuda_mparticles_sort{cmprts.n_cells()};
    
    sort_by_cell.find_indices_ids(cmprts);
    sort_by_cell.stable_sort_cidx();
    sort_by_cell.find_offsets();

    dim3 dimGrid(cmprts.n_cells());
    k_collide<<<dimGrid, THREADS_PER_BLOCK>>>(sort_by_cell.d_off.data().get());
  }

private:
  int interval_;
  double nu_;
  int nicell_;
  double dt_;
};

