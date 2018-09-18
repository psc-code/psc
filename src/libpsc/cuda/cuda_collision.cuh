
#pragma once

#include "cuda_mparticles.h"
#include "cuda_mparticles_sort.cuh"
// FIXME, horrible hack...
#define DEVICE __device__
#include "binary_collision.hxx"

#include <curand_kernel.h>

#define THREADS_PER_BLOCK 128

template<typename cuda_mparticles, typename RngState>
struct CudaCollision;

// ======================================================================
// RngStateFake

struct RngStateFake
{
  void init(dim3 dim_grid) {}

  __device__
  RngFake  operator[](int id) const { return rng_; }

  __device__
  RngFake& operator[](int id)       { return rng_; }

private:
  RngFake rng_;
};
  
// ======================================================================
// RngStateCuda

struct RngStateCuda
{
  // ======================================================================
  // RngStateCuda::Rng

  struct Rng
  {
    // ----------------------------------------------------------------------
    // uniform
    //
    // returns random number in ]0:1]
    
    __device__
    float uniform()
    {
      return curand_uniform(&curand_state);
    }
    
    curandState curand_state;
  };

  void init(dim3 dim_grid);

  __device__
  Rng  operator[](int id) const { return rngs_[id]; }

  __device__
  Rng& operator[](int id)       { return rngs_[id]; }

  Rng* rngs_;
};

// ----------------------------------------------------------------------
// k_curand_setup

__global__
static void k_curand_setup(RngStateCuda rng_state)
{
  int id = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;

  curand_init(1234, id % 1024, 0, &rng_state[id].curand_state); // FIXME, % 1024 hack
}

void RngStateCuda::init(dim3 dim_grid)
{
  int n_threads = dim_grid.x * THREADS_PER_BLOCK;
  
  cudaError_t ierr;
  ierr = cudaMalloc(&rngs_, n_threads * sizeof(*rngs_));
  cudaCheck(ierr);
  
  k_curand_setup<<<dim_grid, THREADS_PER_BLOCK>>>(*this);
  cuda_sync_if_enabled();
};

template<typename cuda_mparticles, typename RngState>
__global__ static void
k_collide(DMparticlesCuda<typename cuda_mparticles::BS> dmprts, uint* d_off, uint* d_id, float nudt0,
	  RngState rng_state, uint n_cells)
{
  CudaCollision<cuda_mparticles, RngState>::d_collide(dmprts, d_off, d_id, nudt0, rng_state, n_cells);
}

// ======================================================================
// cuda_collision

template<typename cuda_mparticles, typename RngState>
struct CudaCollision
{
  using real_t = typename cuda_mparticles::real_t;
  using DMparticles = DMparticlesCuda<typename cuda_mparticles::BS>;
  
  struct Particle
  {
    using real_t = real_t;

    __device__
    Particle(DMparticles& dmprts, int n)
      : dmprts_{dmprts},
	n_{n}
    {
      LOAD_PARTICLE_POS(prt_, dmprts_.xi4_ , n_);
      LOAD_PARTICLE_MOM(prt_, dmprts_.pxi4_, n_);
    }

    __device__
    ~Particle()
    {
      // xi4 is not modified
      STORE_PARTICLE_MOM(prt_, dmprts_.pxi4_, n_);
    }
    
    __device__
    real_t q() const
    {
      int kind = __float_as_int(prt_.kind_as_float);
      return dmprts_.q(kind);
    }

    __device__
    real_t m() const
    {
      int kind = __float_as_int(prt_.kind_as_float);
      return dmprts_.m(kind);
    }
    
    __device__
    real_t u(int d) const
    {
      return prt_.pxi[d];
    }

    __device__
    real_t& u(int d)
    {
      return prt_.pxi[d];
    }

  private:
    DMparticles& dmprts_;
    int n_;
    d_particle prt_;
  };
  
  CudaCollision(int interval, double nu, int nicell, double dt)
    : interval_{interval}, nu_{nu}, nicell_(nicell), dt_(dt)
  {}

  int interval() const
  {
    return interval_;
  }
  
  void operator()(cuda_mparticles& cmprts)
  {
    auto sort_by_cell = cuda_mparticles_sort{cmprts.n_cells()};
    
    sort_by_cell.find_indices_ids(cmprts);
    sort_by_cell.stable_sort_cidx();
    sort_by_cell.find_offsets();
    // for (int c = 0; c <= cmprts.n_cells(); c++) {
    //   printf("off[%d] = %d\n", c, int(sort_by_cell.d_off[c]));
    // }

    int blocks = cmprts.n_cells();
    if (blocks > 32768) blocks = 32768;
    dim3 dimGrid(blocks);

    static bool first_time = true;
    if (first_time) {
      rng_state_.init(dimGrid);
      first_time = false;
    }
    
    // all particles need to have same weight!
    real_t wni = 1.; // FIXME, there should at least be some assert to enforce this //prts.prt_wni(prts[n_start]);
    real_t nudt0 = wni / nicell_ * interval_ * dt_ * nu_;

    k_collide<cuda_mparticles><<<dimGrid, THREADS_PER_BLOCK>>>(cmprts, sort_by_cell.d_off.data().get(),
					      sort_by_cell.d_id.data().get(),
							       nudt0, rng_state_,
					      cmprts.n_cells());
    cuda_sync_if_enabled();
  }

  __device__
  static void d_collide(DMparticles dmprts, uint* d_off, uint* d_id, float nudt0,
			RngState rng_state, uint n_cells)
  {
    
    int id = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;
    /* Copy state to local memory for efficiency */
    auto rng = rng_state[id];
    BinaryCollision<Particle> bc;
    
    for (uint bidx = blockIdx.x; bidx < n_cells; bidx += gridDim.x) {
      uint beg = d_off[bidx];
      uint end = d_off[bidx + 1];
      real_t nudt1 = nudt0 * (end - beg & ~1); // somewhat counteract that we don't collide the last particle if odd
      for (uint n = beg + 2*threadIdx.x; n + 1 < end; n += 2*THREADS_PER_BLOCK) {
	//printf("%d/%d: n = %d off %d\n", blockIdx.x, threadIdx.x, n, d_off[blockIdx.x]);
	Particle prt1{dmprts, int(d_id[n  ])};
	Particle prt2{dmprts, int(d_id[n+1])};
	bc(prt1, prt2, nudt1, rng);
      }
    }
    
    rng_state[id] = rng;
  }

private:
  int interval_;
  double nu_;
  int nicell_;
  double dt_;
  RngState rng_state_;
};

