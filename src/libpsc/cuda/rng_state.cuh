
#pragma once

#include <curand_kernel.h>

#define DEVICE __device__
#include "binary_collision.hxx" // for Rng, FIXME
#include "cuda_bits.h"

#include <thrust/device_vector.h>

// ======================================================================
// RngStateFake

class RngStateFake
{
public:
  // ======================================================================
  // RngStateFake::Device

  class Device
  {
  public:
    Device(RngStateFake& rng) : rng_(rng.rng_) {}

    __device__ RngFake operator[](int id) const { return rng_; }
    __device__ RngFake& operator[](int id) { return rng_; }

  private:
    RngFake rng_;
  };

  void resize(int size) { size_ = size; }
  int size() { return size_; }

private:
  RngFake rng_;
  int size_ = 0;
};

// ======================================================================
// RngStateCuda

class RngStateCuda
{
public:
  // ======================================================================
  // RngStateCuda::Rng

  class Rng
  {
  public:
    // ----------------------------------------------------------------------
    // uniform
    //
    // returns random number in ]0:1]

    __device__ float uniform() { return curand_uniform(&curand_state); }

    curandState curand_state;
  };

  // ======================================================================
  // RngStateCuda::Device

  class Device
  {
  public:
    Device(RngStateCuda& rng_state) : rngs_(rng_state.rngs_.data().get()), size_(rng_state.size()) {}

    __device__ Rng operator[](int id) const { return rngs_[id]; }
    __device__ Rng& operator[](int id) { return rngs_[id]; }

    __device__ int size() const { return size_; }

    Rng* rngs_;
    int size_;
  };

  RngStateCuda() = default;
  RngStateCuda(int size) { resize(size); }

  // also reseeds!
  void resize(int size);
  int size() { return rngs_.size(); }

private:
  psc::device_vector<Rng> rngs_;
};

// ----------------------------------------------------------------------
// k_curand_setup

#define THREADS_PER_BLOCK 256

__global__ static void k_curand_setup(RngStateCuda::Device rng_state, int size)
{
  int id = threadIdx.x + blockIdx.x * THREADS_PER_BLOCK;

  if (id < size) {
    curand_init(0, id, 0, &rng_state[id].curand_state);
  }
}

inline void RngStateCuda::resize(int size)
{
  rngs_.resize(size);

  dim3 dim_grid = (size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  k_curand_setup<<<dim_grid, THREADS_PER_BLOCK>>>(*this, size);
  cuda_sync_if_enabled();
};

#undef THREADS_PER_BLOCK
