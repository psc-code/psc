
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

  using value_type = Device;

  void resize(int size) { size_ = size; }
  int size() const { return size_; }
  int capacity() const { return size_; }

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
    __device__ float normal() { return curand_normal(&curand_state); }
    __device__ float2 normal2() { return curand_normal2(&curand_state); }

    curandState curand_state;
  };

  using value_type = Rng;

  // ======================================================================
  // RngStateCuda::Device

  class Device
  {
  public:
    Device(RngStateCuda& rng_state)
      : rngs_(rng_state.rngs_.data().get()), size_(rng_state.size())
    {}

    __device__ Rng operator[](int id) const { return rngs_[id]; }
    __device__ Rng& operator[](int id) { return rngs_[id]; }

    __device__ int size() const { return size_; }

    Rng* rngs_;
    int size_;
  };

  RngStateCuda() = default;
  RngStateCuda(const RngStateCuda&) = delete;
  RngStateCuda(int size) { resize(size); }

  ~RngStateCuda() { mem_rnd -= allocated_bytes(rngs_); }

  // also reseeds!
  void resize(int size);
  int size() const { return rngs_.size(); }
  int capacity() const { return rngs_.capacity(); }

private:
  psc::device_vector<Rng> rngs_;
};

// ----------------------------------------------------------------------
// k_curand_setup

__global__ static void k_curand_setup(RngStateCuda::Device rng_state, int size)
{
  int id = threadIdx.x + blockIdx.x * blockDim.x;

  if (id < size) {
    curand_init(0, id, 0, &rng_state[id].curand_state);
  }
}

inline void RngStateCuda::resize(int size)
{
  if (size > rngs_.size()) {
    const int threads_per_block = 256;

    std::cout << "RngState resize " << size << "\n";
    mem_rnd -= allocated_bytes(rngs_);
    rngs_.resize(size);
    mem_rnd += allocated_bytes(rngs_);

    dim3 dim_grid = (size + threads_per_block - 1) / threads_per_block;
    k_curand_setup<<<dim_grid, threads_per_block>>>(*this, size);
    cuda_sync_if_enabled();
  }
};
