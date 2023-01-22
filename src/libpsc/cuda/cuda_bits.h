
#ifndef CUDA_BITS_H
#define CUDA_BITS_H

#include "PscConfig.h"

#include "gpu/gpu_runtime.h"

#include <iostream>
#include <fstream>

#ifdef PSC_HAVE_RMM
#include <rmm/device_vector.hpp>
#include <rmm/exec_policy.hpp>

namespace psc
{
template <typename T>
using device_vector = rmm::device_vector<T>;
}

#else
#include <thrust/device_vector.h>

namespace psc
{
template <typename T>
using device_vector = thrust::device_vector<T>;
}
#endif

#define cudaCheck(ierr)                                                        \
  do {                                                                         \
    if (ierr != hipSuccess) {                                                  \
      fprintf(stderr, "IERR = %d (%s) %s:%d\n", ierr, hipGetErrorName(ierr),   \
              __FILE__, __LINE__);                                             \
    }                                                                          \
    if (ierr != hipSuccess)                                                    \
      abort();                                                                 \
  } while (0)

static bool CUDA_SYNC = true;

#define cuda_sync_if_enabled()                                                 \
  do {                                                                         \
    hipError_t ierr = hipGetLastError();                                       \
    cudaCheck(ierr);                                                           \
    if (CUDA_SYNC) {                                                           \
      hipError_t ierr = hipDeviceSynchronize();                                \
      cudaCheck(ierr);                                                         \
    }                                                                          \
  } while (0)

__host__ __device__ static inline float cuda_int_as_float(int i)
{
  union
  {
    int i;
    float f;
  } u;
  u.i = i;
  return u.f;
};

__host__ __device__ static inline int cuda_float_as_int(float f)
{
  union
  {
    int i;
    float f;
  } u;
  u.f = f;
  return u.i;
};

extern std::size_t mem_particles;
extern std::size_t mem_collisions;
extern std::size_t mem_randomize_sort;
extern std::size_t mem_sort_by_block;
extern std::size_t mem_bnd;
extern std::size_t mem_heating;
extern std::size_t mem_bndp;
extern std::size_t mem_rnd;

template <typename V>
std::size_t allocated_bytes(const V& v)
{
  return v.capacity() * sizeof(typename V::value_type);
}

void mem_stats_csv_header(std::ostream& of);
void mem_stats_csv(std::ostream& of, int timestep, int n_patches, int n_prts);

void mem_stats(std::string file, int line, std::ostream& of);

void mem_pool_print();

#endif
