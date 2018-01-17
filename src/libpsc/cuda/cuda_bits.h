
#ifndef CUDA_BITS_H
#define CUDA_BITS_H

#define cudaCheck(ierr) do {						\
    if (ierr != cudaSuccess)						\
      fprintf(stderr, "IERR = %d (%s)\n", ierr, cudaGetErrorName(ierr)); \
    assert(ierr == cudaSuccess);					\
  } while(0)

static bool CUDA_SYNC = true;

#define cuda_sync_if_enabled() do {					\
    cudaError ierr = cudaGetLastError(); cudaCheck(ierr);		\
    if (CUDA_SYNC) {							\
      cudaError_t ierr = cudaDeviceSynchronize(); cudaCheck(ierr);	\
    }									\
  } while(0)

static inline float cuda_int_as_float(int i)
{
  union { int i; float f; } u;
  u.i = i;
  return u.f;
};

static inline int cuda_float_as_int(float f)
{
  union { int i; float f; } u;
  u.f = f;
  return u.i;
};

#endif
