
#ifndef CACHE_SHAPE_ARRAYS
#define CACHE_SHAPE_ARRAYS 2
#endif

// ----------------------------------------------------------------------

__shared__ real sdata[(2*SW + 1) * THREADS_PER_BLOCK];

#define SDATA_(sdata, tid, j) (sdata[((j)+SW) * THREADS_PER_BLOCK + (tid)])
#define SDATA(tid,j) SDATA_(sdata, tid, j)

__device__ static real
find_shape_coeff_d_shift(int j, real h, short int shift)
{
  if (j ==  0 + shift) {
    return real(.75) - sqr(h);
  } else if (j == -1 + shift) {
    return real(.5) * sqr(real(1.5)-fabsr(h-real(1.)));
  } else if (j == +1 + shift) {
    return real(.5) * sqr(real(1.5)-fabsr(h+real(1.)));
  } else {
    return real(0.);
  }
}

// ----------------------------------------------------------------------
#if CACHE_SHAPE_ARRAYS == 1 // broken

struct shape_info {
  real s0[2*5], s1[2*5];
};

#define s0(si, d, j) ARR_3_off2((si)->s0, d, j)
#define s1(si, d, j) ARR_3_off2((si)->s1, d, j)

__device__ static void
cache_shape_arrays(struct shape_info *si, real *h0, real *h1, short int *shift)
{
  find_shape_coeff_d(si->s0, 1, h0[1], 0);
  find_shape_coeff_d(si->s0, 2, h0[2], 0);

  find_shape_coeff_d(si->s1, 1, h1[1], shift[1]);
  find_shape_coeff_d(si->s1, 2, h1[2], shift[2]);

  for (int d = 1; d < 3; d++) {
    for (int i = -2; i <= 2; i++) {
      s1(si, d, i) -= s0(si, d, i);
    }
  }
}

// ----------------------------------------------------------------------
#elif CACHE_SHAPE_ARRAYS == 2

#if DIM == DIM_Z

#define DECLARE_SHAPE_INFO			\
  real __s0z[2], __s1z[2];			\
  short int __shift0[1], __shift1[1]		\

#define SHAPE_INFO_PARAMS __s0z, __s1z, __shift0, __shift1
#define SHAPE_INFO_ARGS real *__s0z, real *__s1z, short int *__shift0, short int *__shift1

#define SI_SHIFT0Z __shift0[0]
#define SI_SHIFT1Z __shift1[0]
#define SI_SHIFT10Z (__shift1[0] - __shift0[0])

#elif DIM == DIM_YZ

#define DECLARE_SHAPE_INFO			\
  real __s0y[2], __s1y[2];			\
  real __s0z[2], __s1z[2];			\
  short int __shift0[2], __shift1[2]		\

#define SHAPE_INFO_PARAMS __s0y, __s1y, __s0z, __s1z, __shift0, __shift1
#define SHAPE_INFO_ARGS real *__s0y, real *__s1y, real *__s0z, real *__s1z, short int *__shift0, short int *__shift1

#define SI_SHIFT0Y __shift0[0]
#define SI_SHIFT1Y __shift1[0]
#define SI_SHIFT10Y (__shift1[0] - __shift0[0])
#define SI_SHIFT0Z __shift0[1]
#define SI_SHIFT1Z __shift1[1]
#define SI_SHIFT10Z (__shift1[1] - __shift0[1])

#endif

__device__ static inline void
calc_shape_coeff(real *__s, real h)
{
  __s[0] = find_shape_coeff_d_shift(-1, h, 0);
  __s[1] = find_shape_coeff_d_shift( 0, h, 0);
}

#if DIM == DIM_Z

__device__ static void
cache_shape_arrays(SHAPE_INFO_ARGS, real *h0, real *h1,
		   short int shift0z, short int shift1z)
{
  __shift0[0] = shift0z;
  __shift1[0] = shift1z;
  calc_shape_coeff(__s0z, h0[2]);
  calc_shape_coeff(__s1z, h1[2]);
}

#elif DIM == DIM_YZ

__device__ static void
cache_shape_arrays(SHAPE_INFO_ARGS, real *h0, real *h1,
		   short int shift0y, short int shift0z,
		   short int shift1y, short int shift1z)
{
  __shift0[0] = shift0y;
  __shift1[0] = shift1y;
  __shift0[1] = shift0z;
  __shift1[1] = shift1z;
  calc_shape_coeff(__s0y, h0[1]);
  calc_shape_coeff(__s1y, h1[1]);
  calc_shape_coeff(__s0z, h0[2]);
  calc_shape_coeff(__s1z, h1[2]);
}

#endif

#define pick_shape_coeff(t, comp, j, shift)	  \
  __pick_shape_coeff(__s ## t ## comp, j, shift)  \

__device__ static real
__pick_shape_coeff(const real *__s, int j, int shift)
{
  real s;
  if (j == shift - 1) {
    s = __s[0];
  } else if (j == shift + 0) {
    s = __s[1];
  } else if (j == shift + 1) {
    s = real(1.) - __s[0] - __s[1];
  } else {
    s = real(0.);
  }
  return s;
}

// ----------------------------------------------------------------------
#elif CACHE_SHAPE_ARRAYS == 3

#if DIM == DIM_Z

#define DECLARE_SHAPE_INFO			\
  real __hz[2];					\
  real __s0z[2], __s1z[2];			\
  short int __shift0[1], __shift1[1]		\

#define SHAPE_INFO_PARAMS __hz,			\
    __s0z, __s1z,				\
    __shift0, __shift1
#define SHAPE_INFO_ARGS real *__hz,			\
    real *__s0z, real *__s1z,				\
    short int *__shift0, short int *__shift1

#define SI_SHIFT0Z __shift0[0]
#define SI_SHIFT1Z __shift1[0]
#define SI_SHIFT10Z (__shift1[0] - __shift0[0])

__device__ static void
cache_shape_arrays(SHAPE_INFO_ARGS, real *h0, real *h1,
		   short int shift0z, short int shift1z)
{
  __hz[0] = h0[2];
  __hz[1] = h1[2];
  SI_SHIFT0Z = shift0z;
  SI_SHIFT1Z = shift1z - shift0z;
}

#else

#define DECLARE_SHAPE_INFO			\
  real __hy[2], __hz[2];			\
  real __s0y[2], __s1y[2];			\
  real __s0z[2], __s1z[2];			\
  short int __shift0[2], __shift1[2]		\

#define SHAPE_INFO_PARAMS __hy, __hz,			\
    __s0y, __s1y, __s0z, __s1z,				\
    __shift0, __shift1
#define SHAPE_INFO_ARGS real *__hy, real *__hz,			\
    real *__s0y, real *__s1y, real *__s0z, real *__s1z,		\
    short int *__shift0, short int *__shift1

#define SI_SHIFT0Y __shift0[0]
#define SI_SHIFT1Y __shift1[0]
#define SI_SHIFT10Y (__shift1[0] - __shift0[0])
#define SI_SHIFT0Z __shift0[1]
#define SI_SHIFT1Z __shift1[1]
#define SI_SHIFT10Z (__shift1[1] - __shift0[1])

__device__ static void
cache_shape_arrays(SHAPE_INFO_ARGS, real *h0, real *h1,
		   short int shift0y, short int shift0z,
		   short int shift1y, short int shift1z)
{
  __hy[0] = h0[1];
  __hy[1] = h1[1];
  __hz[0] = h0[2];
  __hz[1] = h1[2];
  SI_SHIFT0Y = shift0y;
  SI_SHIFT1Y = shift1y - shift0y;
  SI_SHIFT0Z = shift0z;
  SI_SHIFT1Z = shift1z - shift0z;
}

#endif

#define pick_shape_coeff(t, comp, j, shift) ({				\
      const int __y __attribute__((unused)) = 1;			\
      const int __z __attribute__((unused)) = 2;			\
      __pick_shape_coeff(j, shift, __ ## comp, __h ## comp[t]);		\
  })


__device__ static real
__pick_shape_coeff(int j, int shift, int d, real h)
{
  real s;
  if (j == shift - 1) {
    s = find_shape_coeff_d_shift(-1, h, 0);
  } else if (j == shift + 0) {
    s = find_shape_coeff_d_shift( 0, h, 0);
  } else if (j == shift + 1) {
    s = find_shape_coeff_d_shift(+1, h, 0);
  } else {
    s = real(0.);
  }
  return s;
}

// ----------------------------------------------------------------------
#elif CACHE_SHAPE_ARRAYS == 4

#if DIM == DIM_Z

#define DECLARE_SHAPE_INFO			\
  real __s0z[2], __s1z[2];			\
  short int __shift0[1], __shift1[1]		\

#define SHAPE_INFO_PARAMS __s0z, __s1z, __shift0, __shift1
#define SHAPE_INFO_ARGS real *__s0z, real *__s1z, short int *__shift0, short int *__shift1

#define SI_SHIFT0Z __shift0[0]
#define SI_SHIFT1Z __shift1[0]
#define SI_SHIFT10Z (__shift1[0] - __shift0[0])

#elif DIM == DIM_YZ

__shared__ real ___s0y[2][THREADS_PER_BLOCK];
__shared__ real ___s1y[2][THREADS_PER_BLOCK];
__shared__ real ___s0z[2][THREADS_PER_BLOCK];
__shared__ real ___s1z[2][THREADS_PER_BLOCK];
__shared__ short int ___shift0[2][THREADS_PER_BLOCK];
__shared__ short int ___shift1[2][THREADS_PER_BLOCK];

#define DECLARE_SHAPE_INFO	do {} while (0)

#define SHAPE_INFO_PARAMS 0
#define SHAPE_INFO_ARGS int *unused

#define SI_SHIFT0Y ___shift0[0][threadIdx.x]
#define SI_SHIFT1Y ___shift1[0][threadIdx.x]
#define SI_SHIFT10Y (___shift1[0][threadIdx.x] - ___shift0[0][threadIdx.x])
#define SI_SHIFT0Z ___shift0[1][threadIdx.x]
#define SI_SHIFT1Z ___shift1[1][threadIdx.x]
#define SI_SHIFT10Z (___shift1[1][threadIdx.x] - ___shift0[1][threadIdx.x])

#endif

__device__ static inline void
calc_shape_coeff(real *__s, real h)
{
  __s[0] = find_shape_coeff_d_shift(-1, h, 0);
  __s[1] = find_shape_coeff_d_shift( 0, h, 0);
}

#if DIM == DIM_Z

__device__ static void
cache_shape_arrays(SHAPE_INFO_ARGS, real *h0, real *h1,
		   short int shift0z, short int shift1z)
{
  __shift0[0] = shift0z;
  __shift1[0] = shift1z;
  calc_shape_coeff(__s0z, h0[2]);
  calc_shape_coeff(__s1z, h1[2]);
}

#elif DIM == DIM_YZ

__device__ static void
cache_shape_arrays(SHAPE_INFO_ARGS, real *h0, real *h1,
		   short int shift0y, short int shift0z,
		   short int shift1y, short int shift1z)
{
  int tid = threadIdx.x;
  ___shift0[0][tid] = shift0y;
  ___shift1[0][tid] = shift1y;
  ___shift0[1][tid] = shift0z;
  ___shift1[1][tid] = shift1z;
  ___s0y[0][tid] = find_shape_coeff_d_shift(-1, h0[1], 0);
  ___s0y[1][tid] = find_shape_coeff_d_shift( 0, h0[1], 0);
  ___s1y[0][tid] = find_shape_coeff_d_shift(-1, h1[1], 0);
  ___s1y[1][tid] = find_shape_coeff_d_shift( 0, h1[1], 0);
  ___s0z[0][tid] = find_shape_coeff_d_shift(-1, h0[2], 0);
  ___s0z[1][tid] = find_shape_coeff_d_shift( 0, h0[2], 0);
  ___s1z[0][tid] = find_shape_coeff_d_shift(-1, h1[2], 0);
  ___s1z[1][tid] = find_shape_coeff_d_shift( 0, h1[2], 0);
}

#endif

#define pick_shape_coeff(t, comp, j, shift)	  \
  pick_shape_coeff_s ## t ## comp(j, shift)

#define MAKE_PICK_SHAPE_COEFF(s0y)				\
__device__ static real					        \
pick_shape_coeff_ ## s0y(int j, int shift)			\
{								\
  int tid = threadIdx.x;					\
  real s;							\
  if (j == shift - 1) {						\
    s = ___ ## s0y[0][tid];					\
  } else if (j == shift + 0) {					\
    s = ___ ## s0y[1][tid];					\
  } else if (j == shift + 1) {					\
    s = real(1.) - ___ ## s0y[0][tid] - ___ ## s0y[1][tid];	\
  } else {							\
    s = real(0.);						\
  }								\
  return s;							\
}

MAKE_PICK_SHAPE_COEFF(s0y)
MAKE_PICK_SHAPE_COEFF(s1y)
MAKE_PICK_SHAPE_COEFF(s0z)
MAKE_PICK_SHAPE_COEFF(s1z)

#else

#error unknown CACHE_SHAPE_ARRAYS

#endif

// FIXME use 2*SW
#define BLOCKSIZE ((BLOCKSIZE_X) * (BLOCKSIZE_Y + 6) * (BLOCKSIZE_Z + 6))
#define BLOCKSTRIDE (((BLOCKSIZE + 3) / 4) * 4)

#define scratch(m,jy,jz) (scratch[(m)*BLOCKSTRIDE +		\
				  ((jz)+3) * (BLOCKSIZE_Y+6) +	\
				  (jy)+3])

#define forall_j(x)	do {				\
    for (int j = -SW; j <= SW; j++) {			\
      x							\
     }						        \
} while (0)

#include "common_reduce.c"
#undef forall_j

#define forall_j(x)	do {				\
    for (int j = -2; j <= 1; j++) {			\
      x							\
     }						        \
} while (0)

#define reduce_sum_sdata reduce_sum_sdata4
#include "common_reduce.c"
#undef forall_j
#undef reduce_sum_sdata

// ======================================================================

__shared__ real sdata1[THREADS_PER_BLOCK];

__device__ static float
reduce_sum(float mySum)
{
  unsigned int tid = threadIdx.x;

  sdata1[tid] = mySum;
  __syncthreads();

  // do reduction in shared mem
  if (THREADS_PER_BLOCK >= 512) {
    if (tid < 256) {
      sdata1[tid] = mySum = mySum + sdata1[tid + 256];
    }
    __syncthreads();
  }
  if (THREADS_PER_BLOCK >= 256) {
    if (tid < 128) {
      sdata1[tid] = mySum = mySum + sdata1[tid + 128];
    }
    __syncthreads();
  }
  if (THREADS_PER_BLOCK >= 128) {
    if (tid < 64) {
      sdata1[tid] = mySum = mySum + sdata1[tid + 64];
    }
    __syncthreads();
  }

  if (tid < 32) {
    // now that we are using warp-synchronous programming (below)
    // we need to declare our shared memory volatile so that the compiler
    // doesn't reorder stores to it and induce incorrect behavior.
    volatile float* smem = sdata1;
    if (THREADS_PER_BLOCK >= 64) {
      smem[tid] = mySum = mySum + smem[tid + 32];
    }
    if (THREADS_PER_BLOCK >= 32) {
      smem[tid] = mySum = mySum + smem[tid + 16];
    }
    if (THREADS_PER_BLOCK >= 16) {
      smem[tid] = mySum = mySum + smem[tid + 8];
    }
    if (THREADS_PER_BLOCK >= 8) {
      smem[tid] = mySum = mySum + smem[tid + 4];
    }
    if (THREADS_PER_BLOCK >= 4) {
      smem[tid] = mySum = mySum + smem[tid + 2];
    }
    if (THREADS_PER_BLOCK >= 2) {
      smem[tid] = mySum = mySum + smem[tid + 1];
    }
  }

  return mySum;
}

__device__ static float
reduce_sum_warp(float mySum)
{
  unsigned int tid = threadIdx.x;

  // now that we are using warp-synchronous programming (below)
  // we need to declare our shared memory volatile so that the compiler
  // doesn't reorder stores to it and induce incorrect behavior.
  volatile float* smem = sdata1;
  smem[tid] = mySum;
  smem[tid] = mySum = mySum + smem[tid + 16];
  smem[tid] = mySum = mySum + smem[tid + 8];
  smem[tid] = mySum = mySum + smem[tid + 4];
  smem[tid] = mySum = mySum + smem[tid + 2];
  smem[tid] = mySum = mySum + smem[tid + 1];
  
  return mySum;
}
