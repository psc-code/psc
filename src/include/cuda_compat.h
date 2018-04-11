
#ifndef CUDA_COMPAT_H
#define CUDA_COMPAT_H

#include <cmath>


// ----------------------------------------------------------------------
// cuda compatibility stuff

#ifndef __CUDACC__

#define __host__
#define __device__

typedef struct { float x; float y; float z; float w; } float4;

#endif

template<typename T>
T rsqrt(T x);


#ifndef __CUDA_ARCH__

// on the host, provide default implementation as 1/sqrt()

template<typename T>
T rsqrt(T x)
{
  return T(1.) / std::sqrt(x);
}

#endif

#ifdef __CUDA_ARCH__

// on the device, provide only specialized implementation for float,
// using rsqrtf() intrinsic

template<>
float rsqrt(float x)
{
  return rsqrtf(x);
}

#endif


#endif

