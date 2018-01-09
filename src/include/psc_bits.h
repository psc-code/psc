
#ifndef PSC_BITS_H
#define PSC_BITS_H

#include <cmath>

#ifdef __CUDACC__
#endif

template<typename T>
int fint(T val)
{
#ifdef __CUDA_ARCH__
  return __float2int_rd(val);
#else
  return (int) std::floor(val);
#endif
}


#endif

