
#ifndef PSC_BITS_H
#define PSC_BITS_H

#include <cmath>

template<typename T>
int fint(T val)
{
#ifdef __CUDA_ARCH__
  return __float2int_rd(val);
#else
  return (int) std::floor(val);
#endif
}

template<typename T>
int nint(T val)
{
  return (int) std::round(val);
}

#define sqr(a) ((a) * (a))

using uint = unsigned int;


#endif

