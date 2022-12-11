
#pragma once

#ifdef PSC_HAVE_RMM
#include <rmm/mr/device/thrust_allocator_adaptor.hpp>

#define GTENSOR_DEFAULT_DEVICE_ALLOCATOR(T) rmm::mr::thrust_allocator<T>
#endif

#include <gtensor/gtensor.h>

namespace
{

template <typename T, gt::size_type N>
GT_INLINE gt::sarray<T, N> operator*(T a, const gt::sarray<T, N>& x)
{
  gt::sarray<T, N> res;
  for (gt::size_type i = 0; i < N; i++) {
    res[i] = a * x[i];
  }
  return res;
}

template <typename T, gt::size_type N>
GT_INLINE gt::sarray<T, N> operator/(const gt::sarray<T, N>& x, T a)
{
  gt::sarray<T, N> res;
  for (gt::size_type i = 0; i < N; i++) {
    res[i] = x[i] / a;
  }
  return res;
}

template <typename T, gt::size_type N>
GT_INLINE gt::sarray<T, N> operator/(T a, const gt::sarray<T, N>& x)
{
  gt::sarray<T, N> res;
  for (gt::size_type i = 0; i < N; i++) {
    res[i] = a / x[i];
  }
  return res;
}

template <typename T, gt::size_type N>
GT_INLINE gt::sarray<T, N> operator-(const gt::sarray<T, N>& x)
{
  gt::sarray<T, N> res;
  for (gt::size_type i = 0; i < N; i++) {
    res[i] = -x[i];
  }
  return res;
}

template <typename T, gt::size_type N>
GT_INLINE gt::sarray<T, N> operator+(const gt::sarray<T, N>& x,
                                     const gt::sarray<T, N>& y)
{
  gt::sarray<T, N> res;
  for (gt::size_type i = 0; i < N; i++) {
    res[i] = x[i] + y[i];
  }
  return res;
}

template <typename T, gt::size_type N>
GT_INLINE gt::sarray<T, N> operator-(const gt::sarray<T, N>& x,
                                     const gt::sarray<T, N>& y)
{
  gt::sarray<T, N> res;
  for (gt::size_type i = 0; i < N; i++) {
    res[i] = x[i] - y[i];
  }
  return res;
}

template <typename T, gt::size_type N>
GT_INLINE gt::sarray<T, N> operator*(const gt::sarray<T, N>& x,
                                     const gt::sarray<T, N>& y)
{
  gt::sarray<T, N> res;
  for (gt::size_type i = 0; i < N; i++) {
    res[i] = x[i] * y[i];
  }
  return res;
}

template <typename T, gt::size_type N>
GT_INLINE gt::sarray<T, N> operator/(const gt::sarray<T, N>& x,
                                     const gt::sarray<T, N>& y)
{
  gt::sarray<T, N> res;
  for (gt::size_type i = 0; i < N; i++) {
    res[i] = x[i] / y[i];
  }
  return res;
}

} // namespace
