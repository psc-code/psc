
#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>

#include "cuda_compat.h"
#include <kg/Macros.h>
#include <psc/gtensor.h>

namespace kg
{

// ======================================================================
// Vec

template <typename T, std::size_t N>
struct Vec : gt::sarray<T, N>
{
  using base_type = gt::sarray<T, N>;
  using value_type = T;
  using size_t = std::size_t;

  using base_type::base_type;

  // ----------------------------------------------------------------------
  // construct from pointer to values

  KG_INLINE static Vec fromPointer(const T* p)
  {
    Vec res;
    for (auto& val : res) {
      val = *p++;
    }
    return res;
  }

  // ----------------------------------------------------------------------
  // converting to Vec of different type (e.g., float -> double)

  template <typename U>
  KG_INLINE explicit operator Vec<U, N>() const
  {
    Vec<U, N> res;
    for (size_t i = 0; i < N; i++) {
      res[i] = U((*this)[i]);
    }
    return res;
  }

  // ----------------------------------------------------------------------
  // arithmetic

  KG_INLINE Vec operator-() const
  {
    Vec res;
    for (size_t i = 0; i < N; i++) {
      res[i] = -(*this)[i];
    }
    return res;
  }

  KG_INLINE Vec& operator+=(const Vec& w)
  {
    for (size_t i = 0; i < N; i++) {
      (*this)[i] += w[i];
    }
    return *this;
  }

  KG_INLINE Vec& operator-=(const Vec& w)
  {
    for (size_t i = 0; i < N; i++) {
      (*this)[i] -= w[i];
    }
    return *this;
  }

  KG_INLINE Vec& operator*=(const Vec& w)
  {
    for (size_t i = 0; i < N; i++) {
      (*this)[i] *= w[i];
    }
    return *this;
  }

  KG_INLINE Vec& operator*=(T s)
  {
    for (size_t i = 0; i < N; i++) {
      (*this)[i] *= s;
    }
    return *this;
  }

  KG_INLINE Vec& operator/=(const Vec& w)
  {
    for (size_t i = 0; i < N; i++) {
      (*this)[i] /= w[i];
    }
    return *this;
  }

  // conversion to pointer

  KG_INLINE operator const T*() const { return this->data(); }

  KG_INLINE operator T*() { return this->data(); }
};

template <typename T, std::size_t N>
bool operator==(const Vec<T, N>& x, const Vec<T, N>& y)
{
  return std::equal(x.begin(), x.end(), y.begin());
}

template <typename T, std::size_t N>
bool operator!=(const Vec<T, N>& x, const Vec<T, N>& y)
{
  return !(x == y);
}

template <typename T, std::size_t N>
KG_INLINE Vec<T, N> operator+(const Vec<T, N>& v, const Vec<T, N>& w)
{
  Vec<T, N> res = v;
  res += w;
  return res;
}

template <typename T, std::size_t N>
KG_INLINE Vec<T, N> operator-(const Vec<T, N>& v, const Vec<T, N>& w)
{
  Vec<T, N> res = v;
  res -= w;
  return res;
}

template <typename T, std::size_t N>
KG_INLINE Vec<T, N> operator*(const Vec<T, N>& v, const Vec<T, N>& w)
{
  Vec<T, N> res = v;
  res *= w;
  return res;
}

template <typename T, std::size_t N>
KG_INLINE Vec<T, N> operator*(T s, const Vec<T, N>& v)
{
  Vec<T, N> res = v;
  res *= s;
  return res;
}

template <typename T, std::size_t N>
KG_INLINE Vec<T, N> operator/(const Vec<T, N>& v, const Vec<T, N>& w)
{
  Vec<T, N> res = v;
  res /= w;
  return res;
}

template <typename T, std::size_t N>
KG_INLINE std::ostream& operator<<(std::ostream& os, const Vec<T, N>& v)
{
  os << "Vec<T," << N << ">{";
  for (std::size_t n = 0; n < N; n++) {
    os << v[n];
    if (n < N - 1) {
      os << ", ";
    }
  }
  os << "}";
  return os;
}

} // namespace kg

template <typename T>
using Vec3 = kg::Vec<T, 3>;

using Int3 = Vec3<int>;
using UInt3 = Vec3<unsigned int>;
using Float3 = Vec3<float>;
using Double3 = Vec3<double>;
