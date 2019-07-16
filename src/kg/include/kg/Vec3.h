
#ifndef KG_VEC3_H
#define KG_VEC3_H

#include <algorithm>
#include <cassert>
#include <memory>

#include "cuda_compat.h"

// ======================================================================
// array
//
// Basically the same (or at least, a subset) of std::array, but with CUDA
// support

template <typename T, std::size_t N>
struct array
{
  using value_type = T;
  using size_t = std::size_t;

  T arr[N];

  __host__ __device__ T operator[](size_t i) const { return arr[i]; }

  __host__ __device__ T& operator[](size_t i) { return arr[i]; }

  __host__ __device__ const T* data() const { return arr; }

  __host__ __device__ T* data() { return arr; }

  operator T*() { return arr; } // FIXME, should be data()
};

template <typename T, std::size_t N>
bool operator==(const array<T, N>& x, const array<T, N>& y)
{
  return std::equal(x.arr, x.arr + N, y.arr);
}

template <typename T, std::size_t N>
bool operator!=(const array<T, N>& x, const array<T, N>& y)
{
  return !(x == y);
}

// ======================================================================
// Vec3

template <typename T>
struct Vec3
{
  static const int N = 3;
  using value_type = T;
  using size_t = std::size_t;

  T arr[N];

  __host__ __device__ T operator[](size_t i) const { return arr[i]; }

  __host__ __device__ T& operator[](size_t i) { return arr[i]; }

  __host__ __device__ const T* data() const { return arr; }

  __host__ __device__ T* data() { return arr; }

  // ----------------------------------------------------------------------
  // construct from pointer to values

  __host__ __device__ static Vec3 fromPointer(const T* p)
  {
    return {p[0], p[1], p[2]};
  }

  // ----------------------------------------------------------------------
  // converting to Vec3 of different type (e.g., float -> double)

  template <typename U>
  __host__ __device__ explicit operator Vec3<U>() const
  {
    return {U((*this)[0]), U((*this)[1]), U((*this)[2])};
  }

  // ----------------------------------------------------------------------
  // arithmetic

  __host__ __device__ Vec3 operator-() const
  {
    Vec3 res;
    for (int i = 0; i < 3; i++) {
      res[i] = -(*this)[i];
    }
    return res;
  }

  __host__ __device__ Vec3& operator+=(const Vec3& w)
  {
    for (int i = 0; i < 3; i++) {
      (*this)[i] += w[i];
    }
    return *this;
  }

  __host__ __device__ Vec3& operator-=(const Vec3& w)
  {
    for (int i = 0; i < 3; i++) {
      (*this)[i] -= w[i];
    }
    return *this;
  }

  __host__ __device__ Vec3& operator*=(const Vec3& w)
  {
    for (int i = 0; i < 3; i++) {
      (*this)[i] *= w[i];
    }
    return *this;
  }

  __host__ __device__ Vec3& operator*=(T s)
  {
    for (int i = 0; i < 3; i++) {
      (*this)[i] *= s;
    }
    return *this;
  }

  __host__ __device__ Vec3& operator/=(const Vec3& w)
  {
    for (int i = 0; i < 3; i++) {
      (*this)[i] /= w[i];
    }
    return *this;
  }

  // conversion to pointer

  __host__ __device__ operator const T*() const { return data(); }

  __host__ __device__ operator T*() { return data(); }
};

template <typename T>
bool operator==(const Vec3<T>& x, const Vec3<T>& y)
{
  static const int N = 3;
  return std::equal(x.arr, x.arr + N, y.arr);
}

template <typename T, std::size_t N>
bool operator!=(const Vec3<T>& x, const Vec3<T>& y)
{
  return !(x == y);
}

template <typename T>
__host__ __device__ Vec3<T> operator+(const Vec3<T>& v, const Vec3<T>& w)
{
  Vec3<T> res = v;
  res += w;
  return res;
}

template <typename T>
__host__ __device__ Vec3<T> operator-(const Vec3<T>& v, const Vec3<T>& w)
{
  Vec3<T> res = v;
  res -= w;
  return res;
}

template <typename T>
__host__ __device__ Vec3<T> operator*(const Vec3<T>& v, const Vec3<T>& w)
{
  Vec3<T> res = v;
  res *= w;
  return res;
}

template <typename T>
__host__ __device__ Vec3<T> operator*(T s, const Vec3<T>& v)
{
  Vec3<T> res = v;
  res *= s;
  return res;
}

template <typename T>
__host__ __device__ Vec3<T> operator/(const Vec3<T>& v, const Vec3<T>& w)
{
  Vec3<T> res = v;
  res /= w;
  return res;
}

using Int3 = Vec3<int>;
using UInt3 = Vec3<unsigned int>;

#endif
