
#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>

#include "cuda_compat.h"
#include <kg/Macros.h>
#include <psc/gtensor.h>
#include "psc_bits.h"

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

  KG_INLINE Vec& operator+=(T s)
  {
    for (size_t i = 0; i < N; i++) {
      (*this)[i] += s;
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

  KG_INLINE Vec& operator-=(T s)
  {
    for (size_t i = 0; i < N; i++) {
      (*this)[i] -= s;
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

  KG_INLINE Vec& operator/=(T s)
  {
    for (size_t i = 0; i < N; i++) {
      (*this)[i] /= s;
    }
    return *this;
  }

  /**
   * @brief Calculates the sum of this vec's elements.
   * @return the sum
   */
  KG_INLINE T sum() const
  {
    T sum{0};
    for (size_t i = 0; i < N; i++) {
      sum += (*this)[i];
    }
    return sum;
  }

  /**
   * @brief Calculates the product of this vec's elements.
   * @return the product
   */
  KG_INLINE T prod() const
  {
    T prod{1};
    for (size_t i = 0; i < N; i++) {
      prod *= (*this)[i];
    }
    return prod;
  }

  /**
   * @brief Calculates the dot product of this and another vec.
   * @param w the other vec
   * @return the dot product
   */
  KG_INLINE T dot(const Vec& w) const
  {
    T sum{0};
    for (size_t i = 0; i < N; i++) {
      sum += (*this)[i] * w[i];
    }
    return sum;
  }

  /**
   * @brief Calculates the cross product of this and another vec. This method
   * assumes `N` is 3.
   * @param w the right-hand-side vec
   * @return the cross product
   */
  KG_INLINE Vec cross(const Vec& w) const
  {
    return {(*this)[1] * w[2] - (*this)[2] * w[1],
            (*this)[2] * w[0] - (*this)[0] * w[2],
            (*this)[0] * w[1] - (*this)[1] * w[0]};
  }

  /**
   * @brief Calculates the square of the magnitude of this vec.
   * @return the squared magnitude
   */
  KG_INLINE T mag2() const { return this->dot(*this); }

  /**
   * @brief Calculates the magnitude of this vec by taking the square root of
   * the squared magnitude.
   * @return the magnitude
   */
  KG_INLINE T mag() const { return sqrt(this->mag2()); }

  /**
   * @brief Determines the maximal value of any of this vec's elements.
   * @return the maximal value
   */
  KG_INLINE T max() const
  {
    return *std::max_element(this->begin(), this->end());
  }

  /**
   * @brief Determines the minimal value of any of this vec's elements.
   * @return the minimal value
   */
  KG_INLINE T min() const
  {
    return *std::min_element(this->begin(), this->end());
  }

  /**
   * @brief Calculates the elementwise maximum between two vecs.
   * @param v the first vec
   * @param w the second vec
   * @return a vec containing the maximal values
   */
  static KG_INLINE Vec max(const Vec& v, const Vec& w)
  {
    Vec res;
    for (int i = 0; i < N; i++) {
      res[i] = std::max(v[i], w[i]);
    }
    return res;
  }

  /**
   * @brief Calculates the elementwise minimum between two vecs.
   * @param v the first vec
   * @param w the second vec
   * @return a vec containing the minimal values
   */
  static KG_INLINE Vec min(const Vec& v, const Vec& w)
  {
    Vec res;
    for (int i = 0; i < N; i++) {
      res[i] = std::min(v[i], w[i]);
    }
    return res;
  }

  /**
   * @brief Calculates the (multiplicative) inverse of each of this vec's
   * elements.
   * @return a vec of the inverses
   */
  KG_INLINE Vec inv() const { return T(1) / *this; }

  /**
   * @brief Calculates the floor of each of this vec's elements.
   * @return an integer vec of the floors
   */
  KG_INLINE Vec<int, N> fint() const
  {
    Vec<int, N> res;
    for (size_t i = 0; i < N; i++) {
      res[i] = ::fint((*this)[i]);
    }
    return res;
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
KG_INLINE Vec<T, N> operator+(T s, const Vec<T, N>& v)
{
  Vec<T, N> res = v;
  res += s;
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
KG_INLINE Vec<T, N> operator-(T s, const Vec<T, N>& v)
{
  Vec<T, N> res;
  for (size_t i = 0; i < N; i++) {
    res[i] = s - v[i];
  }
  return res;
}
template <typename T, std::size_t N>
KG_INLINE Vec<T, N> operator-(const Vec<T, N>& v, T s)
{
  Vec<T, N> res = v;
  res -= s;
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
KG_INLINE Vec<T, N> operator/(T s, const Vec<T, N>& v)
{
  Vec<T, N> res;
  for (size_t i = 0; i < N; i++) {
    res[i] = s / v[i];
  }
  return res;
}

template <typename T, std::size_t N>
KG_INLINE Vec<T, N> operator/(const Vec<T, N>& v, T s)
{
  Vec<T, N> res = v;
  res /= s;
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
