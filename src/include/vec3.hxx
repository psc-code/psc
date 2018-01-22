
#ifndef VEC3_HXX
#define VEC3_HXX

#include <array>
#include <memory>
#include <cassert>

// ======================================================================
// Vec3

template<typename T>
struct Vec3 : std::array<T, 3>
{
  using Base = std::array<T, 3>;

  using Base::data;
  
  Vec3() = default;

  Vec3(const Vec3&) = default;

  Vec3(const T *p)
  {
    for (int i = 0; i < 3; i++) {
      new (&(*this)[i])	T(p[i]); // placement new -- not really necessary
    }
  }

  Vec3(std::initializer_list<T> l)
  {
    assert(l.size() == 3);
    std::uninitialized_copy(l.begin(), l.end(), data());
  }

  // convert from different type (e.g., float -> double)
  template<typename U>
  Vec3(const Vec3<U>& u)
  {
    for (int i = 0; i < 3; i++) {
      new (&(*this)[i])	T(u[i]); // placement new -- not really necessary
    }
  }

  // arithmetic

  Vec3& operator/=(const Vec3& w) {
    for (int i = 0; i < 3; i++) {
      (*this)[i] /= w[i];
    }
    return *this;
  }
  
  Vec3 operator/(const Vec3& w) const {
    Vec3 res;
    for (int i = 0; i < 3; i++) {
      res[i] = (*this)[i] / w[i];
    }
    return res;
  }
  
  // conversion to pointer
  
  operator const T* () const { return data(); }
  operator T* ()             { return data(); }

  bool operator==(const Vec3& other) const
  {
    return *static_cast<const Base*>(this) == other;
  }
};

using Int3 = Vec3<int>;

#endif

