
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
  using value_type = typename Base::value_type;

  using Base::data;
  
  Vec3() = default;

  // ----------------------------------------------------------------------
  // copy ctor

  Vec3(const Vec3&) = default;

  // ----------------------------------------------------------------------
  // construct by broadcasting single value

  Vec3(T val)
  {
    for (int i = 0; i < 3; i++) {
      new (&(*this)[i])	T(val); // placement new -- not really necessary
    }    
  }

  // ----------------------------------------------------------------------
  // construct from pointer to values
  
  Vec3(const T *p)
  {
    for (int i = 0; i < 3; i++) {
      new (&(*this)[i])	T(p[i]); // placement new -- not really necessary
    }
  }

  // ----------------------------------------------------------------------
  // construct from initializer list
  
  Vec3(std::initializer_list<T> l)
  {
    assert(l.size() == 3);
    std::uninitialized_copy(l.begin(), l.end(), data());
  }

  // ----------------------------------------------------------------------
  // construct by converting from different type (e.g., float -> double)

  template<typename U>
  explicit Vec3(const Vec3<U>& u)
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
    Vec3 res = *this;
    res /= w;
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

