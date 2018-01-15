
#ifndef VEC3_HXX
#define VEC3_HXX

#include <array>
#include <cassert>

// ======================================================================
// Vec3

template<typename T>
struct Vec3 : std::array<T, 3>
{
  using Base = std::array<T, 3>;

  using Base::data;
  
  Vec3() = default;
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
  
  operator const T* () const { return data(); }
  operator T* ()             { return data(); }

  bool operator==(const Vec3& other) const
  {
    return *static_cast<const Base*>(this) == other;
  }
};

using Int3 = Vec3<int>;

#endif

