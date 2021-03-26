
#ifndef FIELDS_HXX
#define FIELDS_HXX

#include "dim.hxx"
#include "kg/Vec3.h"

#include <gtensor/gtensor.h>

template <typename F, typename D = dim_xyz>
class Fields3d
{
public:
  using fields_t = F;
  using real_t = typename fields_t::real_t;
  using dim = D;

  Fields3d(const fields_t& f)
    : data_(
        const_cast<typename fields_t::real_t*>(f.storage().data())), // FIXME
      n_comp_(f.n_comps()),
      ib(f.ib()),
      im(f.im())
  {}

  GT_INLINE const real_t operator()(int m, int i, int j, int k) const
  {
    return data_[index(m, i, j, k)];
  }

  GT_INLINE real_t& operator()(int m, int i, int j, int k)
  {
    return data_[index(m, i, j, k)];
  }

private:
  GT_INLINE int index(int m, int i_, int j_, int k_) const
  {
    int i = dim::InvarX::value ? 0 : i_;
    int j = dim::InvarY::value ? 0 : j_;
    int k = dim::InvarZ::value ? 0 : k_;

#if defined(BOUNDS_CHECK) && !defined(__CUDACC__)
    assert(m >= 0 && m < n_comp_);
    assert(i >= ib[0] && i < ib[0] + im[0]);
    assert(j >= ib[1] && j < ib[1] + im[1]);
    assert(k >= ib[2] && k < ib[2] + im[2]);
#endif

    return ((((((m)) * im[2] + (k - ib[2])) * im[1] + (j - ib[1])) * im[0] +
             (i - ib[0])));
  }

private:
  real_t* data_;
  Int3 ib, im;
  int n_comp_;
};

#endif
