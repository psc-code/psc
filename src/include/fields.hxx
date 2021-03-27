
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
  using value_type = typename fields_t::real_t;
  using dim = D;

  Fields3d(const fields_t& f, const Int3& ib)
    : data_(const_cast<value_type*>(f.storage().data())), // FIXME
      n_comp_(f.storage().shape(3)),
      ib(ib),
      im({f.storage().shape(0), f.storage().shape(1), f.storage().shape(2)})
  {}

  const value_type operator()(int m, int i, int j, int k) const
  {
    return data_[index(m, i, j, k)];
  }

  GT_INLINE value_type& operator()(int m, int i, int j, int k)
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
  value_type* data_;
  Int3 ib, im;
  int n_comp_;
};

template <typename F, typename D = dim_xyz>
class _Fields3d
{
public:
  using fields_t = F;
  using value_type = typename fields_t::value_type;
  using dim = D;

  _Fields3d(const fields_t& f, const Int3& ib)
    : data_(f.data()), shape_(f.shape()), ib_(ib)
  {}

  const value_type& operator()(int m, int i, int j, int k) const
  {
    return data_[index(m, i, j, k)];
  }

  value_type& operator()(int m, int i, int j, int k)
  {
    return data_[index(m, i, j, k)];
  }

private:
  int index(int m, int i_, int j_, int k_) const
  {
    int i = dim::InvarX::value ? 0 : i_ - ib_[0];
    int j = dim::InvarY::value ? 0 : j_ - ib_[1];
    int k = dim::InvarZ::value ? 0 : k_ - ib_[2];

#ifdef BOUNDS_CHECK
    assert(m >= 0 && m < shape_[3]);
    assert(i >= 0 && i < shape_[0]);
    assert(j >= 0 && j < shape_[1]);
    assert(k >= 0 && k < shape_[2]);
#endif

    return ((m * shape_[2] + k) * shape_[1] + j) * shape_[0] + i;
  }

private:
  value_type* data_;
  gt::shape_type<4> shape_;
  Int3 ib_;
};

#endif
