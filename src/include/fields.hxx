
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
  using value_type = typename fields_t::value_type;
  using shape_type = typename fields_t::shape_type;
  using dim = D;

  Fields3d(const fields_t& e, const Int3& ib) : e_(e), data_(e.data()), ib_(ib)
  {}

  shape_type shape() const { return e_.shape(); }
  int shape(int d) const { return e_.shape(d); }

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
    assert(m >= 0 && m < shape(3));
    assert(i >= 0 && i < shape(0));
    assert(j >= 0 && j < shape(1));
    assert(k >= 0 && k < shape(2));
#endif

    return ((m * shape(2) + k) * shape(1) + j) * shape(0) + i;
  }

private:
  fields_t e_;
  value_type* data_;
  Int3 ib_;
};

#endif
