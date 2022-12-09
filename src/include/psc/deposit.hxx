
#pragma once

#include <kg/Vec3.h>
#include <dim.hxx>
#include <psc_bits.h>
#include <psc/gtensor.h>

namespace psc
{

// ---------------------------------------------------------------------------
// DepositNc

template <typename R, typename D>
class DepositNc
{
public:
  using real_t = R;
  using real3_t = gt::sarray<R, 3>;
  using dim_t = D;

  template <typename F>
  void operator()(F& flds, const gt::shape_type<3>& ib, real3_t x, real_t val)
  {
    Int3 l;
    real3_t h;
    for (int d = 0; d < 3; d++) {
      l[d] = fint(x[d]);
      h[d] = x[d] - l[d];
      l[d] -= ib[d];
    }
    (*this)(flds, l, h, val);
  }

  template <typename F>
  void operator()(F& flds, Int3 l, real3_t h, real_t val, dim_yz tag)
  {
    flds(0, l[1] + 0, l[2] + 0) += val * (1.f - h[1]) * (1.f - h[2]);
    flds(0, l[1] + 1, l[2] + 0) += val * h[1] * (1.f - h[2]);
    flds(0, l[1] + 0, l[2] + 1) += val * (1.f - h[1]) * h[2];
    flds(0, l[1] + 1, l[2] + 1) += val * h[1] * h[2];
  }

  template <typename F>
  void operator()(F& flds, Int3 l, real3_t h, real_t val, dim_xyz tag)
  {
    // clang-format off
    flds(l[0] + 0, l[1] + 0, l[2] + 0) += val * (1.f - h[0]) * (1.f - h[1]) * (1.f - h[2]);
    flds(l[0] + 1, l[1] + 0, l[2] + 0) += val *        h[0]  * (1.f - h[1]) * (1.f - h[2]);
    flds(l[0] + 0, l[1] + 1, l[2] + 0) += val * (1.f - h[0]) *        h[1]  * (1.f - h[2]);
    flds(l[0] + 1, l[1] + 1, l[2] + 0) += val *        h[0]  *        h[1]  * (1.f - h[2]);
    flds(l[0] + 0, l[1] + 0, l[2] + 1) += val * (1.f - h[0]) * (1.f - h[1]) *        h[2];
    flds(l[0] + 1, l[1] + 0, l[2] + 1) += val *        h[0]  * (1.f - h[1]) *        h[2];
    flds(l[0] + 0, l[1] + 1, l[2] + 1) += val * (1.f - h[0]) *        h[1]  *        h[2];
    flds(l[0] + 1, l[1] + 1, l[2] + 1) += val *        h[0]  *        h[1]  *        h[2];
    // clang-format on
  }

  template <typename F>
  void operator()(F& flds, Int3 l, real3_t h, real_t val)
  {
    (*this)(flds, l, h, val, dim_t{});
  }
};

} // namespace psc
