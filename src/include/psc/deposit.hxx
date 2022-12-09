
#pragma once

namespace psc
{

// ---------------------------------------------------------------------------
// DepositNc

template <typename R, typename D>
class DepositNc
{
public:
  using real_t = R;
  using real3_t = Vec3<R>;
  using dim_t = D;

  template <typename F>
  void operator()(F& flds, real3_t x)
  {
    Int3 l;
    real3_t h;
    for (int d = 0; d < 3; d++) {
      l[d] = fint(x[d]);
      h[d] = x[d] - l[d];
    }
    (*this)(flds, l, h);
  }

  template <typename F>
  void operator()(F& flds, Int3 l, real3_t h, dim_yz tag)
  {
    flds(0, l[1] + 0, l[2] + 0) += (1.f - h[1]) * (1.f - h[2]);
    flds(0, l[1] + 1, l[2] + 0) += h[1] * (1.f - h[2]);
    flds(0, l[1] + 0, l[2] + 1) += (1.f - h[1]) * h[2];
    flds(0, l[1] + 1, l[2] + 1) += h[1] * h[2];
  }

  template <typename F>
  void operator()(F& flds, Int3 l, real3_t h, dim_xyz tag)
  {
    flds(l[0] + 0, l[1] + 0, l[2] + 0) +=
      (1.f - h[0]) * (1.f - h[1]) * (1.f - h[2]);
    flds(l[0] + 1, l[1] + 0, l[2] + 0) += h[0] * (1.f - h[1]) * (1.f - h[2]);
    flds(l[0] + 0, l[1] + 1, l[2] + 0) += (1.f - h[0]) * h[1] * (1.f - h[2]);
    flds(l[0] + 1, l[1] + 1, l[2] + 0) += h[0] * h[1] * (1.f - h[2]);
    flds(l[0] + 0, l[1] + 0, l[2] + 1) += (1.f - h[0]) * (1.f - h[1]) * h[2];
    flds(l[0] + 1, l[1] + 0, l[2] + 1) += h[0] * (1.f - h[1]) * h[2];
    flds(l[0] + 0, l[1] + 1, l[2] + 1) += (1.f - h[0]) * h[1] * h[2];
    flds(l[0] + 1, l[1] + 1, l[2] + 1) += h[0] * h[1] * h[2];
  }

  template <typename F>
  void operator()(F& flds, Int3 l, real3_t h)
  {
    (*this)(flds, l, h, dim_t{});
  }
};

} // namespace psc
