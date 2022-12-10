
#pragma once

#include <kg/Vec3.h>
#include <dim.hxx>
#include <psc_bits.h>
#include <psc/gtensor.h>

namespace psc
{
namespace deposit
{
namespace norm
{

// ---------------------------------------------------------------------------
// DepositNorm1st
//
// 1st order deposition given cell index l and normalized
// in-cell position h in [0,1[

template <typename R, typename D>
class Deposit1st
{
public:
  using real_t = R;
  using real3_t = gt::sarray<R, 3>;
  using dim_t = D;

private:
  template <typename F>
  void operator()(F& flds, const gt::sarray<int, 3>& l, const real3_t& h,
                  real_t val, dim_yz tag)
  {
    flds(0, l[1] + 0, l[2] + 0) += val * (1.f - h[1]) * (1.f - h[2]);
    flds(0, l[1] + 1, l[2] + 0) += val * h[1] * (1.f - h[2]);
    flds(0, l[1] + 0, l[2] + 1) += val * (1.f - h[1]) * h[2];
    flds(0, l[1] + 1, l[2] + 1) += val * h[1] * h[2];
  }

  template <typename F>
  void operator()(F& flds, const gt::sarray<int, 3>& l, const real3_t& h,
                  real_t val, dim_xyz tag)
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

public:
  template <typename F>
  void operator()(F& flds, const gt::sarray<int, 3>& l, const real3_t& h,
                  real_t val)
  {
    (*this)(flds, l, h, val, dim_t{});
  }
};

// ----------------------------------------------------------------------------
// deposit directly, assuming grid spacing == 1, pass in patch-relative position
// x

template <typename T, typename D>
class Deposit1stNc
{
public:
  template <typename F>
  void operator()(F& flds, const gt::sarray<int, 3>& ib,
                  const gt::sarray<T, 3>& x, T val)

  {
    gt::sarray<int, 3> l;
    gt::sarray<T, 3> h;
    for (int d = 0; d < 3; d++) {
      l[d] = fint(x[d]);
      h[d] = x[d] - l[d];
      l[d] -= ib[d];
    }
    Deposit1st<T, D> deposit;
    deposit(flds, l, h, val);
  }
};

template <typename T, typename D>
class Deposit1stCc
{
public:
  template <typename F>
  void operator()(F& flds, const gt::sarray<int, 3>& ib,
                  const gt::sarray<T, 3>& x, T val)

  {
    gt::sarray<int, 3> l;
    gt::sarray<T, 3> h;
    for (int d = 0; d < 3; d++) {
      l[d] = fint(x[d] - .5f);
      h[d] = x[d] - .5f - l[d];
      l[d] -= ib[d];
    }
    Deposit1st<T, D> deposit;
    deposit(flds, l, h, val);
  }
};

template <typename D, typename F, typename T>
void nc(F& flds, const gt::sarray<int, 3>& ib, const gt::sarray<T, 3>& x, T val)
{
  Deposit1stNc<T, D> deposit;
  deposit(flds, ib, x, val);
}

template <typename D, typename F, typename T>
void cc(F& flds, const gt::sarray<int, 3>& ib, const gt::sarray<T, 3>& x, T val)
{
  Deposit1stCc<T, D> deposit;
  deposit(flds, ib, x, val);
}

} // namespace norm

namespace code // deposition in code units
{

// ----------------------------------------------------------------------------
// Deposit1st
//
// Wrapper around Deposit to n.c. / c.c. grid that translates from code units

template <typename R, typename D, template <typename, typename> class Deposit>
class Deposit1st
{
public:
  using real_t = R;
  using dim_t = D;
  using real3_t = gt::sarray<real_t, 3>;

  Deposit1st(const real3_t& dx, real_t fnqs)
    : dxi_{real_t(1.) / dx}, fnqs_{fnqs}
  {}

  template <typename P, typename F>
  void operator()(const P& prt, const F& flds, const gt::shape_type<3>& ib,
                  real_t val)
  {
    real3_t x = prt.x() * dxi_;
    real_t value = fnqs_ * val;

    Deposit<real_t, dim_t> deposit;
    deposit(flds, ib, x, value);
  }

  real3_t dxi_;
  real_t fnqs_;
};

// ----------------------------------------------------------------------------
// Deposit1stNc
//
// Deposition to NC grid in code units

template <typename R, typename D>
using Deposit1stNc = Deposit1st<R, D, psc::deposit::norm::Deposit1stNc>;

// ----------------------------------------------------------------------------
// Deposit1stCc
//
// Deposition to CC grid in code units

template <typename R, typename D>
using Deposit1stCc = Deposit1st<R, D, psc::deposit::norm::Deposit1stCc>;

} // namespace code
} // namespace deposit
} // namespace psc
