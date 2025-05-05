
#pragma once

#include "centering.hxx"

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
// Deposit1st
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

// ---------------------------------------------------------------------------
// Deposit2nd
//
// 2nd order deposition given cell index l and normalized
// negated in-cell position h in [-1/2,1/2]

template <typename R, typename D>
class Deposit2nd
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
    real_t gmy = .5f * (.5f + h[1]) * (.5f + h[1]);
    real_t gmz = .5f * (.5f + h[2]) * (.5f + h[2]);
    real_t g0y = .75f - h[1] * h[1];
    real_t g0z = .75f - h[2] * h[2];
    real_t g1y = .5f * (.5f - h[1]) * (.5f - h[1]);
    real_t g1z = .5f * (.5f - h[2]) * (.5f - h[2]);

    // clang-format off
    flds(0, l[1] - 1, l[2] - 1) += gmy * gmz * val;
    flds(0, l[1]    , l[2] - 1) += g0y * gmz * val;
    flds(0, l[1] + 1, l[2] - 1) += g1y * gmz * val;
    flds(0, l[1] - 1, l[2]    ) += gmy * g0z * val;
    flds(0, l[1]    , l[2]    ) += g0y * g0z * val;
    flds(0, l[1] + 1, l[2]    ) += g1y * g0z * val;
    flds(0, l[1] - 1, l[2] + 1) += gmy * g1z * val;
    flds(0, l[1]    , l[2] + 1) += g0y * g1z * val;
    flds(0, l[1] + 1, l[2] + 1) += g1y * g1z * val;
    // clang-format on
  }

  template <typename F>
  void operator()(F& flds, const gt::sarray<int, 3>& l, const real3_t& h,
                  real_t val, dim_xyz tag)
  {
    real_t gmx = .5f * (.5f + h[0]) * (.5f + h[0]);
    real_t gmy = .5f * (.5f + h[1]) * (.5f + h[1]);
    real_t gmz = .5f * (.5f + h[2]) * (.5f + h[2]);
    real_t g0x = .75f - h[0] * h[0];
    real_t g0y = .75f - h[1] * h[1];
    real_t g0z = .75f - h[2] * h[2];
    real_t g1x = .5f * (.5f - h[0]) * (.5f - h[0]);
    real_t g1y = .5f * (.5f - h[1]) * (.5f - h[1]);
    real_t g1z = .5f * (.5f - h[2]) * (.5f - h[2]);

    // clang-format off
    flds(l[0] - 1, l[1] - 1, l[2] - 1) += gmx * gmy * gmz * val;
    flds(l[0]    , l[1] - 1, l[2] - 1) += g0x * gmy * gmz * val;
    flds(l[0] + 1, l[1] - 1, l[2] - 1) += g1x * gmy * gmz * val;
    flds(l[0] - 1, l[1]    , l[2] - 1) += gmx * g0y * gmz * val;
    flds(l[0]    , l[1]    , l[2] - 1) += g0x * g0y * gmz * val;
    flds(l[0] + 1, l[1]    , l[2] - 1) += g1x * g0y * gmz * val;
    flds(l[0] - 1, l[1] + 1, l[2] - 1) += gmx * g1y * gmz * val;
    flds(l[0]    , l[1] + 1, l[2] - 1) += g0x * g1y * gmz * val;
    flds(l[0] + 1, l[1] + 1, l[2] - 1) += g1x * g1y * gmz * val;
    flds(l[0] - 1, l[1] - 1, l[2]    ) += gmx * gmy * g0z * val;
    flds(l[0]    , l[1] - 1, l[2]    ) += g0x * gmy * g0z * val;
    flds(l[0] + 1, l[1] - 1, l[2]    ) += g1x * gmy * g0z * val;
    flds(l[0] - 1, l[1]    , l[2]    ) += gmx * g0y * g0z * val;
    flds(l[0]    , l[1]    , l[2]    ) += g0x * g0y * g0z * val;
    flds(l[0] + 1, l[1]    , l[2]    ) += g1x * g0y * g0z * val;
    flds(l[0] - 1, l[1] + 1, l[2]    ) += gmx * g1y * g0z * val;
    flds(l[0]    , l[1] + 1, l[2]    ) += g0x * g1y * g0z * val;
    flds(l[0] + 1, l[1] + 1, l[2]    ) += g1x * g1y * g0z * val;
    flds(l[0] - 1, l[1] - 1, l[2] + 1) += gmx * gmy * g1z * val;
    flds(l[0]    , l[1] - 1, l[2] + 1) += g0x * gmy * g1z * val;
    flds(l[0] + 1, l[1] - 1, l[2] + 1) += g1x * gmy * g1z * val;
    flds(l[0] - 1, l[1]    , l[2] + 1) += gmx * g0y * g1z * val;
    flds(l[0]    , l[1]    , l[2] + 1) += g0x * g0y * g1z * val;
    flds(l[0] + 1, l[1]    , l[2] + 1) += g1x * g0y * g1z * val;
    flds(l[0] - 1, l[1] + 1, l[2] + 1) += gmx * g1y * g1z * val;
    flds(l[0]    , l[1] + 1, l[2] + 1) += g0x * g1y * g1z * val;
    flds(l[0] + 1, l[1] + 1, l[2] + 1) += g1x * g1y * g1z * val;
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

template <centering::Centering C>
class DepositCentering
{
public:
  static const centering::Centering CENTERING = C;
};

// ----------------------------------------------------------------------------
// deposit directly, assuming grid spacing == 1, pass in patch-relative position
// x

template <typename T, typename D>
class Deposit1stNc : public DepositCentering<centering::NC>
{
public:
  static std::string suffix() { return "_1st_nc"; }

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
class Deposit1stCc : public DepositCentering<centering::CC>
{
public:
  static std::string suffix() { return "_1st_cc"; }

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

// ----------------------------------------------------------------------------
// deposit directly, assuming grid spacing == 1, pass in patch-relative position
// x

template <typename T, typename D>
class Deposit2ndNc : public DepositCentering<centering::NC>
{
public:
  static std::string suffix() { return "_2nd_nc"; }

  template <typename F>
  void operator()(F& flds, const gt::sarray<int, 3>& ib,
                  const gt::sarray<T, 3>& x, T val)

  {
    gt::sarray<int, 3> l;
    gt::sarray<T, 3> h;
    for (int d = 0; d < 3; d++) {
      l[d] = nint(x[d]);
      h[d] = l[d] - x[d]; // negated!
      l[d] -= ib[d];
    }
    Deposit2nd<T, D> deposit;
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
// Deposit
//
// Wrapper around DepositNorm to n.c. / c.c. grid that translates from code
// units

template <typename R, typename D,
          template <typename, typename> class DEPOSIT_NORM>
class Deposit
{
public:
  using real_t = R;
  using dim_t = D;
  using real3_t = gt::sarray<real_t, 3>;
  using DepositNorm = DEPOSIT_NORM<real_t, dim_t>;
  static const centering::Centering CENTERING = DepositNorm::CENTERING;

  static std::string suffix() { return DepositNorm::suffix(); }

  Deposit(const real3_t& dx, real_t fnqs) : dxi_{real_t(1.) / dx}, fnqs_{fnqs}
  {}

  template <typename F>
  void operator()(const F& flds, const gt::shape_type<3>& ib, const real3_t& xi,
                  real_t val)
  {
    real3_t x = xi * dxi_;
    real_t value = fnqs_ * val;

    DepositNorm deposit;
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
using Deposit1stNc = Deposit<R, D, psc::deposit::norm::Deposit1stNc>;

// ----------------------------------------------------------------------------
// Deposit1stCc
//
// Deposition to CC grid in code units

template <typename R, typename D>
using Deposit1stCc = Deposit<R, D, psc::deposit::norm::Deposit1stCc>;

// ----------------------------------------------------------------------------
// Deposit2ndNc
//
// Deposition to NC grid in code units

template <typename R, typename D>
using Deposit2ndNc = Deposit<R, D, psc::deposit::norm::Deposit2ndNc>;

} // namespace code
} // namespace deposit
} // namespace psc
