
#pragma once

#include <math.h>

#include <psc/moment.hxx>

// ======================================================================
// rho

template <typename S, typename D>
struct Moment_rho_1st_nc : ItemMomentCRTP<Moment_rho_1st_nc<S, D>, S>
{
  using Base = ItemMomentCRTP<Moment_rho_1st_nc<S, D>, S>;
  using moment_type =
    psc::moment::moment_rho<psc::deposit::code::Deposit1stNc, D>;

  using Base::Base;
};
