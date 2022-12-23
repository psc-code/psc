
#pragma once

#include <math.h>

#include <psc/moment.hxx>

// ======================================================================
// rho

template <typename MF, typename D>
struct Moment_rho_1st_nc : ItemMomentCRTP<Moment_rho_1st_nc<MF, D>, MF>
{
  using Base = ItemMomentCRTP<Moment_rho_1st_nc<MF, D>, MF>;
  using moment_type =
    psc::moment::moment_rho<psc::deposit::code::Deposit1stNc, D>;

  using Base::Base;
};
