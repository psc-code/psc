
#pragma once

#include <math.h>

#include <psc/moment.hxx>

// ======================================================================
// rho

template <typename S, typename D>
using Moment_rho_1st_nc =
  ItemMoment<psc::moment::moment_rho<psc::deposit::code::Deposit1stNc, D>, S>;
