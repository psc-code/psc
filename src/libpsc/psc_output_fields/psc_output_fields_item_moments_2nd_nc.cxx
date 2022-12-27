
#include "fields.hxx"
#include <psc/moment.hxx>

#include <math.h>

template <typename S, typename D>
using Moment_rho_2nd_nc =
  ItemMoment<psc::moment::moment_rho<psc::deposit::code::Deposit2ndNc, D>, S>;
