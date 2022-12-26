
#include "fields.hxx"
#include <psc/moment.hxx>

#include <math.h>

// ======================================================================
// rho

template <typename MF, typename D>
class Moment_rho_2nd_nc
  : public ItemMomentCRTP<Moment_rho_2nd_nc<MF, D>, typename MF::Storage>
{
public:
  using Base = ItemMomentCRTP<Moment_rho_2nd_nc<MF, D>, typename MF::Storage>;
  using moment_type =
    psc::moment::moment_rho<psc::deposit::code::Deposit2ndNc, D>;

  using Base::Base;
};
