
#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "fields_item_dive_cuda.hxx"

#define BND 2
#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

// ======================================================================

void cuda_mfields_calc_dive_yz(MfieldsStateCuda& mflds, MfieldsCuda& mdive)
{
  auto dx = mflds.grid().domain.dx;

  auto bnd = mflds.ibn();
  auto flds = mflds.gt().view(_s(bnd[0], -bnd[0]), _s(-1 + bnd[1], -bnd[1]),
                              _s(-1 + bnd[2], -bnd[2]));

  auto bd = mdive.ibn();
  auto dive =
    mdive.gt().view(_s(bd[0], -bd[0]), _s(bd[1], -bd[1]), _s(bd[2], -bd[2]), 0);

  auto ey = flds.view(_all, _s(1, _), _s(1, _), EY);
  auto eym = flds.view(_all, _s(_, -1), _s(1, _), EY);
  auto ez = flds.view(_all, _s(1, _), _s(1, _), EZ);
  auto ezm = flds.view(_all, _s(1, _), _s(_, -1), EZ);
  dive = (ey - eym) / dx[1] + (ez - ezm) / dx[2];

  cuda_sync_if_enabled();
}
