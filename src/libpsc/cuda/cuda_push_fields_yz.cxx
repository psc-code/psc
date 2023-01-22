
#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "fields_item_dive_cuda.hxx"
#include "cuda_bits.h"

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
  auto dive = mdive.gt().view(0, _s(bd[1], -bd[1]), _s(bd[2], -bd[2]), 0);

  auto s0 = _s(1, _);
  auto sm = _s(_, -1);
  dive = (flds.view(0, s0, s0, EY) - flds.view(0, sm, s0, EY)) / dx[1] +
         (flds.view(0, s0, s0, EZ) - flds.view(0, s0, sm, EZ)) / dx[2];

  cuda_sync_if_enabled();
}

void cuda_mfields_calc_dive_xyz(MfieldsStateCuda& mflds, MfieldsCuda& mdive)
{
  auto dx = mflds.grid().domain.dx;

  auto bnd = mflds.ibn();
  auto flds =
    mflds.gt().view(_s(-1 + bnd[0], -bnd[0]), _s(-1 + bnd[1], -bnd[1]),
                    _s(-1 + bnd[2], -bnd[2]));

  auto bd = mdive.ibn();
  auto dive =
    mdive.gt().view(_s(bd[0], -bd[0]), _s(bd[1], -bd[1]), _s(bd[2], -bd[2]), 0);

  auto s0 = _s(1, _);
  auto sm = _s(_, -1);
  dive = (flds.view(s0, s0, s0, EX) - flds.view(sm, s0, s0, EX)) / dx[0] +
         (flds.view(s0, s0, s0, EY) - flds.view(s0, sm, s0, EY)) / dx[1] +
         (flds.view(s0, s0, s0, EZ) - flds.view(s0, s0, sm, EZ)) / dx[2];

  cuda_sync_if_enabled();
}
