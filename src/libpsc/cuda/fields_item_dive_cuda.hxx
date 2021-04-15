
#include <psc_fields_cuda.h>
#include <fields_item.hxx>

#include "cuda_iface.h"

void cuda_mfields_calc_dive_yz(MfieldsStateCuda& mflds, MfieldsCuda& mf, int p);

template <>
struct Item_dive<MfieldsStateCuda>
{
  using Mfields = MfieldsCuda;
  using MfieldsState = MfieldsStateCuda;
  constexpr static const char* name = "dive";
  constexpr static int n_comps = 1;
  static std::vector<std::string> fld_names() { return {"dive"}; } // FIXME

  Item_dive(MfieldsState& mflds)
    : mres_{mflds.grid(), n_comps, mflds.grid().ibn}
  {
    assert(mflds.grid().isInvar(0));
    for (int p = 0; p < mflds.n_patches(); p++) {
      cuda_mfields_calc_dive_yz(mflds, mres_, p);
    }
  }

  Mfields& result() { return mres_; }

  auto gt()
  {
    auto bnd = mres_.ibn();
    return mres_.gt().view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]),
                           _s(bnd[2], -bnd[2]));
  }

private:
  Mfields mres_;
};
