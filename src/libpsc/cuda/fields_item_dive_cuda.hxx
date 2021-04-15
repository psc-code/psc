
#include <psc_fields_cuda.h>
#include <fields_item.hxx>

#include "cuda_iface.h"

struct Item_dive_cuda
{
  using Mfields = MfieldsCuda;
  using MfieldsState = MfieldsStateCuda;
  constexpr static const char* name = "dive";
  constexpr static int n_comps = 1;
  static std::vector<std::string> fld_names() { return {"dive"}; } // FIXME

  Item_dive_cuda(MfieldsState& mflds)
    : mres_{mflds.grid(), n_comps, mflds.grid().ibn}
  {
    for (int p = 0; p < mflds.n_patches(); p++) {
      cuda_mfields_calc_dive_yz(mflds.cmflds(), mres_.cmflds(), p);
    }
  }

  Mfields& result() { return mres_; }

private:
  Mfields mres_;
};
