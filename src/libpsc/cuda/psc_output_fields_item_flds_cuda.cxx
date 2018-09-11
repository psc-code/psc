
#include "psc_output_fields_item_private.h"
#include "psc_fields_cuda.h"
#include "fields_item.hxx"

#include "cuda_iface.h"

// ======================================================================

struct Item_dive_cuda
{
  using Mfields = MfieldsCuda;
  using MfieldsState = MfieldsStateCuda;
  constexpr static const char* name = "dive";
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "dive" }; } // FIXME

  static void run(MfieldsState& mflds, Mfields& mres)
  {
    for (int p = 0; p < mres.n_patches(); p++) {
      cuda_mfields_calc_dive_yz(mflds.cmflds(), mres.cmflds(), p);
    }
  }
};

FieldsItemOps<FieldsItemFields<Item_dive_cuda>> psc_output_fields_item_dive_cuda_ops;
