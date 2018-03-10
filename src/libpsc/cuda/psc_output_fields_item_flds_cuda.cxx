
#include "psc_output_fields_item_private.h"
#include "psc_fields_cuda.h"
#include "fields_item.hxx"

#include "cuda_iface.h"

// ======================================================================

struct Item_dive_cuda
{
  using mfields_t = PscMfieldsCuda;
  constexpr static const char* name = "dive_cuda";
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "dive" }; } // FIXME

  static void run(mfields_t mflds, mfields_t mres)
  {
    cuda_mfields *cmflds = mflds->cmflds;
    cuda_mfields *cmres = mres->cmflds;
    
    for (int p = 0; p < mres->n_patches(); p++) {
      cuda_mfields_calc_dive_yz(cmflds, cmres, p);
    }
  }
};

FieldsItemOps<FieldsItemFields<Item_dive_cuda>> psc_output_fields_item_dive_cuda_ops;
