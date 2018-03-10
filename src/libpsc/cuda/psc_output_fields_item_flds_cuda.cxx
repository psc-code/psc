
#include "psc_output_fields_item_private.h"
#include "psc_fields_cuda.h"
#include "fields_item.hxx"

#include "cuda_iface.h"

// ======================================================================

struct FieldsItem_dive_cuda : FieldsItemCRTP<FieldsItem_dive_cuda>
{
  using Base = FieldsItemCRTP<FieldsItem_dive_cuda>;
  using Base::Base;
  
  static const char* name() { return "dive_cuda"; }
  constexpr static int n_comps = 1;
  constexpr static fld_names_t fld_names() { return { "dive" }; } // FIXME
  constexpr static int flags = 0;

  void run(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) override
  {
    assert(ppsc->domain.gdims[0] == 1);

    PscMfieldsCuda mf = mflds_base.get_as<PscMfieldsCuda>(EX, EX+3);
    PscMfieldsCuda mf_res = mres_base_->get_as<PscMfieldsCuda>(0, 0);
    cuda_mfields *cmflds = mf->cmflds;
    cuda_mfields *cmres = mf_res->cmflds;
    
    for (int p = 0; p < mf_res->n_patches(); p++) {
      cuda_mfields_calc_dive_yz(cmflds, cmres, p);
    }
    
    mf.put_as(mflds_base, 0, 0);
    mf_res.put_as(mres_base_, 0, 1);
  }
};

FieldsItemOps<FieldsItem_dive_cuda> psc_output_fields_item_dive_cuda_ops;

