
#include "psc.h"
#include "psc_fields_as_single.h"
#include "psc_fields_c.h"
#include "fields.hxx"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define PFX(x) psc_fields_single_ ## x
#define MPFX(x) psc_mfields_single_ ## x
#define MFIELDS MfieldsSingle

using Fields = Fields3d<MfieldsSingle::fields_t>;
using FieldsC = Fields3d<fields_c_t>;

// ======================================================================
// convert to c

template<typename MfieldsBase, typename MfieldsSingle, typename MfieldsC>
static void psc_mfields_single_copy_from_c(MfieldsBase& mflds, MfieldsBase& mflds_c, int mb, int me)
{
  auto& mf = dynamic_cast<MfieldsSingle&>(mflds);
  auto& mf_c = dynamic_cast<MfieldsC&>(mflds_c);
  for (int p = 0; p < mf.n_patches(); p++) {
    auto flds = mf[p];
    Fields F(flds);
    FieldsC F_c(mf_c[p]);
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib_[2]; jz < flds.ib_[2] + flds.im_[2]; jz++) {
	for (int jy = flds.ib_[1]; jy < flds.ib_[1] + flds.im_[1]; jy++) {
	  for (int jx = flds.ib_[0]; jx < flds.ib_[0] + flds.im_[0]; jx++) {
	    F(m, jx,jy,jz) = F_c(m, jx,jy,jz);
	  }
	}
      }
    }
  }
}

template<typename MfieldsBase, typename MfieldsSingle, typename MfieldsC>
static void psc_mfields_single_copy_to_c(MfieldsBase& mflds, MfieldsBase& mflds_c, int mb, int me)
{
  auto& mf = dynamic_cast<MfieldsSingle&>(mflds);
  auto& mf_c = dynamic_cast<MfieldsC&>(mflds_c);
  for (int p = 0; p < mf.n_patches(); p++) {
    auto flds = mf[p];
    Fields F(flds);
    FieldsC F_c(mf_c[p]);
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib_[2]; jz < flds.ib_[2] + flds.im_[2]; jz++) {
	for (int jy = flds.ib_[1]; jy < flds.ib_[1] + flds.im_[1]; jy++) {
	  for (int jx = flds.ib_[0]; jx < flds.ib_[0] + flds.im_[0]; jx++) {
	    F_c(m, jx,jy,jz) = F(m, jx,jy,jz);
	  }
	}
      }
    }
  }
}

template<> const MfieldsBase::Convert MfieldsSingle::convert_to_ = {
  { std::type_index(typeid(MfieldsC)), psc_mfields_single_copy_to_c<MfieldsBase, MfieldsSingle, MfieldsC> },
};

template<> const MfieldsBase::Convert MfieldsSingle::convert_from_ = {
  { std::type_index(typeid(MfieldsC)), psc_mfields_single_copy_from_c<MfieldsBase, MfieldsSingle, MfieldsC> },
};

template<> const MfieldsStateBase::Convert MfieldsStateSingle::convert_to_ = {
  { std::type_index(typeid(MfieldsStateDouble)), psc_mfields_single_copy_to_c<MfieldsStateBase, MfieldsStateSingle, MfieldsStateDouble> },
};

template<> const MfieldsStateBase::Convert MfieldsStateSingle::convert_from_ = {
  { std::type_index(typeid(MfieldsStateDouble)), psc_mfields_single_copy_from_c<MfieldsStateBase, MfieldsStateSingle, MfieldsStateDouble> },
};

#include "psc_fields_common.cxx"

