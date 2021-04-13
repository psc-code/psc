
#include "psc.h"
#include "psc_fields_as_single.h"
#include "psc_fields_c.h"
#include "fields.hxx"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define PFX(x) psc_fields_single_##x
#define MPFX(x) psc_mfields_single_##x
#define MFIELDS MfieldsSingle

// ======================================================================
// convert to c

template <typename MfieldsBase, typename MfieldsSingle, typename MfieldsC>
static void psc_mfields_single_copy_from_c(MfieldsBase& mflds,
                                           MfieldsBase& mflds_c, int mb, int me)
{
  auto& mf = dynamic_cast<MfieldsSingle&>(mflds);
  auto& mf_c = dynamic_cast<MfieldsC&>(mflds_c);
  for (int p = 0; p < mf.n_patches(); p++) {
    auto flds = make_Fields3d<dim_xyz>(mf[p]);
    auto flds_c = make_Fields3d<dim_xyz>(mf_c[p]);
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib()[2]; jz < flds.ib()[2] + flds.shape(2); jz++) {
        for (int jy = flds.ib()[1]; jy < flds.ib()[1] + flds.shape(1); jy++) {
          for (int jx = flds.ib()[0]; jx < flds.ib()[0] + flds.shape(0); jx++) {
            flds(m, jx, jy, jz) = flds_c(m, jx, jy, jz);
          }
        }
      }
    }
  }
}

template <typename MfieldsBase, typename MfieldsSingle, typename MfieldsC>
static void psc_mfields_single_copy_to_c(MfieldsBase& mflds,
                                         MfieldsBase& mflds_c, int mb, int me)
{
  auto& mf = dynamic_cast<MfieldsSingle&>(mflds);
  auto& mf_c = dynamic_cast<MfieldsC&>(mflds_c);
  for (int p = 0; p < mf.n_patches(); p++) {
    auto flds = make_Fields3d<dim_xyz>(mf[p]);
    auto flds_c = make_Fields3d<dim_xyz>(mf_c[p]);
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib()[2]; jz < flds.ib()[2] + flds.shape(2); jz++) {
        for (int jy = flds.ib()[1]; jy < flds.ib()[1] + flds.shape(1); jy++) {
          for (int jx = flds.ib()[0]; jx < flds.ib()[0] + flds.shape(0); jx++) {
            flds_c(m, jx, jy, jz) = flds(m, jx, jy, jz);
          }
        }
      }
    }
  }
}

template <>
const MfieldsBase::Convert MfieldsSingle::convert_to_ = {
  {std::type_index(typeid(MfieldsC)),
   psc_mfields_single_copy_to_c<MfieldsBase, MfieldsSingle, MfieldsC>},
};

template <>
const MfieldsBase::Convert MfieldsSingle::convert_from_ = {
  {std::type_index(typeid(MfieldsC)),
   psc_mfields_single_copy_from_c<MfieldsBase, MfieldsSingle, MfieldsC>},
};

template <>
const MfieldsStateBase::Convert MfieldsStateSingle::convert_to_ = {
  {std::type_index(typeid(MfieldsStateDouble)),
   psc_mfields_single_copy_to_c<MfieldsStateBase, MfieldsStateSingle,
                                MfieldsStateDouble>},
};

template <>
const MfieldsStateBase::Convert MfieldsStateSingle::convert_from_ = {
  {std::type_index(typeid(MfieldsStateDouble)),
   psc_mfields_single_copy_from_c<MfieldsStateBase, MfieldsStateSingle,
                                  MfieldsStateDouble>},
};

#include "psc_fields_common.cxx"
