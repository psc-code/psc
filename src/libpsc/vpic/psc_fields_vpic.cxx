
#include "psc_fields_c.h"
#include "psc_fields_single.h"
#include "psc.h"
#include "fields.hxx"

#include "mrc_domain.h"
#include "mrc_bits.h"

#include "vpic_iface.h"

#if 0

using Fields = Fields3d<MfieldsStateVpic::fields_t>;
using FieldsS = Fields3d<fields_single_t>;

static const int map_psc2vpic[MfieldsStateVpic::N_COMP] = {
  [JXI] = MfieldsStateVpic::JFX, [JYI] = MfieldsStateVpic::JFY, [JZI] = MfieldsStateVpic::JFZ,
  [EX]  = MfieldsStateVpic::EX , [EY]  = MfieldsStateVpic::EY , [EZ]  = MfieldsStateVpic::EZ,
  [HX]  = MfieldsStateVpic::BX , [HY]  = MfieldsStateVpic::BY , [HZ]  = MfieldsStateVpic::BZ,

  [9]   = MfieldsStateVpic::TCAX, [10]  = MfieldsStateVpic::TCAY, [11]  = MfieldsStateVpic::TCAZ,
  [12]  = MfieldsStateVpic::DIV_E_ERR, [13] = MfieldsStateVpic::DIV_B_ERR,
  [14]  = MfieldsStateVpic::RHOB     , [15] = MfieldsStateVpic::RHOF,
  [16]  = 16, [17] = 17, [18] = 18, [19] = 19,
};

// ======================================================================
// convert from/to "single"

static void MfieldsStateVpic_copy_from_single(MfieldsStateBase& mflds, MfieldsStateBase& mflds_single, int mb, int me)
{
  auto& mf_single = dynamic_cast<MfieldsSingle&>(mflds_single);
  auto& mf = dynamic_cast<MfieldsStateVpic&>(mflds);
  for (int p = 0; p < mf.n_patches(); p++) {
    fields_vpic_t flds = mf[p];
    FieldsS F_s(mf_single[p]);

    assert(mf.n_comps() == MfieldsStateVpic::N_COMP);
    for (int m = mb; m < me; m++) {
      int m_vpic = map_psc2vpic[m];
      for (int jz = flds.ib_[2]; jz < flds.ib_[2] + flds.im_[2]; jz++) {
	for (int jy = flds.ib_[1]; jy < flds.ib_[1] + flds.im_[1]; jy++) {
	  for (int jx = flds.ib_[0]; jx < flds.ib_[0] + flds.im_[0]; jx++) {
	    flds(m_vpic, jx,jy,jz) = F_s(m, jx,jy,jz);
	  }
	}
      }
    }
  }
}

static void MfieldsStateVpic_copy_to_single(MfieldsStateBase& mflds, MfieldsStateBase& mflds_single, int mb, int me)
{
  auto& mf_single = dynamic_cast<MfieldsSingle&>(mflds_single);
  auto& mf = dynamic_cast<MfieldsStateVpic&>(mflds);
  for (int p = 0; p < mf.n_patches(); p++) {
    fields_vpic_t flds = mf[p];
    fields_single_t flds_single = mf_single[p];
    FieldsS F_s(flds_single);

    int ib[3], ie[3];
    for (int d = 0; d < 3; d++) {
      ib[d] = MAX(flds.ib_[d], flds_single.ib_[d]);
      ie[d] = MIN(flds.ib_[d] + flds.im_[d], flds_single.ib_[d] + flds_single.im_[d]);
    }

    assert(mf.n_comps() == MfieldsStateVpic::N_COMP);
    for (int m = mb; m < me; m++) {
      int m_vpic = map_psc2vpic[m];
      for (int jz = ib[2]; jz < ie[2]; jz++) {
	for (int jy = ib[1]; jy < ie[1]; jy++) {
	  for (int jx = ib[0]; jx < ie[0]; jx++) {
	    F_s(m, jx,jy,jz) = flds(m_vpic, jx,jy,jz);
	  }
	}
      }
    }
  }
}

// ======================================================================

const MfieldsStateBase::Convert MfieldsStateVpic::convert_to_ = {
  { std::type_index(typeid(MfieldsStateSingle)), MfieldsStateVpic_copy_to_single },
};

const MfieldsStateBase::Convert MfieldsStateVpic::convert_from_ = {
  { std::type_index(typeid(MfieldsStateSingle)), MfieldsStateVpic_copy_from_single },
};

#endif
