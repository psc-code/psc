
#include "psc_fields_vpic.h"

#include "psc_fields_c.h"
#include "psc_fields_single.h"
#include "psc.h"
#include "psc_method.h"
#include "psc_particles_vpic.h"
#include "fields.hxx"

#include "mrc_domain.h"
#include "mrc_bits.h"

#include "vpic_iface.h"

using Fields = Fields3d<fields_vpic_t>;
using FieldsS = Fields3d<fields_single_t>;

static const int map_psc2vpic[VPIC_MFIELDS_N_COMP] = {
  [JXI] = VPIC_MFIELDS_JFX, [JYI] = VPIC_MFIELDS_JFY, [JZI] = VPIC_MFIELDS_JFZ,
  [EX]  = VPIC_MFIELDS_EX , [EY]  = VPIC_MFIELDS_EY , [EZ]  = VPIC_MFIELDS_EZ,
  [HX]  = VPIC_MFIELDS_BX , [HY]  = VPIC_MFIELDS_BY , [HZ]  = VPIC_MFIELDS_BZ,

  [9]   = VPIC_MFIELDS_TCAX, [10]  = VPIC_MFIELDS_TCAY, [11]  = VPIC_MFIELDS_TCAZ,
  [12]  = VPIC_MFIELDS_DIV_E_ERR, [13] = VPIC_MFIELDS_DIV_B_ERR,
  [14]  = VPIC_MFIELDS_RHOB     , [15] = VPIC_MFIELDS_RHOF,
  [16]  = 16, [17] = 17, [18] = 18, [19] = 19,
};

// ======================================================================
// convert from/to "single"

static void MfieldsVpic_copy_from_single(MfieldsBase& mflds, MfieldsBase& mflds_single, int mb, int me)
{
  auto& mf_single = dynamic_cast<MfieldsSingle&>(mflds_single);
  auto& mf = dynamic_cast<MfieldsVpic&>(mflds);
  for (int p = 0; p < mf.n_patches(); p++) {
    fields_vpic_t flds = mf[p];
    FieldsS F_s(mf_single[p]);

    
    // FIXME, hacky way to distinguish whether we want
    // to copy the field into the standard PSC component numbering
    if (mf.n_comps() == VPIC_MFIELDS_N_COMP) {
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
    } else {
      if (mb != me) {
	assert(0);
      }
    }
  }
}

static void MfieldsVpic_copy_to_single(MfieldsBase& mflds, MfieldsBase& mflds_single, int mb, int me)
{
  auto& mf_single = dynamic_cast<MfieldsSingle&>(mflds_single);
  auto& mf = dynamic_cast<MfieldsVpic&>(mflds);
  for (int p = 0; p < mf.n_patches(); p++) {
    fields_vpic_t flds = mf[p];
    fields_single_t flds_single = mf_single[p];
    FieldsS F_s(flds_single);

    int ib[3], ie[3];
    for (int d = 0; d < 3; d++) {
      ib[d] = MAX(flds.ib_[d], flds_single.ib_[d]);
      ie[d] = MIN(flds.ib_[d] + flds.im_[d], flds_single.ib_[d] + flds_single.im_[d]);
    }

    // FIXME, hacky way to distinguish whether we want
    // to copy the field into the standard PSC component numbering or,
    // as in this case, just copy one-to-one
    if (mf.n_comps() == VPIC_HYDRO_N_COMP) {
      for (int m = mb; m < me; m++) {
	for (int jz = ib[2]; jz < ie[2]; jz++) {
	  for (int jy = ib[1]; jy < ie[1]; jy++) {
	    for (int jx = ib[0]; jx < ie[0]; jx++) {
	      F_s(m, jx,jy,jz) = flds(m, jx,jy,jz);
	    }
	  }
	}
      }
    } else if (mf.n_comps() == VPIC_MFIELDS_N_COMP) {
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
    } else {
      assert(0);
    }
  }
}

static void MfieldsStateVpic_copy_from_single(MfieldsBase& mflds, MfieldsBase& mflds_single, int mb, int me)
{
  auto& mf_single = dynamic_cast<MfieldsSingle&>(mflds_single);
  auto& mf = dynamic_cast<MfieldsStateVpic&>(mflds);
  for (int p = 0; p < mf.n_patches(); p++) {
    fields_vpic_t flds = mf[p];
    FieldsS F_s(mf_single[p]);

    
    // FIXME, hacky way to distinguish whether we want
    // to copy the field into the standard PSC component numbering
    if (mf.n_comps() == VPIC_MFIELDS_N_COMP) {
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
    } else {
      if (mb != me) {
	assert(0);
      }
    }
  }
}

static void MfieldsStateVpic_copy_to_single(MfieldsBase& mflds, MfieldsBase& mflds_single, int mb, int me)
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

    // FIXME, hacky way to distinguish whether we want
    // to copy the field into the standard PSC component numbering or,
    // as in this case, just copy one-to-one
    if (mf.n_comps() == VPIC_HYDRO_N_COMP) {
      for (int m = mb; m < me; m++) {
	for (int jz = ib[2]; jz < ie[2]; jz++) {
	  for (int jy = ib[1]; jy < ie[1]; jy++) {
	    for (int jx = ib[0]; jx < ie[0]; jx++) {
	      F_s(m, jx,jy,jz) = flds(m, jx,jy,jz);
	    }
	  }
	}
      }
    } else if (mf.n_comps() == VPIC_MFIELDS_N_COMP) {
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
    } else {
      assert(0);
    }
  }
}

// ======================================================================

static int ref_count_fields, ref_count_hydro;

// ----------------------------------------------------------------------
// MfieldsVpic ctor

MfieldsVpic::MfieldsVpic(const Grid_t& grid, int n_fields, Int3 ibn)
  : MfieldsBase(grid, n_fields, ibn)
{
  assert(grid.n_patches() == 1);
  assert((ibn == Int3{ 1, 1, 1 }));

  psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim_);

  if (n_fields == VPIC_MFIELDS_N_COMP) {
    // make sure we notice if we create a second psc_mfields
    // which would share its memory with the first
    if (ref_count_fields != 0) {
      MHERE;
    }
    ref_count_fields++;

    vmflds_fields_ = sim_->field_array_;
  } else if (n_fields == VPIC_HYDRO_N_COMP) {
    // make sure we notice if we create a second psc_mfields
    // which would share its memory with the first
    assert(ref_count_hydro == 0);
    ref_count_hydro++;

    vmflds_hydro = sim_->hydro_array_;
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// MfieldsVpic dtor

MfieldsVpic::~MfieldsVpic()
{
  if (n_fields_ == VPIC_MFIELDS_N_COMP) {
    ref_count_fields--;
  } else if (n_fields_ == VPIC_HYDRO_N_COMP) {
    ref_count_hydro--;
  } else {
    assert(0);
  }
}

const MfieldsBase::Convert MfieldsVpic::convert_to_ = {
  { std::type_index(typeid(MfieldsSingle)), MfieldsVpic_copy_to_single },
};

const MfieldsBase::Convert MfieldsVpic::convert_from_ = {
  { std::type_index(typeid(MfieldsSingle)), MfieldsVpic_copy_from_single },
};

const MfieldsBase::Convert MfieldsStateVpic::convert_to_ = {
  { std::type_index(typeid(MfieldsSingle)), MfieldsStateVpic_copy_to_single },
};

const MfieldsBase::Convert MfieldsStateVpic::convert_from_ = {
  { std::type_index(typeid(MfieldsSingle)), MfieldsStateVpic_copy_from_single },
};

