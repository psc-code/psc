
#pragma once

#include "push_fields.hxx"
#include "vpic_iface.h"

// ======================================================================
// PushFieldsVpic

template<typename PushFieldsOps>
struct PushFieldsVpic_ : PushFieldsBase
{
  using MfieldsState = typename PushFieldsOps::MfieldsState;
  
  void push_E(MfieldsState& mflds, double dt_fac)
  {
    TIC PushFieldsOps::advance_e(mflds, dt_fac); TOC(advance_e, 1);
#if 0
    user_field_injection();
#endif
  }

  void push_H(MfieldsState& mflds, double dt_fac)
  {
    TIC PushFieldsOps::advance_b(mflds, dt_fac); TOC(advance_b, 1);
  }
};

template<typename MfieldsState>
using PushFieldsVpicWrap = PushFieldsVpic_<VpicPushFieldsOps<MfieldsState>>;

template<typename MfieldsState>
using PushFieldsVpic = PushFieldsVpic_<PscPushFieldsOps<MfieldsState, PscFieldArrayLocalOps<MfieldsState>,
							PscFieldArrayRemoteOps<MfieldsState>>>;

