
#pragma once

#include "push_fields.hxx"
#include "vpic_iface.h"

// ======================================================================
// PushFieldsVpicWrap

template<typename MfieldsState>
struct PushFieldsVpicWrap
{
  using PushFieldsOps = VpicPushFieldsOps<MfieldsState>;

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

// ======================================================================
// PushFieldsVpic

template<typename MfieldsState>
struct PushFieldsVpic
{
  using PushFieldsOps = PscPushFieldsOps<MfieldsState,
					 PscFieldArrayLocalOps<MfieldsState>,
					 PscFieldArrayRemoteOps<MfieldsState>>;

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

