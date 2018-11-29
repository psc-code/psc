
#pragma once

#include "push_fields.hxx"
#include "vpic_iface.h"

#ifdef USE_VPIC

// ======================================================================
// PushFieldsVpicWrap

template<typename MfieldsState>
struct PushFieldsVpicWrap
{
  void push_E(MfieldsState& mflds, double dt_fac)
  {
    field_array_t* fa = mflds;
    TIC return fa->kernel->advance_e(fa, dt_fac); TOC(advance_e, 1);
#if 0
    user_field_injection();
#endif
  }

  void push_H(MfieldsState& mflds, double dt_fac)
  {
    field_array_t* fa = mflds;
    TIC return fa->kernel->advance_b(fa, dt_fac); TOC(advance_b, 1);
  }
};

#endif

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

