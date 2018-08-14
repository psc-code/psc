
#pragma once

#include "push_fields.hxx"
#include "vpic_iface.h"
#include "psc_method.h"
#include "psc_fields_vpic.h"

// ======================================================================
// PushFieldsVpic

struct PushFieldsVpic : PushFieldsBase
{
  PushFieldsVpic()
  {
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim_);
  }

  void push_E(MfieldsState& mflds, double dt_fac)
  {
    TIC PushFieldsOps::advance_e(mflds.vmflds(), dt_fac); TOC(advance_e, 1);
    sim_->field_injection();
  }

  void push_H(MfieldsState& mflds, double dt_fac)
  {
    TIC PushFieldsOps::advance_b(mflds.vmflds(), dt_fac); TOC(advance_b, 1);
  }

  void push_E(MfieldsStateBase& mflds_base, double dt_fac) override
  {
    assert(0);
#if 0
    // needs J, E, B, TCA, material
    auto& mflds = mflds_base.get_as<MfieldsStateVpic>(JXI, MfieldsStateVpic::N_COMP);
    push_E(mflds, dt_fac);
    // updates E, TCA, and B ghost points FIXME 9 == TCAX
    mflds_base.put_as(mflds, EX, 9 + 3);
#endif
  }

  void push_H(MfieldsStateBase& mflds_base, double dt_fac) override
  {
    assert(0);
#if 0
    // needs E, B
    auto& mflds = mflds_base.get_as<MfieldsStateVpic>(EX, HX + 6);
    push_H(mflds, dt_fac);
    // updates B
    mflds_base.put_as(mflds, HX, HX + 3);
#endif
  }

private:
  Simulation *sim_;
};

