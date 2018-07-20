
#pragma once

#include "push_particles.hxx"
#include "vpic_iface.h"
#include "psc_method.h"

#include "psc_fields_vpic.h"
#include "psc_particles_vpic.h"

// ======================================================================
// PushParticlesVpic

struct PushParticlesVpic : PushParticlesBase
{
  PushParticlesVpic()
  {
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim_);
  }

  void prep(MparticlesBase& _mprts_base, MfieldsBase& mflds_base) override
  {
    // needs E, B
    auto& mflds = mflds_base.get_as<MfieldsVpic>(EX, HX + 6);
    sim_->push_mprts_prep(*mflds.vmflds_fields);
    mflds_base.put_as(mflds, 0, 0);
  }
  
  void push_mprts(MparticlesBase& mprts_base, MfieldsBase& mflds_base) override
  {
    // needs E, B (not really, because they're already in interpolator), rhob?
    auto& mflds = mflds_base.get_as<MfieldsVpic>(EX, HX + 6);
    auto& mprts = mprts_base.get_as<MparticlesVpic>();
    
    sim_->push_mprts(*mprts.vmprts, *mflds.vmflds_fields);
    
    // update jf FIXME: rhob too, probably, depending on b.c.
    mflds_base.put_as(mflds, JXI, JXI + 3);
    mprts_base.put_as(mprts);
  }

private:
  Simulation* sim_;
};

