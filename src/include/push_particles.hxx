
#pragma once

#include "particles.hxx"
#include "fields3d.hxx"

extern int pr_time_step_no_comm; // FIXME

// ======================================================================
// PushParticlesBase

class PushParticlesBase
{
public:
  virtual void prep(MparticlesBase& mprts_base, MfieldsStateBase& mflds_base)
  { assert(0); }
  
  virtual void push_mprts(MparticlesBase& mprts_base, MfieldsStateBase& mflds_base)
  {
    const auto& grid = mprts_base.grid();
    using Bool3 = Vec3<bool>;
    Bool3 invar{grid.isInvar(0), grid.isInvar(1), grid.isInvar(2)};

    // if this function has not been overriden, call the pusher for the appropriate dims
    if (invar == Bool3{false, true, true}) { // x
      push_mprts_x(mprts_base, mflds_base);
    } else if (invar == Bool3{true, false, true}) { // y
      push_mprts_y(mprts_base, mflds_base);
    } else if (invar == Bool3{true, true, false}) { // z
      push_mprts_z(mprts_base, mflds_base);
    } else if (invar == Bool3{false, false, true}) { // xy
      push_mprts_xy(mprts_base, mflds_base);
    } else if (invar == Bool3{false, true, false}) { // xz
      push_mprts_xz(mprts_base, mflds_base);
    } else if (invar == Bool3{true, false, false}) { // yz
      push_mprts_yz(mprts_base, mflds_base);
    } else if (invar == Bool3{false, false, false}) { // xyz
      push_mprts_xyz(mprts_base, mflds_base);
    } else {
      push_mprts_1(mprts_base, mflds_base);
    }
  }

  virtual void push_mprts_xyz(MparticlesBase& mprts, MfieldsStateBase& mflds_base)
  { assert(0); }

  virtual void push_mprts_xy(MparticlesBase& mprts, MfieldsStateBase& mflds_base)
  { assert(0); }

  virtual void push_mprts_xz(MparticlesBase& mprts, MfieldsStateBase& mflds_base)
  { assert(0); }

  virtual void push_mprts_yz(MparticlesBase& mprts, MfieldsStateBase& mflds_base)
  { assert(0); }

  virtual void push_mprts_x(MparticlesBase& mprts, MfieldsStateBase& mflds_base)
  { assert(0); }

  virtual void push_mprts_y(MparticlesBase& mprts, MfieldsStateBase& mflds_base)
  { assert(0); }

  virtual void push_mprts_z(MparticlesBase& mprts, MfieldsStateBase& mflds_base)
  { assert(0); }

  virtual void push_mprts_1(MparticlesBase& mprts, MfieldsStateBase& mflds_base)
  { assert(0); }

  virtual void stagger_mprts(MparticlesBase& mprts_base, MfieldsStateBase& mflds_base)
  {
    const auto& grid = mprts_base.grid();
    using Bool3 = Vec3<bool>;
    Bool3 invar{grid.isInvar(0), grid.isInvar(1), grid.isInvar(2)};

    if (invar == Bool3{true, false, false}) { // yz
      stagger_mprts_yz(mprts_base, mflds_base);
    } else if (invar == Bool3{true, true, true}) { // 1
      stagger_mprts_1(mprts_base, mflds_base);
    } else {
      mprintf("WARNING: no stagger_mprts() case!\n");
    }
  }
  
  virtual void stagger_mprts_yz(MparticlesBase& mprts, MfieldsStateBase& mflds_base)
  { mprintf("WARNING: %s not implemented\n", __func__); }

  virtual void stagger_mprts_1(MparticlesBase& mprts, MfieldsStateBase& mflds_base)
  { mprintf("WARNING: %s not implemented\n", __func__); }

};

