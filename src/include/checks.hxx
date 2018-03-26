
#pragma once

#include "psc_checks_private.h"

// ======================================================================
// class ChecksBase

class ChecksBase
{
public:
  virtual void continuity_before_particle_push(psc* psc) = 0;
  virtual void continuity_after_particle_push(psc* psc) = 0;
  virtual void gauss(psc* psc) = 0;
};

// ======================================================================
// PscChecks

template<typename S>
struct PscChecks
{
  using sub_t = S;

  static_assert(std::is_convertible<sub_t*, ChecksBase*>::value,
  		"sub classes used in PscChecksParticles must derive from ChecksBase");
  
  explicit PscChecks(psc_checks* checks)
    : checks_(checks)
  {}

  void continuity_before_particle_push(psc* psc)
  {
    struct psc_checks_ops *ops = psc_checks_ops(checks_);
    assert(ops && ops->continuity_before_particle_push);
    ops->continuity_before_particle_push(checks_, psc);
  }

  void continuity_after_particle_push(psc* psc)
  { 
    struct psc_checks_ops *ops = psc_checks_ops(checks_);
    assert(ops && ops->continuity_after_particle_push);
    ops->continuity_after_particle_push(checks_, psc);
  }

  void gauss(psc* psc)
  {
    struct psc_checks_ops *ops = psc_checks_ops(checks_);
    assert(ops && ops->gauss);
    ops->gauss(checks_, psc);
  }

  sub_t* operator->() { return sub(); }

  sub_t* sub() { return mrc_to_subobj(checks_, sub_t); }

private:
  psc_checks* checks_;
};

using PscChecksBase = PscChecks<ChecksBase>;

