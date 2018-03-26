
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
    sub()->continuity_before_particle_push(psc);
  }

  void continuity_after_particle_push(psc* psc)
  { 
    sub()->continuity_after_particle_push(psc);
  }

  void gauss(psc* psc)
  {
    sub()->gauss(psc);
  }

  sub_t* operator->() { return sub(); }

  sub_t* sub() { return mrc_to_subobj(checks_, sub_t); }

private:
  psc_checks* checks_;
};

using PscChecksBase = PscChecks<ChecksBase>;

// ======================================================================
// ChecksWrapper

template<typename Checks_t>
class ChecksWrapper
{
public:
  const static size_t size = sizeof(Checks_t);
  
  static void setup(struct psc_checks *_checks)
  {
    PscChecks<Checks_t> checks(_checks);
    new(checks.sub()) Checks_t{psc_checks_comm(_checks), _checks->params};
  }

  static void destroy(struct psc_checks *_checks)
  {
    PscChecks<Checks_t> checks(_checks);
    checks->~Checks_t();
  }
};

