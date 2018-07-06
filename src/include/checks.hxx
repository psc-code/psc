
#pragma once

#include "psc_checks_private.h"

// ======================================================================
// class ChecksBase

class ChecksBase
{
public:
  virtual void continuity_before_particle_push(psc* psc, MparticlesBase& mprts) = 0;
  virtual void continuity_after_particle_push(psc* psc, MparticlesBase& mprts, MfieldsBase& mflds) = 0;
  virtual void gauss(psc* psc, MparticlesBase& mprts, MfieldsBase& mflds) = 0;
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

  void continuity_before_particle_push(psc* psc, MparticlesBase& mprts)
  {
    sub()->continuity_before_particle_push(psc, mprts);
  }

  void continuity_after_particle_push(psc* psc, MparticlesBase& mprts, MfieldsBase& mflds)
  { 
    sub()->continuity_after_particle_push(psc, mprts, mflds);
  }

  void gauss(psc* psc, MparticlesBase& mprts, MfieldsBase& mflds)
  {
    sub()->gauss(psc, mprts, mflds);
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
    new(checks.sub()) Checks_t{ppsc->grid(), psc_checks_comm(_checks), _checks->params};
  }

  static void destroy(struct psc_checks *_checks)
  {
    PscChecks<Checks_t> checks(_checks);
    checks->~Checks_t();
  }
};

