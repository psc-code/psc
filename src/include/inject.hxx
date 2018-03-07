
#pragma once

#include "psc_inject_private.h"

// ======================================================================
// InjectBase

struct InjectBase
{
  virtual void fill_ghosts(PscMfieldsBase mflds, int mb, int me) = 0;
  virtual void add_ghosts(PscMfieldsBase mflds, int mb, int me) = 0;
  virtual void reset() = 0;
};

// ======================================================================
// PscInject

template<typename S>
struct PscInject
{
  using sub_t = S;
  
  static_assert(std::is_convertible<sub_t*, InjectBase*>::value,
  		"sub classes used in PscInjectParticles must derive from InjectBase");
  
  explicit PscInject(psc_inject *inject)
    : inject_(inject)
  {}

  void operator()(PscMparticlesBase mprts, PscMfieldsBase mflds)
  {
    psc_inject_run(inject_, mprts.mprts(), mflds.mflds());
  }

private:
  psc_inject* inject_;
};

using PscInjectBase = PscInject<InjectBase>;

// ======================================================================
// PscInjectWrapper

template<typename Inject>
class PscInjectWrapper
{
public:
  const static size_t size = sizeof(Inject);
  
  static void setup(psc_inject* _inject)
  {
    PscInject<Inject> inject(_inject);
    new(inject.sub()) Inject();
  }

  static void destroy(psc_inject* _inject)
  {
    PscInject<Inject> inject(_inject);
    inject->~Inject();
  }
};

