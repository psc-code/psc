
#pragma once

#include "psc_marder_private.h"

// ======================================================================
// MarderBase

struct MarderBase
{
  virtual void run(PscMparticlesBase mprts_base) = 0;
};

// ======================================================================
// PscMarder

template<typename S>
struct PscMarder
{
  using sub_t = S;
  
  static_assert(std::is_convertible<sub_t*, MarderBase*>::value,
  		"sub classes used in PscMarder must derive from MarderBase");
  
  explicit PscMarder(psc_marder *marder)
    : marder_(marder)
  {}
  
  void operator()(PscMparticlesBase mprts, PscMfieldsBase mflds)
  {
    psc_marder_run(marder_, mflds.mflds(), mprts.mprts());
  }
  
  sub_t* sub() { return mrc_to_subobj(marder_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_marder *marder_;
};

using PscMarderBase = PscMarder<MarderBase>;

