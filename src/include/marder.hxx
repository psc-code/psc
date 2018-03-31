
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
  
  // On ghost cells:
  // It is possible (variant = 1) that ghost cells are set before this is called
  // and the subsequent code expects ghost cells to still be set on return.
  // We're calling fill_ghosts at the end of each iteration, so that's fine.
  // However, for variant = 0, ghost cells aren't set on entry, and they're not
  // expected to be set on return (though we do that, anyway.)

  void operator()(PscMparticlesBase mprts, PscMfieldsBase mflds)
  {
    static int pr;
    if (!pr) {
      pr   = prof_register("psc_marder_run", 1., 0, 0);
    }
    
    if (marder_->every_step < 0 || ppsc->timestep % marder_->every_step != 0) 
      return;
    
    prof_start(pr);
    struct psc_marder_ops *ops = psc_marder_ops(marder_);
    assert(ops && ops->run);
    ops->run(marder_, mflds, mprts);
    prof_stop(pr);
  }
  
  sub_t* sub() { return mrc_to_subobj(marder_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_marder *marder_;
};

using PscMarderBase = PscMarder<MarderBase>;

// ======================================================================
// PscMarderWrapper

template<typename Marder>
class MarderWrapper
{
public:
  const static size_t size = sizeof(Marder);
  
  static void setup(psc_marder* _marder)
  {
    new(_marder->obj.subctx) Marder{};
  }

  static void destroy(psc_marder* _marder)
  {
    PscMarder<Marder> marder(_marder);
    marder->~Marder();
  }
};

