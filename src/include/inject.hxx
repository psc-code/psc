
#pragma once

#include "psc_inject_private.h"

#include "particles.hxx"
#include "fields3d.hxx"

#include <mrc_profile.h>

// ======================================================================
// InjectBase

struct InjectBase
{
  struct psc_mfields *mflds_n;
  struct psc_output_fields_item *item_n;
  struct psc_bnd *item_n_bnd;
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
    static int pr;
    if (!pr) {
      pr = prof_register("inject_run", 1., 0, 0);
    }  
    
    if (!inject_->do_inject ||
	ppsc->timestep % inject_->every_step != 0) {
      return;
    }
    
    struct psc_inject_ops *ops = psc_inject_ops(inject_);
    assert(ops && ops->run);
    
    prof_start(pr);
    ops->run(inject_, mprts.mprts(), mflds.mflds());
    prof_stop(pr);
  }

  sub_t* sub() { return mrc_to_subobj(inject_, sub_t); }
  sub_t* operator->() { return sub(); }

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
    new(inject.sub()) Inject(psc_inject_comm(_inject));
    psc_inject_setup_super(_inject);
  }

  static void destroy(psc_inject* _inject)
  {
    PscInject<Inject> inject(_inject);
    inject->~Inject();
  }
};

