
#pragma once

#include "psc_inject_private.h"

#include "particles.hxx"
#include "fields3d.hxx"

#include <mrc_profile.h>

// ======================================================================
// InjectBase

struct InjectBase
{
  InjectBase(bool do_inject, int every_step, int tau, int kind_n, psc_target* target)
    : do_inject(do_inject),
      every_step(every_step), tau(tau), kind_n(kind_n),
      target(target)
  {}

  virtual void run(PscMparticlesBase mprts_base, PscMfieldsBase mflds_base) = 0;
  
  // param
  const bool do_inject; // whether to inject particles at all
  const int every_step; // inject every so many steps
  const int tau; // in steps
  const int kind_n; // the number of particles to inject are based on this kind's density
  psc_target* const target;
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
    
    prof_start(pr);
    sub()->run(mprts, mflds);
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
    new(inject.sub()) Inject(psc_inject_comm(_inject),
			     _inject->do_inject, _inject->every_step, _inject->tau, _inject->kind_n,
			     _inject->target);

    psc_inject_setup_super(_inject);
  }

  static void destroy(psc_inject* _inject)
  {
    PscInject<Inject> inject(_inject);
    inject->~Inject();
  }
};

