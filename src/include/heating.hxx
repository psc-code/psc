
#pragma once

#include "psc_heating_private.h"

#include "particles.hxx"
#include "fields3d.hxx"

#include <mrc_profile.h>

// ======================================================================
// HeatingBase

struct HeatingBase
{
  virtual void run(PscMparticlesBase mprts_base) = 0;
};

// ======================================================================
// PscHeating

template<typename S>
struct PscHeating
{
  using sub_t = S;
  
  static_assert(std::is_convertible<sub_t*, HeatingBase*>::value,
  		"sub classes used in PscHeatingParticles must derive from HeatingBase");
  
  explicit PscHeating(psc_heating *heating)
    : heating_(heating)
  {}

  void operator()(PscMparticlesBase mprts)
  {
    static int pr;
    if (!pr) {
      pr = prof_register("heating_run", 1., 0, 0);
    }  
    
    prof_start(pr);
    sub()->run(mprts);
    prof_stop(pr);
  }
  
  sub_t* sub() { return mrc_to_subobj(heating_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_heating* heating_;
};

using PscHeatingBase = PscHeating<HeatingBase>;

// ======================================================================
// PscHeatingWrapper

template<typename Heating>
class PscHeatingWrapper
{
public:
  const static size_t size = sizeof(Heating);
  
  static void setup(psc_heating* _heating)
  {
    PscHeating<Heating> heating(_heating);
    new(heating.sub()) Heating(_heating->every_step,
			       _heating->tb, _heating->te,
			       _heating->kind,
			       *_heating->spot);

    psc_heating_setup_super(_heating);
  }

  static void destroy(psc_heating* _heating)
  {
    PscHeating<Heating> heating(_heating);
    heating->~Heating();
  }
};

