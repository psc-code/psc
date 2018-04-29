
#pragma once

#include "heating.hxx"
#include "cuda_heating_iface.hxx"

// ======================================================================
// psc_heating subclass "cuda"

template<typename BS>
struct HeatingCuda : HeatingBase
{
  template<typename FUNC>
  HeatingCuda(int interval, int kind, FUNC get_H);

  ~HeatingCuda();

  void operator()(MparticlesCuda<BS>& mprts);
  
  // ----------------------------------------------------------------------
  // run

  void run(PscMparticlesBase mprts_base) override
  {
    auto& mprts = mprts_base->get_as<MparticlesCuda<BS>>();
    (*this)(mprts);
    mprts_base->put_as(mprts);
  }

private:
  cuda_heating_foil* foil_;
};

