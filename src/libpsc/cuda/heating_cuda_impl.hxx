
#pragma once

#include "heating.hxx"
#include "cuda_heating_iface.hxx"

// ======================================================================
// psc_heating subclass "cuda"

template<typename BS>
struct HeatingCuda : HeatingBase
{
  // ----------------------------------------------------------------------
  // ctor

  template<typename FUNC>
  HeatingCuda(int interval, int kind, FUNC get_H)
    : kind_(kind)
  {
    double val;

    cuda_heating_foil foil{get_H, kind, interval * ppsc->dt};
  }

  // ----------------------------------------------------------------------
  // operator()

  void operator()(MparticlesCuda<BS>& mprts)
  {
    cuda_heating_run_foil(mprts.cmprts());
  }
  
  // ----------------------------------------------------------------------
  // run

  void run(PscMparticlesBase mprts_base) override
  {
    auto& mprts = mprts_base->get_as<MparticlesCuda<BS>>();
    (*this)(mprts);
    mprts_base->put_as(mprts);
  }

private:
  int kind_;
};

