
#pragma once

#include "heating.hxx"

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
    struct cuda_heating_foil foil;
    double val;

    auto spot = get_H;
    foil.zl = spot.zl;
    foil.zh = spot.zh;
    foil.xc = spot.xc;
    foil.yc = spot.yc;
    foil.rH = spot.rH;
    foil.T  = spot.T;
    foil.Mi = spot.Mi;
    foil.kind = kind;
    foil.heating_dt = interval * ppsc->dt;
    cuda_heating_setup_foil(&foil);
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

