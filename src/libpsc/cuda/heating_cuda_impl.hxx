
#pragma once

#include "heating.hxx"

// ======================================================================
// psc_heating subclass "cuda"

struct HeatingCuda : HeatingBase
{
  // ----------------------------------------------------------------------
  // ctor

  template<typename FUNC>
  HeatingCuda(int every_step, int tb, int te, int kind, FUNC get_H)
    : every_step_(every_step),
      tb_(tb), te_(te),
      kind_(kind)
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
    foil.heating_dt = every_step * ppsc->dt;
    cuda_heating_setup_foil(&foil);
  }

  // ----------------------------------------------------------------------
  // operator()

  void operator()(MparticlesCuda& mprts)
  {
    cuda_heating_run_foil(mprts.cmprts());
  }
  
  // ----------------------------------------------------------------------
  // run

  void run(PscMparticlesBase mprts_base) override
  {
    auto& mprts = mprts_base->get_as<MparticlesCuda>();
    (*this)(mprts);
    mprts_base->put_as(mprts);
  }

private:
  int every_step_;
  int tb_, te_;
  int kind_;
};

