
#pragma once

#include "heating.hxx"

// ======================================================================
// psc_heating subclass "cuda"

struct cuda_heating_foil;

template<typename BS>
struct HeatingCuda : HeatingBase
{
  template<typename FUNC>
  HeatingCuda(const Grid_t& grid, int interval, int kind, FUNC get_H);

  ~HeatingCuda();

  void operator()(MparticlesCuda<BS>& mprts);
  
  // ----------------------------------------------------------------------
  // run

  void run(MparticlesBase& mprts_base) override
  {
    auto& mprts = mprts_base.get_as<MparticlesCuda<BS>>();
    (*this)(mprts);
    mprts_base.put_as(mprts);
  }

private:
  cuda_heating_foil* foil_;
};

