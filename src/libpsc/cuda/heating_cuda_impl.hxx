
#pragma once

#include "heating.hxx"
#include "mparticles_cuda.hxx"
#include "psc_fields_cuda.h"
#include "heating_spot_foil.hxx"

// ======================================================================
// psc_heating subclass "cuda"

template <typename HS>
struct cuda_heating_foil;

template <typename BS>
struct HeatingCuda : HeatingBase
{
  template <typename FUNC>
  HeatingCuda(const Grid_t& grid, int interval, FUNC get_H);

  ~HeatingCuda();

  void reset(MparticlesCuda<BS>& mprts);

  void operator()(MparticlesCuda<BS>& mprts);

private:
  cuda_heating_foil<HeatingSpotFoil<dim_xyz>>* foil_;
  int balance_generation_cnt_;
};
