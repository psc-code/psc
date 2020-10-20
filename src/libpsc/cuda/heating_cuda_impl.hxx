
#pragma once

#include "heating.hxx"
#include "mparticles_cuda.hxx"
#include "psc_fields_cuda.h"
#include "heating_spot_foil.hxx"

// ======================================================================
// psc_heating subclass "cuda"

template <typename HS>
struct cuda_heating_foil;

template <typename HS, typename MP>
struct HeatingCuda : HeatingBase
{
  HeatingCuda(const Grid_t& grid, int interval, HS heating_spot);

  ~HeatingCuda();

  void reset(const MP& mprts);

  void operator()(MP& mprts);

private:
  cuda_heating_foil<HS>* foil_;
  int balance_generation_cnt_;
};
