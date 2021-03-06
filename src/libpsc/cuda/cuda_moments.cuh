
#pragma once

#include "psc_fields_cuda.h"

template <typename CudaMparticles, typename dim>
struct CudaMoments1stNcRho
{
  void operator()(CudaMparticles& cmprts, MfieldsCuda& mres);
};

template <typename CudaMparticles, typename dim>
struct CudaMoments1stN
{
  void operator()(CudaMparticles& cmprts, MfieldsCuda& mres);
};

template <typename CudaMparticles, typename dim>
struct CudaMoments1stAll
{
  void operator()(CudaMparticles& cmprts, MfieldsCuda& mres);
};
