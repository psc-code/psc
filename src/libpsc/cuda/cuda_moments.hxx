
#pragma once

#include "psc_fields_cuda.h"

template <typename CudaMparticles, typename dim>
struct CudaMoments1stNcRho
{
  void operator()(CudaMparticles& cmprts, MfieldsCuda::Storage& mres_gt,
                  const Int3& mres_ib);
};

template <typename CudaMparticles, typename dim>
struct CudaMoments1stN
{
  void operator()(CudaMparticles& cmprts, MfieldsCuda::Storage& mres_gt,
                  const Int3& mres_ib);
};

template <typename CudaMparticles, typename dim>
struct CudaMoments1stAll
{
  void operator()(CudaMparticles& cmprts, MfieldsCuda::Storage& mres_gt,
                  const Int3& mres_ib);
};
