
#pragma once

template <typename CudaMparticles, typename dim>
struct CudaMoments1stNcRho
{
  void operator()(CudaMparticles& cmprts, struct cuda_mfields* cmres);
};

template <typename CudaMparticles, typename dim>
struct CudaMoments1stN
{
  void operator()(CudaMparticles& cmprts, struct cuda_mfields* cmres);
};

template <typename CudaMparticles, typename dim>
struct CudaMoments1stAll
{
  void operator()(CudaMparticles& cmprts, struct cuda_mfields* cmres);
};
