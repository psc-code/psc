
#pragma once

template<typename BS>
struct CudaMoments1stNcRho
{
  void operator()(cuda_mparticles<BS>* cmprts, struct cuda_mfields *cmres);
};

template<typename BS>
struct CudaMoments1stNcN
{
  void operator()(cuda_mparticles<BS>* cmprts, struct cuda_mfields *cmres);
};

