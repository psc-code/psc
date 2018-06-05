
#pragma once

template<typename CudaMparticles>
struct CudaMoments1stNcRho
{
  void operator()(CudaMparticles& cmprts, struct cuda_mfields *cmres);

private:
  template<bool REORDER>
  void invoke(CudaMparticles& cmprts, struct cuda_mfields *cmres);
};

template<typename CudaMparticles>
struct CudaMoments1stNcN
{
  void operator()(CudaMparticles& cmprts, struct cuda_mfields *cmres);

private:
  template<bool REORDER>
  void invoke(CudaMparticles& cmprts, struct cuda_mfields *cmres);
};

