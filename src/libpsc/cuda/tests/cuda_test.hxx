
#ifndef CUDA_TEST_HXX
#define CUDA_TEST_HXX

template <typename _CudaMparticles>
struct TestBase
{
  using CudaMparticles = _CudaMparticles;
  using Particle = typename CudaMparticles::Particle;
};

#endif
