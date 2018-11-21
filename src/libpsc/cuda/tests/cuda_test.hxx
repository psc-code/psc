
#ifndef CUDA_TEST_HXX
#define CUDA_TEST_HXX

template<typename _CudaMparticles>
struct TestBase
{
  using CudaMparticles = _CudaMparticles;
  using particle_t = typename CudaMparticles::particle_t;
};


#endif
