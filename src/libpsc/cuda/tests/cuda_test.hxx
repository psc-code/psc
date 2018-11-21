
#ifndef CUDA_TEST_HXX
#define CUDA_TEST_HXX

template<typename _CudaMparticles>
struct TestBase
{
  using CudaMparticles = _CudaMparticles;
  using particle_t = typename CudaMparticles::particle_t;

  CudaMparticles* make_cmprts(Grid_t& grid)
  {
    auto cmprts = new CudaMparticles(grid);
    
    return cmprts;
  }

  CudaMparticles* make_cmprts(Grid_t& grid, const std::vector<particle_t>& prts)
  {
    auto cmprts = new CudaMparticles(grid);

    std::vector<uint> n_prts_by_patch = { uint(prts.size()) };

    cmprts->reserve_all(n_prts_by_patch);
    cmprts->inject_buf(prts.data(), n_prts_by_patch.data());
  
    return cmprts;
  }
  
  template<typename S>
  CudaMparticles* make_cmprts(Grid_t& grid, uint n_prts, S setter)
  {
    std::vector<particle_t> prts;
    prts.reserve(n_prts);
  
    for (int i = 0; i < n_prts; i++) {
      auto prt = setter(i);
      prts.push_back(prt);
    }

    return make_cmprts(grid, prts);
  }
};


#endif
