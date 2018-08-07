
#ifndef CUDA_TEST_HXX
#define CUDA_TEST_HXX

template<typename CMP>
struct TestBase
{
  using CudaMparticles = CMP;

  CudaMparticles* make_cmprts(Grid_t& grid)
  {
    auto cmprts = new CudaMparticles(grid);
    
    return cmprts;
  }

  CudaMparticles* make_cmprts(Grid_t& grid, const std::vector<cuda_mparticles_prt>& prts)
  {
    auto cmprts = new CudaMparticles(grid);

    uint n_prts_by_patch[1];
    n_prts_by_patch[0] = prts.size();

    cmprts->reserve_all(n_prts_by_patch);
    cmprts->inject_buf(prts.data(), n_prts_by_patch);
  
    return cmprts;
  }
  
  template<typename S>
  CudaMparticles* make_cmprts(Grid_t& grid, uint n_prts, S setter)
  {
    std::vector<cuda_mparticles_prt> prts;
    prts.reserve(n_prts);
  
    for (int i = 0; i < n_prts; i++) {
      cuda_mparticles_prt prt = setter(i);
      prts.push_back(prt);
    }

    return make_cmprts(grid, prts);
  }
};


#endif
