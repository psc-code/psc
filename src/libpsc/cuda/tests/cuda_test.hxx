
#ifndef CUDA_TEST_HXX
#define CUDA_TEST_HXX

struct TestBase
{
  cuda_mparticles* make_cmprts(Grid_t& grid)
  {
    auto cmprts = new cuda_mparticles(grid, bs_);
    
    return cmprts;
  }

  template<typename S>
  cuda_mparticles* make_cmprts(Grid_t& grid, uint n_prts, S setter)
  {
    auto cmprts = new cuda_mparticles(grid, bs_);

    uint n_prts_by_patch[1] = { n_prts };
    cmprts->reserve_all(n_prts_by_patch);
  
    std::vector<cuda_mparticles_prt> prts;
    prts.reserve(n_prts);
  
    for (int i = 0; i < n_prts; i++) {
      cuda_mparticles_prt prt = setter(i);
      prts.push_back(prt);
    }

    cmprts->inject(prts.data(), n_prts_by_patch);
  
    return cmprts;
  }

private:
  Int3 bs_ = { 1, 1, 1 };
};


#endif
