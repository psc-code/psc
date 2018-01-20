
#include "grid.hxx"
#include "fields.hxx"
#include "cuda_mfields.h"
#include "cuda_mparticles.h"

#include "mrc_profile.h"

struct prof_globals prof_globals; // FIXME

int
prof_register(const char *name, float simd, int flops, int bytes)
{
  return 0;
}

class TestAccel
{
  enum { // FIXME, duplicated
    JXI, JYI, JZI,
    EX , EY , EZ ,
    HX , HY , HZ ,
    N_FIELDS,
  };

public:
  TestAccel()
  {
    init_grid();
    init_cmflds();
    init_cmprts();
  }

  ~TestAccel()
  {
    delete cmflds_;
    delete cmprts_;
  }

  void init_grid()
  {
    grid_.dt = 1.;

    Grid_t::Patch patch{};
    grid_.ldims = { 1, 1, 1 };
    grid_.dx = { L, L, L };
    
    patch.xb = { 0., 0., 0. };
    patch.xe = { L,  L,  L  };
    grid_.patches.push_back(patch);

    grid_.kinds.push_back(Grid_t::Kind(1.,  1., "test_species"));
  }
  
  void init_cmflds()
  {
    cmflds_ = new cuda_mfields(grid_, N_FIELDS, { 1, 1, 1 });
    fields_single_t flds = cmflds_->get_host_fields();
    Fields3d<fields_single_t> F(flds);
    F(EX, 0,0,0) = 1;
    F(EX, 0,1,0) = 1;
    F(EX, 0,0,1) = 1;
    F(EX, 0,1,1) = 1;
    
    F(EY, 0,0,0) = 2;
    F(EY, 0,0,1) = 2;
    F(EY, 1,0,0) = 2;
    F(EY, 1,0,1) = 2;
    
    F(EZ, 0,0,0) = 3;
    F(EZ, 1,0,0) = 3;
    F(EZ, 0,1,0) = 3;
    F(EZ, 1,1,0) = 3;

    cmflds_->copy_to_device(0, flds, 0, N_FIELDS);
    cmflds_->dump("accel.fld.json");
    flds.dtor();
  };

  void init_cmprts()
  {
    uint n_prts_by_patch[1] = { n_prts };
    
    cmprts_ = new cuda_mparticles(grid_, { 1, 1, 1 });
    cmprts_->reserve_all(n_prts_by_patch);
    
    std::vector<cuda_mparticles_prt> prts;
    prts.reserve(n_prts);
    
    for (int n = 0; n < n_prts; n++) {
      cuda_mparticles_prt prt = {};
      prts.push_back(prt);
      // inject_particle( sp,
      //                uniform( rng(0), 0, L ),
      //                uniform( rng(0), 0, L ),
      //                uniform( rng(0), 0, L ),
      //                0., 0., 0., 1., 0., 0 );
    }
    cmprts_->inject(prts.data(), n_prts_by_patch);
  }
  

  void run()
  {
  };

private:
  double L = 1e10;
  unsigned int n_prts = 13;

  Grid_t grid_;
  cuda_mfields* cmflds_;
  cuda_mparticles* cmprts_;
};

// ----------------------------------------------------------------------
// main

int
main(void)
{
  TestAccel test;
  test.run();
}
