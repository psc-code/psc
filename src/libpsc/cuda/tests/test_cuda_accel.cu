
#include "grid.hxx"
#include "fields.hxx"
#include "cuda_mfields.h"
#include "cuda_mparticles.h"

#include "../vpic/PscRng.h"

#include "gtest/gtest.h"

// Rng hackiness

using Rng = PscRng;
using RngPool = PscRngPool<Rng>;

// enum hackiness

enum IP { // FIXME, dup
  IP_STD, // standard interpolation
  IP_EC,  // energy-conserving interpolation
};

enum DEPOSIT { // FIXME, dup
  DEPOSIT_VB_2D,
  DEPOSIT_VB_3D,
};

enum CURRMEM { // FIXME, dup
  CURRMEM_SHARED,
  CURRMEM_GLOBAL,
};

enum { // FIXME, duplicated
  JXI, JYI, JZI,
  EX , EY , EZ ,
  HX , HY , HZ ,
  N_FIELDS,
};

// profile hackiness

#include "mrc_profile.h"

struct prof_globals prof_globals; // FIXME

int
prof_register(const char *name, float simd, int flops, int bytes)
{
  return 0;
}

// ======================================================================
// class PushMprtsTest

struct PushMprtsTest : ::testing::Test
{
  Grid_t* grid_;
  cuda_mparticles* cmprts_;
  cuda_mfields* cmflds_;

  RngPool rngpool;
  
  const double L = 1e10;
  const Int3 bs_ = { 1, 1, 1 };

  void SetUp()
  {
    grid_ = new Grid_t({ 1, 1, 1 }, { L, L, L });
  }

  template<typename S>
  cuda_mfields *make_cmflds(S set)
  {
    cuda_mfields *cmflds = new cuda_mfields(*grid_, N_FIELDS, { 0, 2, 2 });

    fields_single_t flds = cmflds->get_host_fields();
    Fields3d<fields_single_t> F(flds);

    F(EX, 0,0,0) = set(EX);
    F(EX, 0,1,0) = set(EX);
    F(EX, 0,0,1) = set(EX);
    F(EX, 0,1,1) = set(EX);
    
    F(EY, 0,0,0) = set(EY);
    F(EY, 0,0,1) = set(EY);
    //    F(EY, 1,0,0) = set(EY);
    //    F(EY, 1,0,1) = set(EY);
    
    F(EZ, 0,0,0) = set(EZ);
    //    F(EZ, 1,0,0) = set(EZ);
    F(EZ, 0,1,0) = set(EZ);
    //    F(EZ, 1,1,0) = set(EZ);

    F(HX, 0,0,0) = set(HX);
    F(HX, 1,0,0) = set(HX);

    F(HY, 0,0,0) = set(HY);
    F(HY, 0,1,0) = set(HY);

    F(HZ, 0,0,0) = set(HZ);
    F(HZ, 0,0,1) = set(HZ);

    cmflds->copy_to_device(0, flds, 0, N_FIELDS);
    cmflds->dump("accel.fld.json");
    flds.dtor();
  
    return cmflds;
  }

  cuda_mparticles* make_cmprts(uint n_prts)
  {
    cuda_mparticles* cmprts = new cuda_mparticles(*grid_, bs_);

    uint n_prts_by_patch[1] = { n_prts };
    cmprts->reserve_all(n_prts_by_patch);
  
    return cmprts;
  }
};

// ======================================================================
// Accel test

TEST_F(PushMprtsTest, Accel)
{
  const int n_prts = 131;
  const int n_steps = 10;
  const cuda_mparticles::real_t eps = 1e-5;

  // init fields
  cmflds_ = make_cmflds([&] (int m) -> cuda_mfields::real_t {
      switch(m) {
      case EX: return 1.;
      case EY: return 2.;
      case EZ: return 3.;
      default: return 0.;
      }
    });

  // init particles
  Rng *rng = rngpool[0];

  grid_->kinds.push_back(Grid_t::Kind(1., 1., "test_species"));
  cmprts_ = make_cmprts(n_prts);
  
  std::vector<cuda_mparticles_prt> prts;
  prts.reserve(n_prts);
  
  for (int i = 0; i < n_prts; i++) {
    cuda_mparticles_prt prt = {};
    prt.xi[0] = rng->uniform(0, L);
    prt.xi[1] = rng->uniform(0, L);
    prt.xi[2] = rng->uniform(0, L);
    prt.qni_wni = 1.;
    
    prts.push_back(prt);
  }
  uint n_prts_by_patch[1] = { n_prts };
  cmprts_->inject(prts.data(), n_prts_by_patch);
  //cmprts_->dump();
  
  int n_failed = 0;
  for (int n = 0; n < n_steps; n++) {
    cuda_push_mprts_yz(cmprts_, cmflds_, bs_, IP_EC, DEPOSIT_VB_3D, CURRMEM_GLOBAL);
    cmprts_->get_particles(0, [&] (int i, const cuda_mparticles_prt &prt) {
	if (std::abs(prt.pxi[0] - 1*(n+1)) > eps ||
	    std::abs(prt.pxi[1] - 2*(n+1)) > eps ||
	    std::abs(prt.pxi[2] - 3*(n+1)) > eps) {
	  printf("FAIL: n %d i %d px %g %g %g // exp %g %g %g\n", n, i,
		 prt.pxi[0], prt.pxi[1], prt.pxi[2],
		 1.*(n+1), 2.*(n+1), 3.*(n+1));
	  n_failed++;
	}
      });
    
    //cmprts_->dump();
  }
  EXPECT_EQ(n_failed, 0);
}

// ======================================================================
// Cyclo test

TEST_F(PushMprtsTest, Cyclo)
{
  const int n_prts = 131;
  const int n_steps = 64;
  const cuda_mparticles::real_t eps = 1e-2;

  // init fields
  cmflds_ = make_cmflds([&] (int m) -> cuda_mfields::real_t {
      switch(m) {
      case HZ: return 2. * M_PI / n_steps;
      default: return 0.;
      }
    });

  // init particles
  Rng *rng = rngpool[0];

  grid_->kinds.push_back(Grid_t::Kind(2., 1., "test_species"));
  cmprts_ = make_cmprts(n_prts);
  
  std::vector<cuda_mparticles_prt> prts;
  prts.reserve(n_prts);
  
  for (int i = 0; i < n_prts; i++) {
    cuda_mparticles_prt prt = {};
    prt.xi[0] = rng->uniform(0, L);
    prt.xi[1] = rng->uniform(0, L);
    prt.xi[2] = rng->uniform(0, L);
    prt.pxi[0] = 1.; // gamma = 2
    prt.pxi[1] = 1.;
    prt.pxi[2] = 1.;
    prt.qni_wni = rng->uniform(0, 1.);;
    
    prts.push_back(prt);
  }
  uint n_prts_by_patch[1] = { n_prts };
  cmprts_->inject(prts.data(), n_prts_by_patch);
  //cmprts_->dump();
  
  int n_failed = 0;
  
  for (int n = 0; n < n_steps; n++) {
    cuda_push_mprts_yz(cmprts_, cmflds_, bs_, IP_EC, DEPOSIT_VB_3D, CURRMEM_GLOBAL);
    double ux = (cos(2*M_PI*(0.125*n_steps-(n+1))/(double)n_steps) /
		 cos(2*M_PI*(0.125*n_steps)      /(double)n_steps));
    double uy = (sin(2*M_PI*(0.125*n_steps-(n+1))/(double)n_steps) /
		 sin(2*M_PI*(0.125*n_steps)      /(double)n_steps));
    double uz = 1.;
    cmprts_->get_particles(0, [&] (int i, const cuda_mparticles_prt &prt) {
	if (std::abs(prt.pxi[0] - ux) > eps ||
	    std::abs(prt.pxi[1] - uy) > eps ||
	    std::abs(prt.pxi[2] - uz) > eps) {
	  printf("FAIL: n %d i %d px %g %g %g // exp %g %g %g\n", n, i,
		 prt.pxi[0], prt.pxi[1], prt.pxi[2], ux, uy, uz);
	  n_failed++;
	}
      });
    
    //cmprts_->dump();
  }
  EXPECT_EQ(n_failed, 0);
}

