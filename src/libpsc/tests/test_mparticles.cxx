
#include <gtest/gtest.h>

#include "test_common.hxx"

#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "../libpsc/vpic/mparticles_vpic.hxx"
#include "../libpsc/vpic/PscGridBase.h"
#include "../libpsc/vpic/PscParticleBc.h"
#include "../libpsc/vpic/PscParticlesBase.h"
#include "../libpsc/vpic/vpic_config.h"
#ifdef USE_CUDA
#include "../libpsc/cuda/mparticles_cuda.hxx"
#endif

#ifdef DO_VPIC
using VpicConfig = VpicConfigWrap;
#else
using VpicConfig = VpicConfigPsc;
#endif

using Grid = VpicConfig::Grid;
using MparticlesVpic = VpicConfig::Mparticles;

template<typename _Mparticles, typename _MakeGrid = MakeTestGrid1>
struct Config
{
  using Mparticles = _Mparticles;
  using MakeGrid = _MakeGrid;
};

using MparticlesTestTypes = ::testing::Types<Config<MparticlesSingle>
					    ,Config<MparticlesSingle, MakeTestGridYZ>
					    ,Config<MparticlesDouble>
					    ,Config<MparticlesVpic, MakeTestGridYZ1>
#ifdef USE_CUDA
					    ,Config<MparticlesCuda<BS144>, MakeTestGridYZ1>
					    ,Config<MparticlesCuda<BS144>, MakeTestGridYZ>
#endif
					     >;

TYPED_TEST_CASE(MparticlesTest, MparticlesTestTypes);

// ======================================================================
// MparticlesTest

template<typename T>
struct MparticlesTest : ::testing::Test
{
  using Mparticles = typename T::Mparticles;
  using Particle = typename Mparticles::Particle;
  using MakeGrid = typename T::MakeGrid;

  MparticlesTest()
    : grid_{MakeGrid{}()}
  {
    grid_.kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));
  }

  template<typename tag>
  Mparticles mk_mprts(tag dummy)
  {
    Mparticles mprts(grid_);
    mprts.define_species("test_species", 1., 1., 100, 10,
			 10, 0);
    return mprts;
  }

  Mparticles mk_mprts(MparticlesVpic* dummy)
  {
    // FIXME, vgrid_ is bad, and this kinda belongs to where grid_ is set up

    // Setup basic grid parameters
    auto& domain = grid_.domain;
    double dx[3], xl[3], xh[3];
    for (int d = 0; d < 3; d++) {
      dx[d] = domain.length[d] / domain.gdims[d];
      xl[d] = domain.corner[d];
      xh[d] = xl[d] + domain.length[d];
    }
    vgrid_.setup(dx, grid_.dt, 1., 1.);
  
    // Define the grid
    vgrid_.partition_periodic_box(xl, xh, domain.gdims, domain.np);

    Mparticles mprts(grid_, &vgrid_);
    mprts.define_species("test_species", 1., 1., 100, 10,
			 10, 0);
    return mprts;
  }

  Mparticles mk_mprts()
  {
    return mk_mprts(static_cast<Mparticles*>(nullptr));
  }

  template<typename _Mparticles>
  void inject_test_particles(_Mparticles& mprts, int n_prts)
  {
    auto inj = mprts.injector();
    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto injector = inj[p];
      auto& patch = mprts.grid().patches[p];
      for (int n = 0; n < n_prts; n++) {
	double nn = double(n) / n_prts;
	auto L = patch.xe - patch.xb;
	particle_inject prt = {};
	prt.x[0] = patch.xb[0] + nn * L[0];
	prt.x[1] = patch.xb[1] + nn * L[1];
	prt.x[2] = patch.xb[2] + nn * L[2];
	prt.kind = 0;
	prt.w = 1.;
	injector(prt);
      }
    }
  }

  const Grid_t& grid() { return grid_; }
  
private:
  Grid_t grid_;
  Grid vgrid_; // FIXME
};

// -----------------------------------------------------------------------
// Constructor

TYPED_TEST(MparticlesTest, Constructor)
{
  auto mprts = this->mk_mprts();
}

// ----------------------------------------------------------------------
// Inject

TYPED_TEST(MparticlesTest, Inject)
{
  const int n_prts = 4;

  auto mprts = this->mk_mprts();

  this->inject_test_particles(mprts, n_prts);
}

// ----------------------------------------------------------------------
// ConstAccessor

TYPED_TEST(MparticlesTest, ConstAccessor)
{
  using Base = MparticlesTest<TypeParam>;
  using Real3 = typename Base::Mparticles::Real3;

  auto mprts = this->mk_mprts();

  {
    auto inj = mprts.injector();
    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto& patch = mprts.grid().patches[p];
      auto injector = inj[p];
      auto x = .5 * (patch.xb + patch.xe);
      injector({{x[0], x[1], x[2]}, {1., 2., 3.}, double(p), 0});
      injector({{x[0], x[1], x[2]}, {1., 2., 3.}, double(p), 0});
    }
  }

  // check iterating of particles per patch
  {
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); ++p) {
      for (auto prt: accessor[p]) {
	EXPECT_EQ(prt.u(), (Real3{1., 2., 3.}));
	EXPECT_EQ(prt.w(), p);
      }
    }
  }

  // check indexed access
  {
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); ++p) {
      EXPECT_EQ(accessor[p][0].w(), p);
      EXPECT_EQ(accessor[p][1].w(), p);
    }
  }
}

// -----------------------------------------------------------------------
// Inject2

TYPED_TEST(MparticlesTest, Inject2)
{
  using Mparticles = typename TypeParam::Mparticles;
  const int n_prts = 1;
  
  auto mprts = this->mk_mprts();

  // FIXME, kinda the cheap way out
  if (mprts.n_patches() != 4) {
    return;
  }
  ASSERT_EQ(mprts.n_patches(), 4);

  int nn = 0;
  {
    auto inj = mprts.injector();
    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto injector = inj[p];
      auto& patch = mprts.grid().patches[p];
      for (int n = 0; n < n_prts; n++) {
	particle_inject prt = {};
	auto x = .5 * (patch.xb + patch.xe);
	int kind = 0;
	// use weight to store particle number for testing
	injector({{x[0], x[1], x[2]}, {}, double(nn), kind});
	nn++;
      }
    }
  }
  
  EXPECT_EQ(mprts.get_n_prts(), 4);
  
  auto n_prts_by_patch = mprts.get_size_all();
  EXPECT_EQ(n_prts_by_patch, std::vector<uint>({n_prts, n_prts, n_prts, n_prts}));

  nn = 0;
  auto accessor = mprts.accessor();
  for (int p = 0; p < mprts.n_patches(); ++p) {
    auto& patch = mprts.grid().patches[p];
    for (auto prt: accessor[p]) {
      // xm is patch-relative position
      auto xm = .5 * (patch.xe - patch.xb);
      EXPECT_EQ(prt.x()[0], xm[0]);
      EXPECT_EQ(prt.x()[1], xm[1]);
      EXPECT_EQ(prt.x()[2], xm[2]);
      EXPECT_EQ(prt.qni_wni(), nn);
      EXPECT_EQ(prt.w(), nn);
      EXPECT_EQ(prt.kind(), 0);

      auto x = .5 * (patch.xb + patch.xe);
      EXPECT_EQ(prt.position()[0], x[0]);
      EXPECT_EQ(prt.position()[1], x[1]);
      EXPECT_EQ(prt.position()[2], x[2]);
      nn++;
    }
  }
}

// ----------------------------------------------------------------------
// Injector

TYPED_TEST(MparticlesTest, Injector)
{
  const int n_prts = 4;

  auto mprts = this->mk_mprts();

  this->inject_test_particles(mprts, n_prts);

  auto accessor = mprts.accessor();
  for (int p = 0; p < mprts.n_patches(); ++p) {
    auto prts = accessor[p];
    auto& patch = mprts.grid().patches[p];
    EXPECT_EQ(prts.size(), n_prts);
    int n = 0;
    for (auto prt : prts) {
      double nn = double(n) / n_prts;
      auto L = patch.xe - patch.xb;
      auto x = prt.position();
      EXPECT_EQ(x[0], patch.xb[0] + nn * L[0]);
      EXPECT_EQ(x[1], patch.xb[1] + nn * L[1]);
      EXPECT_EQ(x[2], patch.xb[2] + nn * L[2]);
      EXPECT_EQ(prt.w(), 1.);
      EXPECT_EQ(prt.kind(), 0);
      ++n;
    }
  }
}

// ----------------------------------------------------------------------
// Conversion to MparticlesSingle

template<typename Mparticles, typename MakeTestGrid>
struct TestConversionToMparticlesSingle
{
  void operator()(Mparticles& mprts) {}
};

template<typename Mparticles>
struct TestConversionToMparticlesSingle<Mparticles, MakeTestGridYZ1>
{
  using Double3 = Vec3<double>;

  void operator()(Mparticles& mprts)
  {
    auto&& mprts_single = mprts.template get_as<MparticlesSingle>();
    EXPECT_EQ(mprts_single.get_n_prts(), 2);
    
    {
      auto accessor = mprts_single.accessor();
      EXPECT_EQ(accessor[0][0].position(), (Double3{0., -40., -80.}));
      EXPECT_EQ(accessor[0][1].position(), (Double3{5., 0., 0.}));
    }
  }
};

TYPED_TEST(MparticlesTest, ConversionToMparticlesSingle)
{
  auto mprts = this->mk_mprts();
  this->inject_test_particles(mprts, 2);

#if 0
  {
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      for (auto prt : accessor[p]) {
	mprintf("xyz %g %g %g\n", prt.position()[0], prt.position()[1], prt.position()[2]);
      }
    }
  }
#endif
  auto test = TestConversionToMparticlesSingle<typename TypeParam::Mparticles, typename TypeParam::MakeGrid>{};
  test(mprts);
}

// ----------------------------------------------------------------------
// Conversion from MparticlesSingle

template<typename Mparticles, typename MakeTestGrid>
struct TestConversionFromMparticlesSingle
{
  void operator()(MparticlesSingle& mprts_single) {}
};

template<typename Mparticles>
struct TestConversionFromMparticlesSingle<Mparticles, MakeTestGridYZ1>
{
  using Double3 = Vec3<double>;

  void operator()(MparticlesSingle& mprts_single)
  {
    auto mprts = mprts_single.template get_as<Mparticles>();
    EXPECT_EQ(mprts.get_n_prts(), 2);
    
    {
      auto accessor = mprts.accessor();
      EXPECT_EQ(accessor[0][0].position(), (Double3{0., -40., -80.}));
      EXPECT_EQ(accessor[0][1].position(), (Double3{5., 0., 0.}));
    }
  }
};

// FIXME, MparticlesVpic can't be converted the usual way...
template<>
struct TestConversionFromMparticlesSingle<MparticlesVpic, MakeTestGridYZ1>
{
  using Double3 = Vec3<double>;

  void operator()(MparticlesSingle& mprts_single) {}
};

TYPED_TEST(MparticlesTest, ConversionFromMparticlesSingle)
{
  MparticlesSingle mprts(this->grid());
  this->inject_test_particles(mprts, 2);

#if 0
  {
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      for (auto prt : accessor[p]) {
	mprintf("xyz %g %g %g\n", prt.position()[0], prt.position()[1], prt.position()[2]);
      }
    }
  }
#endif
  auto test = TestConversionFromMparticlesSingle<typename TypeParam::Mparticles, typename TypeParam::MakeGrid>{};
  test(mprts);
}


int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  //#ifdef USE_VPIC FIXME
  MPI_Comm_dup(MPI_COMM_WORLD, &psc_comm_world);
  MPI_Comm_rank(psc_comm_world, &psc_world_rank);
  MPI_Comm_size(psc_comm_world, &psc_world_size);
  //#endif

  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();

  MPI_Finalize();
  return rc;
}
