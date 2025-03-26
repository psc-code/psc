
#include <gtest/gtest.h>

#include "test_common.hxx"

#include "psc_particles_double.h"
#include "psc_particles_single.h"
#include "particle_with_id.h"
#include "setup_particles.hxx"
#ifdef USE_CUDA
#include "../libpsc/cuda/mparticles_cuda.hxx"
#include "../libpsc/cuda/mparticles_cuda.inl"
#endif
#include "particles_simple.inl"
#include <kg/io.h>

template <typename _Mparticles, typename _MakeGrid = MakeTestGrid1>
struct Config
{
  using Mparticles = _Mparticles;
  using MakeGrid = _MakeGrid;
};

using MparticlesInflowTestTypes =
  ::testing::Types<Config<MparticlesSingle>,
                   Config<MparticlesSingle, MakeTestGridYZ>
#ifdef USE_CUDA
                   ,
                   Config<MparticlesCuda<BS144>, MakeTestGridYZ1>,
                   Config<MparticlesCuda<BS144>, MakeTestGridYZ>
#endif
                   >;

TYPED_TEST_SUITE(MparticlesInflowTest, MparticlesInflowTestTypes);

// ======================================================================
// MparticlesInflowTest

template <typename T>
struct MparticlesInflowTest : ::testing::Test
{
  using Mparticles = typename T::Mparticles;
  using Particle = typename Mparticles::Particle;
  using MakeGrid = typename T::MakeGrid;

  MparticlesInflowTest() : grid_{MakeGrid{}()}
  {
    grid_.kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));
  }

  template <typename tag>
  Mparticles mk_mprts(tag dummy)
  {
    Mparticles mprts(grid_);
    mprts.define_species("test_species", 1., 1., 100, 10, 10, 0);
    return mprts;
  }

  Mparticles mk_mprts() { return mk_mprts(static_cast<Mparticles*>(nullptr)); }

  template <typename _Mparticles>
  void inject_test_particles(_Mparticles& mprts, int n_prts)
  {
    auto inj = mprts.injector();
    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto injector = inj[p];
      auto& patch = mprts.grid().patches[p];
      for (int n = 0; n < n_prts; n++) {
        double nn = double(n) / n_prts;
        auto L = patch.xe - patch.xb;
        psc::particle::Inject prt = {{patch.xb[0] + nn * L[0],
                                      patch.xb[1] + nn * L[1],
                                      patch.xb[2] + nn * L[2]},
                                     {},
                                     1.,
                                     0};
        injector(prt);
      }
    }
  }

  const Grid_t& grid() { return grid_; }

private:
  Grid_t grid_;
};

// ----------------------------------------------------------------------
// Inject

TYPED_TEST(MparticlesInflowTest, Inject)
{
  const int n_prts = 4;

  auto mprts = this->mk_mprts();

  this->inject_test_particles(mprts, n_prts);
}

// -----------------------------------------------------------------------
// Inject2

TYPED_TEST(MparticlesInflowTest, Inject2)
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
        psc::particle::Inject prt{
          .5 * (patch.xb + patch.xe), {}, double(nn), 0};
        // use weight to store particle number for testing
        injector(prt);
        nn++;
      }
    }
  }

  EXPECT_EQ(mprts.size(), 4);

  auto n_prts_by_patch = mprts.sizeByPatch();
  EXPECT_EQ(n_prts_by_patch,
            std::vector<uint>({n_prts, n_prts, n_prts, n_prts}));

  nn = 0;
  auto accessor = mprts.accessor();
  for (int p = 0; p < mprts.n_patches(); ++p) {
    auto& patch = mprts.grid().patches[p];
    for (auto prt : accessor[p]) {
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

TYPED_TEST(MparticlesInflowTest, Injector)
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

TEST(TestSetupParticlesInflow, Simple)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, .1};
  Mparticles mprts{grid};

  SetupParticles<Mparticles> setup_particles(grid);

  Double3 pos = {5., 5., 5.};
  Double3 u = Double3{0.0, 0.0, 0.0};
  int tag = 0;
  psc_particle_np np = {0, 1.0, [u]() { return u; }, tag};
  double wni = 1.0;
  auto prt = setup_particles.setupParticle(np, pos, wni);

  EXPECT_EQ(prt.x, pos);
  EXPECT_EQ(prt.u, u);
  EXPECT_EQ(prt.w, wni);
  EXPECT_EQ(prt.kind, 0);
  EXPECT_EQ(prt.tag, tag);
}

TEST(TestSetupParticlesInflow, Maxwellian)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, .1};
  Mparticles mprts{grid};

  SetupParticles<Mparticles> setup_particles(grid);

  Double3 pos = {5., 5., 5.};
  Double3 u = {0.0, 0.5, 0.0};
  Double3 T = {0.0, 0.0, 0.0};
  int tag = 0;

  psc_particle_npt npt = {0, 1.0, u, T, tag};
  auto p = setup_particles.createMaxwellian(npt);
  psc_particle_np np = {0, 1.0, p, tag};
  double wni = 1.0;
  auto prt = setup_particles.setupParticle(np, pos, wni);

  EXPECT_EQ(prt.x, pos);
  EXPECT_EQ(prt.u, u);
  EXPECT_EQ(prt.w, wni);
  EXPECT_EQ(prt.kind, 0);
  EXPECT_EQ(prt.tag, tag);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();

  MPI_Finalize();
  return rc;
}
