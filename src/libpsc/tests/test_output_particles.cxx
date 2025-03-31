
#include "gtest/gtest.h"

#include <dim.hxx>
#include <psc_particles_single.h>
#include <psc_particles_double.h>
#include "../libpsc/psc_output_particles/output_particles_ascii_impl.hxx"
#include "../libpsc/psc_output_particles/output_particles_hdf5_impl.hxx"

template <typename _Dim, typename _Mparticles, typename _OutputParticles>
struct Config
{
  using Dim = _Dim;
  using Mparticles = _Mparticles;
  using OutputParticles = _OutputParticles;
};

template <typename T>
struct OutputParticlesTest : ::testing::Test
{
  using Dim = typename T::Dim;

  const double L = 160;

  void make_psc(const Grid_t::Kinds& kinds)
  {
    Int3 gdims = {16, 16, 16};
    if (Dim::InvarX::value) {
      gdims[0] = 1;
      ibn[0] = 0;
    }
    if (Dim::InvarY::value) {
      gdims[1] = 1;
      ibn[1] = 0;
    }
    if (Dim::InvarZ::value) {
      gdims[2] = 1;
      ibn[2] = 0;
    }

    auto grid_domain =
      Grid_t::Domain{gdims, {L, L, L}, {0., 0., 0.}, {2, 1, 1}};
    auto grid_bc =
      psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                    {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                    {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                    {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 200;
    auto coeff = Grid_t::Normalization{norm_params};

    grid_ = new Grid_t{grid_domain, grid_bc, kinds, coeff, 1.};
  }

  const Grid_t& grid() const
  {
    assert(grid_);
    return *grid_;
  }

private:
  Grid_t* grid_;
  Int3 ibn = {2, 2, 2};
};

using OutputParticlesTestTypes = ::testing::Types<
  Config<dim_xyz, MparticlesSingle, OutputParticlesAscii>
#ifdef H5_HAVE_PARALLEL
  ,
  Config<dim_xyz, MparticlesSingle, OutputParticlesHdf5<ParticleSelectorAll>>,
  Config<dim_xyz, MparticlesDouble, OutputParticlesHdf5<ParticleSelectorAll>>
#endif
  >;

TYPED_TEST_SUITE(OutputParticlesTest, OutputParticlesTestTypes);

// ======================================================================
// Test1

TYPED_TEST(OutputParticlesTest, Test1)
{
  using Mparticles = typename TypeParam::Mparticles;
  using OutputParticles = typename TypeParam::OutputParticles;

  auto kinds = Grid_t::Kinds{{1., 100., "ion"}, {-1., 1., "electron"}};
  this->make_psc(kinds);
  const auto& grid = this->grid();

  // init particle
  auto n_prts_by_patch = std::vector<uint>{1};

  Mparticles mprts{grid};
  {
    auto injector = mprts.injector();
    for (int p = 0; p < grid.n_patches(); p++) {
      auto info = grid.mrc_domain().localPatchInfo(p);
      if (info.global_patch == 0) {
        injector[p]({{1., 0., 0.}, {}, 1., 0});
      }
      if (info.global_patch == 1) {
        injector[p]({{121., 0., 0.}, {}, 1., 1});
      }
    }
  }

  auto params = OutputParticlesParams{};
  params.every_step = 1;
  params.data_dir = ".";
  params.basename = "prt";

  auto outp = OutputParticles{grid, params};
  outp(mprts);
}

// ----------------------------------------------------------------------
// main

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
