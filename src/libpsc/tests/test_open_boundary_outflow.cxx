#include <gtest/gtest.h>

#include "test_common.hxx"

#include "psc.hxx"
#include "output_fields.hxx"
#include "../psc_config.hxx"

// Particles leaving through open boundaries, with no injectors: their virtual
// charge must end up beyond the boundary nodes, so Gauss's law holds on every
// node of the domain.

using Dim = dim_yz;
using PscConfig = PscConfig1vbecDouble<Dim>;

using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;

Grid_t* setupGrid()
{
  auto domain = Grid_t::Domain{{1, 8, 8},  // n grid points
                               {1, 8, 8},  // physical lengths
                               {0, 0, 0},  // location of lower corner
                               {1, 1, 1}}; // n patches

  auto bc = psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_OPEN, BND_FLD_OPEN},
                          {BND_FLD_PERIODIC, BND_FLD_OPEN, BND_FLD_OPEN},
                          {BND_PRT_PERIODIC, BND_PRT_OPEN, BND_PRT_OPEN},
                          {BND_PRT_PERIODIC, BND_PRT_OPEN, BND_PRT_OPEN}};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1.0, 1.0, "e"};
  kinds[KIND_ION] = {1.0, 1.0, "i"};

  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 1;

  double dt = 1;
  Grid_t::Normalization norm{norm_params};

  int n_ghosts = 2;
  Int3 ibn = n_ghosts * Dim::get_noninvariant_mask();

  return new Grid_t{domain, bc, kinds, norm, dt, -1, ibn};
}

// max |div E - rho| over the domain's nodes, lower boundary nodes included
template <typename Mparticles, typename MfieldsState>
double gauss_err(Mparticles& mprts, MfieldsState& mflds)
{
  const auto& grid = mprts.grid();
  using ItemRho = Checks::item_rho_type;
  auto rho = psc::mflds::interior(grid, ItemRho{grid}(mprts));
  auto dive = psc::mflds::interior(grid, Item_dive<MfieldsState>{}(mflds));
  return gt::norm_linf(gt::eval(dive - rho));
}

// Places a neutral pair at `x`, gives each particle its momentum, and steps
// once.
template <typename Check>
void run_pair(Vec3<double> x, Vec3<double> u_e, Vec3<double> u_i, Check check)
{
  PscParams psc_params;
  psc_params.nmax = 1;

  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};
  Checks checks{grid, MPI_COMM_WORLD, ChecksParams{}};
  Balance balance{.1};
  Collision collision{grid, 0, 0.1};

  auto psc = makePscIntegrator<PscConfig>(psc_params, grid, mflds, mprts,
                                          balance, collision, checks);
  {
    auto injector = mprts.injector();
    auto inj = injector[0];
    inj({x, u_e, 1, KIND_ELECTRON});
    inj({x, u_i, 1, KIND_ION});
  }

  psc.pre_first_step();
  psc.step();
  check(grid, mprts, mflds);
}

TEST(OpenBoundaryOutflowTest, Lower)
{
  // the electron stops at y = -.594 and is deleted
  run_pair({.5, .3, 4.5}, {0, -2, 0}, {0, 0, 0},
           [](const Grid_t& grid, Mparticles& mprts, MfieldsState& mflds) {
             ASSERT_EQ(mprts.size(), 1);
             EXPECT_LT(gauss_err(mprts, mflds), 1e-10);
           });
}

TEST(OpenBoundaryOutflowTest, Upper)
{
  // both leave, to different depths in the ghost cell. No real charge is
  // left, so no virtual charge may remain on the boundary node either: E_y
  // vanishes on both sides of it (the upper boundary nodes aren't interior)
  run_pair({.5, 7.7, 4.5}, {0, 2, 0}, {0, 1, 0},
           [](const Grid_t& grid, Mparticles& mprts, MfieldsState& mflds) {
             ASSERT_EQ(mprts.size(), 0);
             auto F = make_Fields3d<Dim>(mflds[0]);
             int ny = grid.ldims[1];
             for (int k = 0; k < grid.ldims[2]; k++) {
               EXPECT_NEAR(F(EY, 0, ny - 1, k), 0., 1e-10);
               EXPECT_NEAR(F(EY, 0, ny, k), 0., 1e-10);
             }
           });
}

TEST(OpenBoundaryOutflowTest, Corner)
{
  // the electron stops at (y, z) = (-.367, -.367), beyond two boundaries
  run_pair({.5, .3, .3}, {0, -2, -2}, {0, 0, 0},
           [](const Grid_t& grid, Mparticles& mprts, MfieldsState& mflds) {
             ASSERT_EQ(mprts.size(), 1);
             EXPECT_LT(gauss_err(mprts, mflds), 1e-10);
           });
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
