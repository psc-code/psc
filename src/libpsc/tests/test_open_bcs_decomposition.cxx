#include "gtest/gtest.h"

#include <psc.hxx>

#include "../libpsc/psc_bnd_fields/radiating.hxx"
#include "../psc_config.hxx"

// ======================================================================
// Results must not depend on how the domain is decomposed into patches. These
// tests run the same setup with 1 and 2 patches along z (periodic), with open
// (radiating) boundaries in y, and compare the resulting fields.

using Dim = dim_yz;
using PscConfig = PscConfig1vbecDouble<Dim>;

using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;
using real_t = MfieldsState::real_t;
using Real3 = Vec3<real_t>;

using RadiatingBC =
  psc::bnd::field::Radiating<Dim, MfieldsState,
                             psc::bnd::field::ConstantPulse<real_t>>;

namespace
{

const Int3 gdims = {1, 8, 8};
const int n_steps = 10;

Grid_t* setupGrid(const PscParams& psc_params, int n_patches_z)
{
  auto domain = Grid_t::Domain{gdims,
                               {1.0, double(gdims[1]), double(gdims[2])},
                               {0.0, 0.0, 0.0},
                               {1, 1, n_patches_z}};

  auto bc = psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_OPEN, BND_FLD_PERIODIC},
                          {BND_FLD_PERIODIC, BND_FLD_OPEN, BND_FLD_PERIODIC},
                          {BND_PRT_PERIODIC, BND_PRT_OPEN, BND_PRT_PERIODIC},
                          {BND_PRT_PERIODIC, BND_PRT_OPEN, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1.0, 1.0, "e"};
  kinds[KIND_ION] = {1.0, 100.0, "i"};

  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 1;
  Grid_t::Normalization norm{norm_params};

  double dt = psc_params.cfl * courant_length(domain);

  Int3 ibn = {0, 2, 2};
  return new Grid_t{domain, bc, kinds, norm, dt, -1, ibn};
}

int find_patch(const Grid_t& grid, const Double3& x)
{
  for (int p = 0; p < grid.n_patches(); p++) {
    const auto& patch = grid.patches[p];
    if (x[1] >= patch.xb[1] && x[1] < patch.xe[1] && x[2] >= patch.xb[2] &&
        x[2] < patch.xe[2]) {
      return p;
    }
  }
  return -1;
}

// global field values, indexed [m][j][k] over the interior cells
using GlobalFields = std::vector<std::vector<std::vector<double>>>;

// Runs a neutral electron-ion pair, starting at `x` with the electron moving
// with momentum `u_electron`, and returns the global fields after n_steps.
GlobalFields run(int n_patches_z, Double3 x, Double3 u_electron)
{
  PscParams psc_params;
  psc_params.nmax = n_steps;
  psc_params.cfl = .75;

  auto grid_ptr = setupGrid(psc_params, n_patches_z);
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};

  ChecksParams checks_params{};
  Checks checks{grid, MPI_COMM_WORLD, checks_params};
  Balance balance{.1};
  Collision collision{grid, 0, 0.1};

  auto psc = makePscIntegrator<PscConfig>(psc_params, grid, mflds, mprts,
                                          balance, collision, checks);

  using psc::Axis;
  using psc::bnd::LoHi;
  using ConstantPulse = psc::bnd::field::ConstantPulse<real_t>;
  Real3 zero = {0.0, 0.0, 0.0};
  psc.add_field_bc(
    new RadiatingBC(ConstantPulse{zero, zero}, Axis::Y, LoHi::Lo));
  psc.add_field_bc(
    new RadiatingBC(ConstantPulse{zero, zero}, Axis::Y, LoHi::Hi));

  {
    int p = find_patch(grid, x);
    EXPECT_GE(p, 0);
    auto injector = mprts.injector();
    auto inj = injector[p];
    inj({x, u_electron, 1, KIND_ELECTRON});
    inj({x, {0.0, 0.0, 0.0}, 1, KIND_ION});
  }

  psc.pre_first_step();
  for (; grid.timestep_ < psc_params.nmax;) {
    psc.step();
  }

  GlobalFields global(NR_FIELDS, std::vector<std::vector<double>>(
                                   gdims[1], std::vector<double>(gdims[2])));
  for (int p = 0; p < grid.n_patches(); p++) {
    auto F = make_Fields3d<dim_xyz>(mflds[p]);
    Int3 off = grid.patches[p].off;
    for (int m = 0; m < NR_FIELDS; m++) {
      for (int k = 0; k < grid.ldims[2]; k++) {
        for (int j = 0; j < grid.ldims[1]; j++) {
          global[m][j + off[1]][k + off[2]] = F(m, 0, j, k);
        }
      }
    }
  }
  return global;
}

void expect_same_fields(const GlobalFields& a, const GlobalFields& b)
{
  const char* names[] = {"jx", "jy", "jz", "ex", "ey", "ez", "hx", "hy", "hz"};
  for (int m = 0; m < NR_FIELDS; m++) {
    double max_abs = 0.0;
    for (int j = 0; j < gdims[1]; j++) {
      for (int k = 0; k < gdims[2]; k++) {
        max_abs = std::max(max_abs, std::abs(a[m][j][k]));
      }
    }
    for (int j = 0; j < gdims[1]; j++) {
      for (int k = 0; k < gdims[2]; k++) {
        EXPECT_NEAR(a[m][j][k], b[m][j][k], 1e-10 * max_abs)
          << names[m] << " at j=" << j << " k=" << k;
      }
    }
  }
}

} // namespace

// Particle in the first cell layer in y, near the patch boundary in z. This
// is the control case.
TEST(OpenBcsDecompositionTest, LowerBoundary)
{
  Double3 x = {0.5, 0.5, 3.5};
  Double3 u = {0.5, 0.0, 0.5};
  auto one_patch = run(1, x, u);
  auto two_patches = run(2, x, u);
  expect_same_fields(one_patch, two_patches);
}

// Same, but mirrored to the last cell layer in y.
TEST(OpenBcsDecompositionTest, UpperBoundary)
{
  Double3 x = {0.5, gdims[1] - 0.5, 3.5};
  Double3 u = {0.5, 0.0, 0.5};
  auto one_patch = run(1, x, u);
  auto two_patches = run(2, x, u);
  expect_same_fields(one_patch, two_patches);
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  psc_init(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  psc_finalize();
  return rc;
}
