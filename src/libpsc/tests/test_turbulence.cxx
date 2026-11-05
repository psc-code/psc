#include <gtest/gtest.h>

#include "psc_fields_c.h"
#include "../psc_turbulence/alfven_turbulence.hxx"
#include "../psc_turbulence/noise_turbulence.hxx"

#include <memory>

using namespace psc::turbulence;

using Dim = dim_yz;
using MfieldsState = MfieldsStateDouble;

const TurbulenceParams default_params{0.01, 1.0, 5};

std::unique_ptr<Grid_t> setupGrid(Int3 n_patches = {1, 1, 1})
{
  double dx = 0.1;
  Int3 gdims{1, 32, 16};
  auto domain =
    Grid_t::Domain{gdims, Double3(gdims) * dx, {0, 0, 0}, n_patches};

  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds{};
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  Grid_t::Normalization norm{norm_params};

  int n_ghosts = 2;
  Int3 ibn = n_ghosts * Dim::get_noninvariant_mask();

  return std::unique_ptr<Grid_t>(
    new Grid_t{domain, bc, kinds, norm, 0.05, -1, ibn});
}

// H over the whole domain, indexed [m][global cell], gathered from all patches
std::vector<std::vector<double>> gather_h(MfieldsState& mflds)
{
  const auto& grid = mflds.grid();
  Int3 gdims = grid.domain.gdims;
  std::vector<std::vector<double>> h(3, std::vector<double>(gdims.prod()));
  for (int p = 0; p < mflds.n_patches(); p++) {
    auto f = make_Fields3d<dim_xyz>(mflds[p]);
    Int3 off = grid.patches[p].off;
    grid.Foreach_3d(0, 0, [&](int jx, int jy, int jz) {
      Int3 j_global = Int3{jx, jy, jz} + off;
      int idx = flatten_index(j_global, gdims);
      for (int m = 0; m < 3; m++) {
        h[m][idx] = f(HX + m, jx, jy, jz);
      }
    });
  }
  return h;
}

double max_abs_h(MfieldsState& mflds)
{
  double max = 0.0;
  for (auto& h_m : gather_h(mflds)) {
    for (double v : h_m) {
      max = std::max(max, std::abs(v));
    }
  }
  return max;
}

template <typename T>
struct TurbulenceTest : ::testing::Test
{
  using Generator = T;

  void generate(MfieldsState& mflds, TurbulenceParams params = default_params)
  {
    Generator{params}.generate(mflds);
  }
};

using Generators = ::testing::Types<AlfvenTurbulence<MfieldsState>,
                                    NoiseTurbulence<MfieldsState>>;
TYPED_TEST_SUITE(TurbulenceTest, Generators);

TYPED_TEST(TurbulenceTest, MeanB2MatchesTarget)
{
  auto grid = setupGrid();
  MfieldsState mflds{*grid};
  this->generate(mflds);

  EXPECT_NEAR(get_mean_b2(mflds), default_params.db2, 1e-12);
}

TYPED_TEST(TurbulenceTest, IsPerturbationInEveryComponent)
{
  auto grid = setupGrid();
  MfieldsState mflds{*grid};
  this->generate(mflds);

  // in yz, HX comes from AY, AZ and HY, HZ from AX, so all should be nonzero
  auto h = gather_h(mflds);
  for (int m = 0; m < 3; m++) {
    double max = 0.0;
    for (double v : h[m]) {
      max = std::max(max, std::abs(v));
    }
    EXPECT_GT(max, 1e-3 * std::sqrt(default_params.db2)) << "m = " << m;
  }
}

TYPED_TEST(TurbulenceTest, DivergenceFree)
{
  auto grid = setupGrid();
  MfieldsState mflds{*grid};
  this->generate(mflds);

  auto dx = grid->domain.dx;
  auto f = make_Fields3d<dim_xyz>(mflds[0]);
  double max_div = 0.0;
  // x is invariant
  grid->Foreach_3d(0, 0, [&](int jx, int jy, int jz) {
    double div = (f(HY, jx, jy + 1, jz) - f(HY, jx, jy, jz)) / dx[1] +
                 (f(HZ, jx, jy, jz + 1) - f(HZ, jx, jy, jz)) / dx[2];
    max_div = std::max(max_div, std::abs(div));
  });

  double b_over_dx = max_abs_h(mflds) / dx[1];
  EXPECT_LT(max_div, 1e-12 * b_over_dx);
}

TYPED_TEST(TurbulenceTest, Periodic)
{
  auto grid = setupGrid();
  MfieldsState mflds{*grid};
  this->generate(mflds);

  // the ghosts H is set in must match the interior a period away
  Int3 ldims = grid->ldims;
  auto f = make_Fields3d<dim_xyz>(mflds[0]);
  double max_diff = 0.0;
  for (int m = HX; m <= HZ; m++) {
    for (int jz = 0; jz < ldims[2]; jz++) {
      max_diff = std::max(
        max_diff, std::abs(f(m, 0, -1, jz) - f(m, 0, ldims[1] - 1, jz)));
      max_diff =
        std::max(max_diff, std::abs(f(m, 0, ldims[1], jz) - f(m, 0, 0, jz)));
    }
    for (int jy = 0; jy < ldims[1]; jy++) {
      max_diff = std::max(
        max_diff, std::abs(f(m, 0, jy, -1) - f(m, 0, jy, ldims[2] - 1)));
      max_diff =
        std::max(max_diff, std::abs(f(m, 0, jy, ldims[2]) - f(m, 0, jy, 0)));
    }
  }

  EXPECT_LT(max_diff, 1e-10 * max_abs_h(mflds));
}

TYPED_TEST(TurbulenceTest, IndependentOfDecomposition)
{
  auto grid_1 = setupGrid({1, 1, 1});
  MfieldsState mflds_1{*grid_1};
  this->generate(mflds_1);

  auto grid_4 = setupGrid({1, 2, 2});
  MfieldsState mflds_4{*grid_4};
  this->generate(mflds_4);
  ASSERT_EQ(mflds_4.n_patches(), 4);

  auto h_1 = gather_h(mflds_1);
  auto h_4 = gather_h(mflds_4);
  double max_diff = 0.0;
  for (int m = 0; m < 3; m++) {
    for (int i = 0; i < h_1[m].size(); i++) {
      max_diff = std::max(max_diff, std::abs(h_1[m][i] - h_4[m][i]));
    }
  }

  EXPECT_LT(max_diff, 1e-12 * max_abs_h(mflds_1));
}

TYPED_TEST(TurbulenceTest, DeterminedBySeed)
{
  auto grid = setupGrid();
  TurbulenceParams params = default_params;

  MfieldsState mflds_a{*grid};
  this->generate(mflds_a, params);
  MfieldsState mflds_b{*grid};
  this->generate(mflds_b, params);
  params.seed += 1;
  MfieldsState mflds_c{*grid};
  this->generate(mflds_c, params);

  EXPECT_EQ(gather_h(mflds_a), gather_h(mflds_b));
  EXPECT_NE(gather_h(mflds_a), gather_h(mflds_c));
}

TYPED_TEST(TurbulenceTest, ScalesWithDb2)
{
  auto grid = setupGrid();
  TurbulenceParams params = default_params;

  MfieldsState mflds_a{*grid};
  this->generate(mflds_a, params);
  params.db2 *= 4.0;
  MfieldsState mflds_b{*grid};
  this->generate(mflds_b, params);

  // same shape, twice the amplitude
  auto h_a = gather_h(mflds_a);
  auto h_b = gather_h(mflds_b);
  double max_diff = 0.0;
  for (int m = 0; m < 3; m++) {
    for (int i = 0; i < h_a[m].size(); i++) {
      max_diff = std::max(max_diff, std::abs(2.0 * h_a[m][i] - h_b[m][i]));
    }
  }
  EXPECT_LT(max_diff, 1e-12 * max_abs_h(mflds_b));
}

// a shorter correlation length puts relatively more power at high k, which
// shows up as a larger mean |curl B|^2 / |B|^2
TYPED_TEST(TurbulenceTest, ShorterCorrelationLengthIsRougher)
{
  auto grid = setupGrid();
  auto dx = grid->domain.dx;

  auto roughness = [&](double correlation_length) {
    TurbulenceParams params = default_params;
    params.correlation_length = correlation_length;
    MfieldsState mflds{*grid};
    this->generate(mflds, params);

    auto f = make_Fields3d<dim_xyz>(mflds[0]);
    double sum = 0.0;
    grid->Foreach_3d(0, 0, [&](int jx, int jy, int jz) {
      // d_y HZ - d_z HY, the x component of curl B, on x edges
      double curl_x = (f(HZ, jx, jy, jz) - f(HZ, jx, jy - 1, jz)) / dx[1] -
                      (f(HY, jx, jy, jz) - f(HY, jx, jy, jz - 1)) / dx[2];
      sum += sqr(curl_x);
    });
    return sum / grid->domain.gdims.prod() / get_mean_b2(mflds);
  };

  EXPECT_GT(roughness(0.1), 1.5 * roughness(3.0));
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
