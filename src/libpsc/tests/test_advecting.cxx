#include <gtest/gtest.h>

#include "psc_fields_c.h"
#include "../psc_bnd_fields/advecting.hxx"

#include <cmath>
#include <memory>
#include <vector>

using Dim = dim_yz;
using real_t = MfieldsC::real_t;
using Real3 = Vec3<real_t>;
using AdvectingPulse = psc::bnd::field::AdvectingPulse<Dim, MfieldsC>;

using psc::Axis;

namespace
{

const Int3 gdims{1, 16, 8};

std::unique_ptr<Grid_t> setupGrid(Int3 n_patches)
{
  // dx = 1, so velocities are in cells per unit time
  auto domain = Grid_t::Domain{gdims, Double3(gdims), {0, 0, 0}, n_patches};

  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds{};
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  Grid_t::Normalization norm{norm_params};

  Int3 ibn = 2 * Dim::get_noninvariant_mask();
  return std::unique_ptr<Grid_t>(
    new Grid_t{domain, bc, kinds, norm, 0.1, -1, ibn});
}

// A smooth snapshot, periodic over the global domain, ghosts included
double snapshot_value(int m, int jy, int jz)
{
  double ky = 2.0 * M_PI / gdims[1], kz = 2.0 * M_PI / gdims[2];
  return std::sin(ky * jy + m) + 0.5 * std::sin(2.0 * kz * jz + m) +
         0.25 * std::cos(3.0 * ky * jy);
}

MfieldsC makeSnapshot(const Grid_t& grid)
{
  MfieldsC snapshot{grid, NR_FIELDS, grid.ibn};
  Int3 ib = -grid.ibn, ie = grid.ldims + grid.ibn;
  for (int p = 0; p < grid.n_patches(); p++) {
    Int3 off = grid.patches[p].off;
    for (int m = HX; m <= HZ; m++) {
      for (int jz = ib[2]; jz < ie[2]; jz++) {
        for (int jy = ib[1]; jy < ie[1]; jy++) {
          for (int jx = ib[0]; jx < ie[0]; jx++) {
            snapshot(m, jx, jy, jz, p) =
              snapshot_value(m, off[1] + jy, off[2] + jz);
          }
        }
      }
    }
  }
  return snapshot;
}

// Samples H where Radiating does: on the inflow boundary plane (y = 0 for
// Lo, y = ldims for Hi), at z = k and k + 1/2. Returns them indexed
// [(m * gdims_z + k) * 2 + half], summed over all ranks.
std::vector<double> sampleBoundary(AdvectingPulse& pulse, const Grid_t& grid,
                                   double t, bool hi)
{
  int n_halves = 2;
  std::vector<double> local(3 * gdims[2] * n_halves);
  for (int p = 0; p < grid.n_patches(); p++) {
    if (hi ? !grid.atBoundaryHi(p, 1) : !grid.atBoundaryLo(p, 1)) {
      continue;
    }
    Int3 off = grid.patches[p].off;
    for (int m = 0; m < 3; m++) {
      for (int k = 0; k < grid.ldims[2]; k++) {
        for (int half = 0; half < n_halves; half++) {
          Real3 x3 = {0.0, hi ? real_t(grid.ldims[1]) : real_t(0.0),
                      real_t(k + 0.5 * half)};
          int i = (m * gdims[2] + off[2] + k) * n_halves + half;
          local[i] = pulse.sample_exterior_field(HX + m, t, p, x3);
        }
      }
    }
  }

  std::vector<double> global(local.size());
  MPI_Allreduce(local.data(), global.data(), local.size(), MPI_DOUBLE, MPI_SUM,
                grid.comm());
  return global;
}

// Advects the same snapshot in on a single patch and on `n_patches`, and
// expects the same boundary samples at each of a series of times, which jump
// by more than a patch and wrap around the domain more than once
void expectSameAsSinglePatch(Int3 n_patches, real_t v)
{
  bool hi = v < 0;
  Real3 background_e = {0.1, 0.2, 0.3};
  Real3 background_h = {1.0, 2.0, 3.0};

  auto grid_ref = setupGrid({1, 1, 1});
  auto snapshot_ref = makeSnapshot(*grid_ref);
  AdvectingPulse pulse_ref{snapshot_ref, Axis::Y, v, background_e,
                           background_h};

  auto grid = setupGrid(n_patches);
  auto snapshot = makeSnapshot(*grid);
  AdvectingPulse pulse{snapshot, Axis::Y, v, background_e, background_h};

  for (double t : {0.0, 0.5, 3.7, 4.4, 11.2, 25.0, 31.9, 70.3}) {
    pulse_ref.tick(t);
    pulse.tick(t);
    auto samples_ref = sampleBoundary(pulse_ref, *grid_ref, t, hi);
    auto samples = sampleBoundary(pulse, *grid, t, hi);
    for (int i = 0; i < samples.size(); i++) {
      EXPECT_NEAR(samples[i], samples_ref[i], 1e-12) << "t " << t << " i " << i;
    }
  }
}

} // namespace

TEST(AdvectingPulseTest, LoInflowMatchesSinglePatch)
{
  expectSameAsSinglePatch({1, 4, 2}, 0.9);
}

TEST(AdvectingPulseTest, LoInflowMatchesSinglePatchMany)
{
  expectSameAsSinglePatch({1, 8, 1}, 0.9);
}

TEST(AdvectingPulseTest, HiInflowMatchesSinglePatch)
{
  expectSameAsSinglePatch({1, 4, 2}, -0.9);
}

// The reference itself: on one patch, the snapshot is sampled where the
// boundary was at time 0, wrapped around the domain
TEST(AdvectingPulseTest, SinglePatchSamplesShiftedSnapshot)
{
  auto grid = setupGrid({1, 1, 1});
  auto snapshot = makeSnapshot(*grid);
  real_t v = 0.9;
  AdvectingPulse pulse{snapshot, Axis::Y, v, {}, {}};

  // shift by a whole number of cells, so that the samples land on the
  // same staggered positions as at time 0
  // on more than one rank, only one of them has the patch, but all must tick
  bool has_patch = grid->n_patches() > 0;
  double t0 = 0.0, t1 = 21.0 / v;
  pulse.tick(t0);
  std::vector<real_t> expected;
  for (int k = 0; has_patch && k < gdims[2]; k++) {
    // y = 0 - 21, wrapped around the domain, is 11
    expected.push_back(
      pulse.sample_exterior_field(HZ, t0, 0, {0.0, 11.0, real_t(k)}));
  }
  pulse.tick(t1);
  for (int k = 0; has_patch && k < gdims[2]; k++) {
    EXPECT_NEAR(pulse.sample_exterior_field(HZ, t1, 0, {0.0, 0.0, real_t(k)}),
                expected[k], 1e-12)
      << "k " << k;
  }
}

TEST(AdvectingPulseTest, Background)
{
  auto grid = setupGrid({1, 4, 2});
  MfieldsC zero{*grid, NR_FIELDS, grid->ibn};
  Real3 background_e = {0.1, 0.2, 0.3};
  Real3 background_h = {1.0, 2.0, 3.0};
  AdvectingPulse pulse{zero, Axis::Y, 0.9, background_e, background_h};

  double t = 5.0;
  pulse.tick(t);
  for (int p = 0; p < grid->n_patches(); p++) {
    Real3 x3 = {0.0, 0.0, 1.5};
    for (int m = 0; m < 3; m++) {
      EXPECT_EQ(pulse.sample_exterior_field(EX + m, t, p, x3), background_e[m]);
      EXPECT_EQ(pulse.sample_exterior_field(HX + m, t, p, x3), background_h[m]);
    }
  }
}

int main(int argc, char** argv)
{
  psc_init(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  psc_finalize();
  return rc;
}
