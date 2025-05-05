#include "gtest/gtest.h"

#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "DiagnosticsDefault.h"
#include "OutputFieldsDefault.h"
#include "add_ghosts_reflecting.hxx"
#include "../psc_config.hxx"

// ======================================================================
// PSC configuration
//
// This sets up compile-time configuration for the code, in particular
// what data structures and algorithms to use
//
// EDIT to change order / floating point type / cuda / 2d/3d

using Dim = dim_yz;
using PscConfig = PscConfig1vbecDouble<Dim>;

// ----------------------------------------------------------------------

using BgkMfields = PscConfig::Mfields;
using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;
using Marder = PscConfig::Marder;
using OutputParticles = PscConfig::OutputParticles;

// ======================================================================
// Global parameters

namespace
{
PscParams psc_params;
} // namespace

// ======================================================================
// setupGrid

Grid_t* setupGrid()
{
  auto domain = Grid_t::Domain{{1, 8, 2},          // n grid points
                               {10.0, 80.0, 20.0}, // physical lengths
                               {0, 0, 0},          // location of lower corner
                               {1, 1, 1}};         // n patches

  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_REFLECTING, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_REFLECTING, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1.0, 1.0, "e"};
  kinds[KIND_ION] = {1.0, 1.0, "i"};

  // --- generic setup
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 1;

  double dt = psc_params.cfl * courant_length(domain);
  Grid_t::Normalization norm{norm_params};

  Int3 ibn = {2, 2, 2};
  if (Dim::InvarX::value) {
    ibn[0] = 0;
  }
  if (Dim::InvarY::value) {
    ibn[1] = 0;
  }
  if (Dim::InvarZ::value) {
    ibn[2] = 0;
  }

  return new Grid_t{domain, bc, kinds, norm, dt, -1, ibn};
}

// ======================================================================
// Unit tests

int ravel_idx(Int3 idx, Int3 shape)
{
  return idx[0] + shape[0] * (idx[1] + shape[1] * (idx[2]));
}

template <typename FE>
void init_mres(FE& mres_gt, Int3 shape, int p)
{
  for (int i = 0; i < shape[0]; i++) {
    for (int j = 0; j < shape[1]; j++) {
      for (int k = 0; k < shape[2]; k++) {
        int cell_id = ravel_idx({i, j, k}, shape);
        mres_gt(i, j, k, 0, p) = cell_id;
      }
    }
  }
}

TEST(ReflectiveBcsTest, AddGhostsReflectingHighCcY)
{
  Grid_t* grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  Int3 ldims = grid.ldims;
  Int3 shape = ldims + 2 * grid.ibn;
  Int3 ib = -grid.ibn;

  auto mres = psc::mflds::zeros<double, gt::space::host_only>(grid, 1, ib);

  EXPECT_EQ(grid.n_patches(), 1);
  int p = 0;

  init_mres(mres, shape, p);

  int dim = 1;
  add_ghosts_reflecting_hi_cc(ldims, mres, ib, p, dim, 0, 1);

  int bx = -ib[0];
  int by = -ib[1];
  int bz = -ib[2];
  for (int x = -bx; x < ldims[0] + bx; x++) {
    for (int y = ldims[1] - by; y < ldims[1]; y++) {
      for (int z = -bz; z < ldims[2] + bz; z++) {
        int cell_id = ravel_idx({bx + x, by + y, bz + z}, shape);

        int reflected_y = 2 * ldims[1] - y - 1;
        int reflected_cell_id =
          ravel_idx({bx + x, by + reflected_y, bz + z}, shape);
        EXPECT_EQ(mres(bx + x, by + y, bz + z, 0, p),
                  cell_id + reflected_cell_id)
          << "xyz " << x << " " << y << " " << z;
      }
    }
  }
}

TEST(ReflectiveBcsTest, AddGhostsReflectingLowCcY)
{
  Grid_t* grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  Int3 ldims = grid.ldims;
  Int3 shape = ldims + 2 * grid.ibn;
  Int3 ib = -grid.ibn;

  auto mres = psc::mflds::zeros<double, gt::space::host_only>(grid, 1, ib);

  EXPECT_EQ(grid.n_patches(), 1);
  int p = 0;

  init_mres(mres, shape, p);

  int dim = 1;
  add_ghosts_reflecting_lo_cc(ldims, mres, ib, p, dim, 0, 1);

  int bx = -ib[0];
  int by = -ib[1];
  int bz = -ib[2];
  for (int x = -bx; x < ldims[0] + bx; x++) {
    for (int y = 0; y < by; y++) {
      for (int z = -bz; z < ldims[2] + bz; z++) {
        int cell_id = ravel_idx({bx + x, by + y, bz + z}, shape);

        int reflected_y = -y - 1;
        int reflected_cell_id =
          ravel_idx({bx + x, by + reflected_y, bz + z}, shape);
        EXPECT_EQ(mres(bx + x, by + y, bz + z, 0, p),
                  cell_id + reflected_cell_id)
          << "xyz " << x << " " << y << " " << z;
      }
    }
  }
}

TEST(ReflectiveBcsTest, AddGhostsReflectingHighCcZ)
{
  Grid_t* grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  Int3 ldims = grid.ldims;
  Int3 shape = ldims + 2 * grid.ibn;
  Int3 ib = -grid.ibn;

  auto mres = psc::mflds::zeros<double, gt::space::host_only>(grid, 1, ib);

  EXPECT_EQ(grid.n_patches(), 1);
  int p = 0;

  init_mres(mres, shape, p);

  int dim = 2;
  add_ghosts_reflecting_hi_cc(ldims, mres, ib, p, dim, 0, 1);

  int bx = -ib[0];
  int by = -ib[1];
  int bz = -ib[2];
  for (int x = -bx; x < ldims[0] + bx; x++) {
    for (int y = -by; y < ldims[1] + by; y++) {
      for (int z = ldims[2] - bz; z < ldims[2]; z++) {
        int cell_id = ravel_idx({bx + x, by + y, bz + z}, shape);

        int reflected_z = 2 * ldims[2] - z - 1;
        int reflected_cell_id =
          ravel_idx({bx + x, by + y, bz + reflected_z}, shape);
        EXPECT_EQ(mres(bx + x, by + y, bz + z, 0, p),
                  cell_id + reflected_cell_id)
          << "xyz " << x << " " << y << " " << z;
      }
    }
  }
}

TEST(ReflectiveBcsTest, AddGhostsReflectingLowCcZ)
{
  Grid_t* grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  Int3 ldims = grid.ldims;
  Int3 shape = ldims + 2 * grid.ibn;
  Int3 ib = -grid.ibn;

  auto mres = psc::mflds::zeros<double, gt::space::host_only>(grid, 1, ib);

  EXPECT_EQ(grid.n_patches(), 1);
  int p = 0;

  init_mres(mres, shape, p);

  int dim = 2;
  add_ghosts_reflecting_lo_cc(ldims, mres, ib, p, dim, 0, 1);

  int bx = -ib[0];
  int by = -ib[1];
  int bz = -ib[2];
  for (int x = -bx; x < ldims[0] + bx; x++) {
    for (int y = -by; y < ldims[1] + by; y++) {
      for (int z = 0; z < bz; z++) {
        int cell_id = ravel_idx({bx + x, by + y, bz + z}, shape);

        int reflected_z = -z - 1;
        int reflected_cell_id =
          ravel_idx({bx + x, by + y, bz + reflected_z}, shape);
        EXPECT_EQ(mres(bx + x, by + y, bz + z, 0, p),
                  cell_id + reflected_cell_id)
          << "xyz " << x << " " << y << " " << z;
      }
    }
  }
}

TEST(ReflectiveBcsTest, AddGhostsReflectingHighNcY)
{
  Grid_t* grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  Int3 ldims = grid.ldims;
  Int3 shape = ldims + 2 * grid.ibn;
  Int3 ib = -grid.ibn;

  auto mres = psc::mflds::zeros<double, gt::space::host_only>(grid, 1, ib);

  EXPECT_EQ(grid.n_patches(), 1);
  int p = 0;

  init_mres(mres, shape, p);

  int dim = 1;
  add_ghosts_reflecting_hi_nc(ldims, mres, ib, p, dim, 0, 1);

  int bx_cc = -ib[0] - !grid.isInvar(0);
  int by_cc = -ib[1] - !grid.isInvar(1);
  int bz_cc = -ib[2] - !grid.isInvar(2);
  Int3 ldims_cc =
    ldims + Int3{!grid.isInvar(0), !grid.isInvar(1), !grid.isInvar(2)};

  for (int x = -bx_cc; x < ldims_cc[0] + bx_cc; x++) {
    for (int y = ldims_cc[1] - by_cc - 1; y < ldims_cc[1] - 1; y++) {
      for (int z = -bz_cc; z < ldims_cc[2] + bz_cc; z++) {
        int cell_id = ravel_idx({x - ib[0], y - ib[1], z - ib[2]}, shape);

        int reflected_y = 2 * (ldims_cc[1] - 1) - y;
        int reflected_cell_id =
          ravel_idx({x - ib[0], reflected_y - ib[1], z - ib[2]}, shape);
        EXPECT_EQ(mres(x - ib[0], y - ib[1], z - ib[2], 0, p),
                  cell_id + reflected_cell_id)
          << "xyz " << x << " " << y << " " << z;
      }
    }
  }
}

TEST(ReflectiveBcsTest, AddGhostsReflectingLowNcY)
{
  Grid_t* grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  Int3 ldims = grid.ldims;
  Int3 shape = ldims + 2 * grid.ibn;
  Int3 ib = -grid.ibn;

  auto mres = psc::mflds::zeros<double, gt::space::host_only>(grid, 1, ib);

  EXPECT_EQ(grid.n_patches(), 1);
  int p = 0;

  init_mres(mres, shape, p);

  int dim = 1;
  add_ghosts_reflecting_lo_nc(ldims, mres, ib, p, dim, 0, 1);

  int bx_cc = -ib[0] - !grid.isInvar(0);
  int by_cc = -ib[1] - !grid.isInvar(1);
  int bz_cc = -ib[2] - !grid.isInvar(2);
  Int3 ldims_cc =
    ldims + Int3{!grid.isInvar(0), !grid.isInvar(1), !grid.isInvar(2)};

  for (int x = -bx_cc; x < ldims_cc[0] + bx_cc; x++) {
    for (int y = 1; y < 1 + by_cc; y++) {
      for (int z = -bz_cc; z < ldims_cc[2] + bz_cc; z++) {
        int cell_id = ravel_idx({x - ib[0], y - ib[1], z - ib[2]}, shape);

        int reflected_y = -y;
        int reflected_cell_id =
          ravel_idx({x - ib[0], reflected_y - ib[1], z - ib[2]}, shape);
        EXPECT_EQ(mres(x - ib[0], y - ib[1], z - ib[2], 0, p),
                  cell_id + reflected_cell_id)
          << "xyz " << x << " " << y << " " << z;
      }
    }
  }
}

TEST(ReflectiveBcsTest, AddGhostsReflectingHighNcZ)
{
  Grid_t* grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  Int3 ldims = grid.ldims;
  Int3 shape = ldims + 2 * grid.ibn;
  Int3 ib = -grid.ibn;

  auto mres = psc::mflds::zeros<double, gt::space::host_only>(grid, 1, ib);

  EXPECT_EQ(grid.n_patches(), 1);
  int p = 0;

  init_mres(mres, shape, p);

  int dim = 2;
  add_ghosts_reflecting_hi_nc(ldims, mres, ib, p, dim, 0, 1);

  int bx_cc = -ib[0] - !grid.isInvar(0);
  int by_cc = -ib[1] - !grid.isInvar(1);
  int bz_cc = -ib[2] - !grid.isInvar(2);
  Int3 ldims_cc =
    ldims + Int3{!grid.isInvar(0), !grid.isInvar(1), !grid.isInvar(2)};

  for (int x = -bx_cc; x < ldims_cc[0] + bx_cc; x++) {
    for (int y = -by_cc; y < ldims_cc[2] + by_cc; y++) {
      for (int z = ldims_cc[2] - bz_cc - 1; z < ldims_cc[2] - 1; z++) {
        int cell_id = ravel_idx({x - ib[0], y - ib[1], z - ib[2]}, shape);

        int reflected_z = 2 * (ldims_cc[2] - 1) - z;
        int reflected_cell_id =
          ravel_idx({x - ib[0], y - ib[1], reflected_z - ib[2]}, shape);
        EXPECT_EQ(mres(x - ib[0], y - ib[1], z - ib[2], 0, p),
                  cell_id + reflected_cell_id)
          << "xyz " << x << " " << y << " " << z;
      }
    }
  }
}

TEST(ReflectiveBcsTest, AddGhostsReflectingLowNcZ)
{
  Grid_t* grid_ptr = setupGrid();
  auto& grid = *grid_ptr;

  Int3 ldims = grid.ldims;
  Int3 shape = ldims + 2 * grid.ibn;
  Int3 ib = -grid.ibn;

  auto mres = psc::mflds::zeros<double, gt::space::host_only>(grid, 1, ib);

  EXPECT_EQ(grid.n_patches(), 1);
  int p = 0;

  init_mres(mres, shape, p);

  int dim = 2;
  add_ghosts_reflecting_lo_nc(ldims, mres, ib, p, dim, 0, 1);

  int bx_cc = -ib[0] - !grid.isInvar(0);
  int by_cc = -ib[1] - !grid.isInvar(1);
  int bz_cc = -ib[2] - !grid.isInvar(2);
  Int3 ldims_cc =
    ldims + Int3{!grid.isInvar(0), !grid.isInvar(1), !grid.isInvar(2)};

  for (int x = -bx_cc; x < ldims_cc[0] + bx_cc; x++) {
    for (int y = -by_cc; y < ldims_cc[1] + by_cc; y++) {
      for (int z = 1; z < 1 + bz_cc; z++) {
        int cell_id = ravel_idx({x - ib[0], y - ib[1], z - ib[2]}, shape);

        int reflected_z = -z;
        int reflected_cell_id =
          ravel_idx({x - ib[0], y - ib[1], reflected_z - ib[2]}, shape);
        EXPECT_EQ(mres(x - ib[0], y - ib[1], z - ib[2], 0, p),
                  cell_id + reflected_cell_id)
          << "xyz " << x << " " << y << " " << z;
      }
    }
  }
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
