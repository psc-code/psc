#pragma once

#include "turbulence_base.hxx"

#include "centering.hxx"
#include "rng.hxx"
#include "kg/VecRange.hxx"

#include <vector>

namespace psc
{
namespace turbulence
{

/// adds to `vector_potential` the potential of a transverse plane wave with
/// wavevector `k_vec` and amplitude `db` in B
template <typename Mfields>
void add_plane_alfven_wave(Mfields& vector_potential, double db,
                           const Double3& k_vec, double polarization,
                           double phase)
{
  const auto& grid = vector_potential.grid();

  double k2 = k_vec.mag2();
  double k = std::sqrt(k2);
  double kxy = std::sqrt(sqr(k_vec[0]) + sqr(k_vec[1]));

  double cos_theta = k_vec[2] / k;
  double sin_theta = kxy / k;

  double cos_phi = k_vec[0] / kxy;
  double sin_phi = k_vec[1] / kxy;

  if (kxy == 0.0) {
    // shouldn't matter given random phase and polarization, just want to
    // avoid nans
    cos_phi = 1.0;
    sin_phi = 0.0;
  }

  Double3 xp_hat{cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta};
  Double3 yp_hat{-sin_phi, cos_phi, 0};

  Double3 a_vec = db * cos(polarization) / k2 * xp_hat.cross(k_vec);
  Double3 b_vec = db * sin(polarization) / k2 * yp_hat.cross(k_vec);

  for (int p = 0; p < vector_potential.n_patches(); ++p) {
    auto& patch = grid.patches[p];
    auto vector_potential_patch = make_Fields3d<dim_xyz>(vector_potential[p]);

    int n_ghosts = vector_potential.ibn().max();
    grid.Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
      Int3 index{jx, jy, jz};

      Double3 pos_ec_x = centering::get_pos(patch, index, centering::EC, 0);
      Double3 pos_ec_y = centering::get_pos(patch, index, centering::EC, 1);
      Double3 pos_ec_z = centering::get_pos(patch, index, centering::EC, 2);

      double arg_x = k_vec.dot(pos_ec_x) + phase;
      double arg_y = k_vec.dot(pos_ec_y) + phase;
      double arg_z = k_vec.dot(pos_ec_z) + phase;

      vector_potential_patch(AX, jx, jy, jz) +=
        a_vec[0] * sin(arg_x) - b_vec[0] * cos(arg_x);
      vector_potential_patch(AY, jx, jy, jz) +=
        a_vec[1] * sin(arg_y) - b_vec[1] * cos(arg_y);
      vector_potential_patch(AZ, jx, jy, jz) +=
        a_vec[2] * sin(arg_z) - b_vec[2] * cos(arg_z);
    });
  }
}

/// @brief Turbulence as a sum of transverse plane waves, one for every
/// wavevector the periodic domain supports, each at a random phase and
/// polarization. Shells of width dk = min(2 pi / L) share shell_power(k) dk
/// equally among their modes.
template <typename MfieldsState>
class AlfvenTurbulence : public TurbulenceGenerator<MfieldsState>
{
public:
  using Base = TurbulenceGenerator<MfieldsState>;
  using Mfields = typename Base::Mfields;
  using Base::Base;
  using Base::params;

  void generate_potential(Mfields& vector_potential) override
  {
    const auto& grid = vector_potential.grid();
    Int3 gdims = grid.domain.gdims;

    // +n and -n are the same real mode, so inject only one of them: the one
    // whose first nonzero component is positive
    auto inject_at_n = [&](Int3 n) {
      for (int d = 0; d < 3; d++) {
        if (n[d] != 0) {
          return n[d] > 0;
        }
      }
      return false;
    };

    // 1. compute values of |k|

    Double3 dk_vec = 2.0 * M_PI / Double3(grid.domain.length);

    Int3 i3_min = (1 - gdims) / 2;
    Int3 i3_max = gdims / 2;

    // skip the half of k-space that inject_at_n rejects anyway
    for (int d = 0; d < 3; d++) {
      if (gdims[d] > 2) {
        i3_min[d] = 0;
        break;
      }
    }

    double dk = dk_vec.min();
    // modes on shell boundaries (e.g. on-axis) give k / dk just below an
    // integer, e.g. 3 * dk -> 2.9999999999999996, so truncating alone would
    // put them in the shell below
    auto shell_index = [&](double k) { return int(k / dk + 1e-9); };
    double kmax = ((Double3)Int3::max(-i3_min, i3_max) * dk_vec).mag();
    int nk = shell_index(kmax) + 1;

    int rank;
    MPI_Comm_rank(grid.comm(), &rank);
    if (rank == 0) {
      LOG_INFO("nk = %d, kmax = %f\n", nk, kmax);
    }

    auto k_mins = std::vector<double>(nk);
    for (int i = 0; i < nk; i++) {
      k_mins[i] = dk * i;
    }

    // 2. compute number of degeneracies of each k-shell

    auto n_cells_per_shell = std::vector<int>(nk);

    for (Int3 i3 : VecRange(i3_min, i3_max)) {
      if (!inject_at_n(i3)) {
        continue;
      }

      Double3 k_vec = Double3(i3) * dk_vec;
      double k = k_vec.mag();
      int idx = shell_index(k);

      if (idx >= nk) {
        LOG_ERROR("k=%f > kmax at idx=(%d, %d, %d)\n", k, i3[0], i3[1], i3[2]);
      }

      n_cells_per_shell[idx] += 1;
    }

    // 3. compute powers for each shell

    auto shell_db2s = std::vector<double>(nk);

    for (int i = 0; i < nk; i++) {
      if (n_cells_per_shell[i] > 0) {
        double k_shell = k_mins[i] + 0.5 * dk;
        shell_db2s[i] = shell_power(k_shell, params.correlation_length) * dk;
      }
    }

    // 4. inject each mode at a random phase and polarization

    auto rng = rng::Uniform<double>(0.0, 1.0, params.seed);

    for (Int3 i3 : VecRange(i3_min, i3_max)) {
      if (!inject_at_n(i3)) {
        continue;
      }

      Double3 k_vec = Double3(i3) * dk_vec;
      int idx = shell_index(k_vec.mag());

      double proportion_of_shell_db2_in_cell = 1.0 / n_cells_per_shell[idx];
      double shell_db2 = shell_db2s[idx];
      double cell_db2 = proportion_of_shell_db2_in_cell * shell_db2;
      double db = sqrt(cell_db2);

      double phase = rng.get(0.0, 2.0 * M_PI);
      double polarization = rng.get(0.0, 2.0 * M_PI);

      add_plane_alfven_wave(vector_potential, db, k_vec, polarization, phase);
    }
  }
};

} // namespace turbulence
} // namespace psc
