#pragma once

#include "turbulence_base.hxx"

#include "centering.hxx"
#include "rng.hxx"
#include "kg/VecRange.hxx"

#include <algorithm>
#include <limits>
#include <random>

namespace psc
{
namespace turbulence
{

/// @brief Isotropic Gabor kernel: a windowed Gaussian times a plane wave of
/// unit wavelength. Positions are in units of the kernel radius.
struct GaborKernel
{
  // varying[d] is 0 for invariant dimensions, else 1
  GaborKernel(Int3 varying) : varying(varying) {}

  void randomize(rng::Uniform<double>& uni)
  {
    // isotropic within the varying dimensions: draw uniformly in the unit
    // ball (or disk, or segment) by rejection, then normalize
    Double3 dir;
    do {
      for (int d = 0; d < 3; d++) {
        dir[d] = varying[d] ? uni.get(-1.0, 1.0) : 0.0;
      }
    } while (dir.mag2() > 1.0 || dir.mag2() == 0.0);
    wave_dir = dir / dir.mag();

    wave_phase = uni.get(0.0, 2.0 * M_PI);
  }

  double sample(Double3 x)
  {
    // invariant dimensions don't count toward distance
    x *= Double3(varying);
    // the Gaussian is still exp(-pi) ~ 4% at r = 1, so window it smoothly to
    // 0 there; then only kernels from neighboring cells contribute, and the
    // noise and its gradient stay continuous across cell edges
    double r2 = x.mag2();
    if (r2 >= 1.0) {
      return 0.0;
    }
    return sqr(1.0 - r2) * std::exp(-M_PI * r2) *
           std::cos(2.0 * M_PI * x.dot(wave_dir) + wave_phase);
  }

  Int3 varying;
  Double3 wave_dir;
  double wave_phase;
};

namespace detail
{

using Int2 = kg::Vec<int, 2>;
using Int5 = kg::Vec<int, 5>;

inline int get_cell_seed(Int5 cell_hyperidx, Int5 hyperdims, int seed)
{
  while (cell_hyperidx.min() < 0) {
    cell_hyperidx += hyperdims;
  }
  int flat_idx = flatten_index(cell_hyperidx % hyperdims, hyperdims);
  // mix, since an LCG's draws for adjacent seeds are strongly correlated
  std::seed_seq seq{static_cast<uint32_t>(flat_idx),
                    static_cast<uint32_t>(seed)};
  uint32_t mixed;
  seq.generate(&mixed, &mixed + 1);
  int cell_seed = mixed & 0x7fffffff;
  if (cell_seed == 0) {
    cell_seed = 1;
  }
  return cell_seed;
}

// Positions and cell_extent are in units of the kernel radius; cells may be
// longer than that.
template <typename KERNEL>
double sample_noise_cell(Int5 cell_hyperidx, Int5 hyperdims,
                         Double3 pos_in_cell, Double3 cell_extent,
                         int n_impulses, int seed, KERNEL& kernel)
{
  rng::Uniform<double> uniform{0.0, 1.0,
                               get_cell_seed(cell_hyperidx, hyperdims, seed)};

  // TODO randomize number of impulses

  double sum = 0.0;
  for (int i = 0; i < n_impulses; i++) {
    Double3 kernel_pos_in_cell =
      Double3{uniform.get(0.0, 1.0), uniform.get(0.0, 1.0),
              uniform.get(0.0, 1.0)} *
      cell_extent;
    kernel.randomize(uniform);
    sum += kernel.sample(pos_in_cell - kernel_pos_in_cell);
  }
  return sum;
}

template <typename KERNEL>
double sample_noise(Int2 d_ilevel, Double3 pos, Double3 cell_extent,
                    int n_impulses, Int5 hyperdims, int seed, KERNEL& kernel)
{
  Int3 i3 = (pos / cell_extent).fint();
  double sum = 0.0;
  // kernels vanish beyond r = 1, which is at most one cell, so only
  // neighboring cells contribute, and only along varying dims
  for (Int3 cell_idx : VecRange(i3 - kernel.varying, i3 + 1 + kernel.varying)) {
    Int5 hyperindex{d_ilevel[0], d_ilevel[1], cell_idx[0], cell_idx[1],
                    cell_idx[2]};
    sum += sample_noise_cell(hyperindex, hyperdims,
                             pos - Double3(cell_idx) * cell_extent, cell_extent,
                             n_impulses, seed, kernel);
  }
  return sum;
}

} // namespace detail

/// @brief Turbulence as a sum of sparse-convolution noise levels evenly spaced
/// in k, each A component drawn independently, then B = curl A.
///
/// Each level's kernels are isotropic, with radius = wavelength = s. They sit
/// in a periodic lattice of cells at least s long, so the noise is periodic on
/// any box. Each level's B is normalized to mean power shell_power(k).
/// @tparam KERNEL constructible from an Int3 of varying dims, with
/// `randomize(rng::Uniform<double>&)`, `sample(Double3)` and a `varying`
/// member; see GaborKernel
template <typename MfieldsState, typename KERNEL = GaborKernel>
class NoiseTurbulence : public TurbulenceGenerator<MfieldsState>
{
public:
  using Base = TurbulenceGenerator<MfieldsState>;
  using Mfields = ::Mfields<typename MfieldsState::real_t>;
  using Kernel = KERNEL;
  using Base::Base;
  using Base::params;

  void generate(MfieldsState& mflds) override
  {
    using detail::Int5;

    const auto& grid = mflds.grid();
    Int3 ibn = mflds.ibn();
    Double3 length = grid.domain.length;
    Int3 gdims = grid.domain.gdims;
    Int3 varying{gdims[0] > 1, gdims[1] > 1, gdims[2] > 1};
    Kernel kernel{varying};

    Mfields vector_potential{grid, 3, ibn};
    MfieldsState level_fields{grid};

    double max_length = 0.0;
    double min_dx = std::numeric_limits<double>::max();
    for (int d = 0; d < 3; d++) {
      if (varying[d]) {
        max_length = std::max(max_length, length[d]);
        min_dx = std::min(min_dx, double(grid.domain.dx[d]));
      }
    }
    // down to wavelengths of 2 grid cells
    int n_levels = std::floor(max_length / (2.0 * min_dx) + 1e-9);
    int n_dims = 3;
    double impulses_per_volume = 10.0;

    for (int i_level = 0; i_level < n_levels; i_level++) {
      double s = max_length / (i_level + 1);
      double k = 2.0 * M_PI / s;

      // Dims shorter than s can't fit this wavelength, so they are invariant
      // at this level (k_d = 0). Otherwise, fit as many cells as possible.
      Int3 n_cells;
      Int3 level_varying;
      for (int d = 0; d < 3; d++) {
        int n = varying[d] ? std::floor(length[d] / s + 1e-9) : 0;
        level_varying[d] = n > 0;
        n_cells[d] = std::max(n, 1);
      }
      Double3 cell_extent = length / Double3(n_cells) / s;
      kernel.varying = level_varying;

      // keep the kernel density, and so the noise statistics, the same at
      // every level, though cells may be longer than s and have fewer varying
      // dims
      double cell_volume = 1.0;
      for (int d = 0; d < 3; d++) {
        if (level_varying[d]) {
          cell_volume *= cell_extent[d];
        }
      }
      int n_impulses =
        std::max(1, int(std::round(impulses_per_volume * cell_volume)));

      Int5 hyperdims{n_dims, n_levels, n_cells[0], n_cells[1], n_cells[2]};
      Mfields level_potential{grid, 3, ibn};

      for (int p = 0; p < grid.n_patches(); ++p) {
        auto vp_patch = make_Fields3d<dim_xyz>(level_potential[p]);
        auto& patch = grid.patches[p];

        for (Int3 i3 : VecRange(-ibn, grid.ldims + ibn)) {
          for (int m : {AX, AY, AZ}) {
            Double3 pos_ec = centering::get_pos(patch, i3, centering::EC, m);
            vp_patch(m, i3[0], i3[1], i3[2]) =
              detail::sample_noise({m, i_level}, pos_ec / s, cell_extent,
                                   n_impulses, hyperdims, params.seed, kernel);
          }
        }
      }

      // Give each level's B a mean power of exactly shell_power, measured
      // with the same discrete curl as the final field. Levels are evenly
      // spaced in k, so dk is absorbed by the final set_mean_b2.
      set_b_from_potential(level_fields, level_potential);
      double scale = std::sqrt(shell_power(k, params.correlation_length) /
                               get_mean_b2(level_fields));

      for (int p = 0; p < grid.n_patches(); ++p) {
        auto vp_patch = make_Fields3d<dim_xyz>(vector_potential[p]);
        auto level_patch = make_Fields3d<dim_xyz>(level_potential[p]);
        for (Int3 i3 : VecRange(-ibn, grid.ldims + ibn)) {
          for (int m = 0; m < 3; m++) {
            vp_patch(m, i3[0], i3[1], i3[2]) +=
              scale * level_patch(m, i3[0], i3[1], i3[2]);
          }
        }
      }
    }

    set_b_from_potential(mflds, vector_potential);

    set_mean_b2(mflds, params.db2);
  }
};

} // namespace turbulence
} // namespace psc
