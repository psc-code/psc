#pragma once

#include "turbulence_base.hxx"

#include "centering.hxx"
#include "rng.hxx"
#include "kg/VecRange.hxx"

#include <algorithm>
#include <limits>
#include <random>
#include <vector>

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

  /// the radial envelope, as a function of r^2
  static double window(double r2)
  {
    // the Gaussian is still exp(-pi) ~ 4% at r = 1, so window it smoothly to
    // 0 there; then only kernels from neighboring cells contribute, and the
    // noise and its gradient stay continuous across cell edges
    if (r2 >= 1.0) {
      return 0.0;
    }
    return sqr(1.0 - r2) * std::exp(-M_PI * r2);
  }

  double sample(Double3 x) const
  {
    // invariant dimensions don't count toward distance
    x *= Double3(varying);
    double r2 = x.mag2();
    if (r2 >= 1.0) {
      return 0.0;
    }
    return window(r2) * std::cos(2.0 * M_PI * x.dot(wave_dir) + wave_phase);
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

/// @brief The impulses of the noise cells [lo, hi) of one A component at one
/// level, each cell's drawn once.
///
/// A cell's impulses depend only on its hyperindex and the seed, so they come
/// out the same no matter which patch draws them. Positions and cell_extent
/// are in units of the kernel radius; cells may be longer than that.
template <typename KERNEL>
class NoiseCells
{
public:
  NoiseCells(Int2 d_ilevel, Int3 lo, Int3 hi, Double3 cell_extent,
             int n_impulses, Int5 hyperdims, int seed, KERNEL kernel)
    : lo_(lo),
      dims_(hi - lo),
      cell_extent_(cell_extent),
      n_impulses_(n_impulses),
      varying_(kernel.varying)
  {
    // TODO randomize number of impulses
    impulses_.reserve(size_t(dims_.prod()) * n_impulses);
    // row-major, like flatten_index
    for (Int3 cell_idx : VecRange(lo, hi)) {
      Int5 hyperindex{d_ilevel[0], d_ilevel[1], cell_idx[0], cell_idx[1],
                      cell_idx[2]};
      rng::Uniform<double> uniform{0.0, 1.0,
                                   get_cell_seed(hyperindex, hyperdims, seed)};
      for (int i = 0; i < n_impulses; i++) {
        Double3 pos_in_cell =
          Double3{uniform.get(0.0, 1.0), uniform.get(0.0, 1.0),
                  uniform.get(0.0, 1.0)} *
          cell_extent;
        kernel.randomize(uniform);
        impulses_.push_back({pos_in_cell, kernel});
      }
    }
  }

  /// noise at `pos`, whose cell and its neighbors along varying dims must be
  /// in [lo, hi)
  double sample(Double3 pos) const
  {
    Int3 i3 = (pos / cell_extent_).fint();
    double sum = 0.0;
    // kernels vanish beyond r = 1, which is at most one cell, so only
    // neighboring cells contribute, and only along varying dims
    for (Int3 cell_idx : VecRange(i3 - varying_, i3 + 1 + varying_)) {
      Double3 pos_in_cell = pos - Double3(cell_idx) * cell_extent_;
      const Impulse* impulses =
        &impulses_[size_t(flatten_index(cell_idx - lo_, dims_)) * n_impulses_];
      double cell_sum = 0.0;
      for (int i = 0; i < n_impulses_; i++) {
        cell_sum +=
          impulses[i].kernel.sample(pos_in_cell - impulses[i].pos_in_cell);
      }
      sum += cell_sum;
    }
    return sum;
  }

private:
  struct Impulse
  {
    Double3 pos_in_cell;
    KERNEL kernel;
  };

  Int3 lo_;
  Int3 dims_;
  Double3 cell_extent_;
  int n_impulses_;
  Int3 varying_;
  std::vector<Impulse> impulses_;
};

/// @brief Wavenumbers of the noise levels, in units of the lowest one: the
/// integers up to floor(x_max), like evenly spaced levels, but above the one
/// where a factor of 2^(1 / levels_per_octave) is a step of at least 1, evenly
/// spaced in log k instead, by the largest factor up to that which lands on
/// floor(x_max). levels_per_octave <= 0 means evenly spaced throughout.
inline std::vector<double> get_level_wavenumbers(double x_max,
                                                 double levels_per_octave)
{
  // the tolerance keeps e.g. x_max = 128 - 1e-14 from losing level 128
  int x_top = std::floor(x_max + 1e-9);
  int x_log = x_top;
  if (levels_per_octave > 0.0) {
    double ratio = std::pow(2.0, 1.0 / levels_per_octave);
    x_log = std::min(x_top, int(std::ceil(1.0 / (ratio - 1.0))));
  }

  std::vector<double> xs;
  for (int x = 1; x <= x_log; x++) {
    xs.push_back(x);
  }
  if (x_log < x_top) {
    int n_steps =
      std::ceil(levels_per_octave * std::log2(double(x_top) / x_log));
    for (int i = 1; i < n_steps; i++) {
      xs.push_back(x_log *
                   std::pow(double(x_top) / x_log, double(i) / n_steps));
    }
    xs.push_back(x_top);
  }
  return xs;
}

/// cells per dim of a level with kernel radius `s`: as many as fit in
/// `length`, and 0 along dims shorter than s, which are invariant at this
/// level, as well as along invariant dims
inline Int3 get_level_n_cells(Double3 length, Int3 varying, double s)
{
  Int3 n_cells;
  for (int d = 0; d < 3; d++) {
    n_cells[d] = varying[d] ? std::floor(length[d] / s + 1e-9) : 0;
  }
  return n_cells;
}

/// Simpson's rule for f over [a, b] with `n` (even) intervals
template <typename F>
double integrate(F f, double a, double b, int n = 128)
{
  double h = (b - a) / n;
  double sum = f(a) + f(b);
  for (int i = 1; i < n; i++) {
    sum += (i % 2 ? 4.0 : 2.0) * f(a + i * h);
  }
  return sum * h / 3.0;
}

/// @brief Expected power spectrum of noise made of KERNEL's impulses, up to a
/// factor per number of dims: |FT|^2 of window(r) cos(2 pi x.dir + phase),
/// averaged over phases and isotropic directions in n_dims dims, as a function
/// of the wavenumber q in cycles per kernel radius.
///
/// It's tabulated up to q_max, beyond which it's < 1e-8 of its peak for
/// GaborKernel, and taken as 0.
template <typename KERNEL>
class KernelPowerSpectrum
{
public:
  static constexpr double q_max = 4.0;

  KernelPowerSpectrum()
  {
    for (int n_dims = 1; n_dims <= 3; n_dims++) {
      // The window's FT along a line is the 1d FT of its projection onto the
      // line, i.e. its integral over the other dims.
      std::vector<double> projection(n_r + 1);
      for (int i = 0; i <= n_r; i++) {
        double x = double(i) / n_r;
        double rest = std::sqrt(std::max(0.0, 1.0 - x * x));
        if (n_dims == 1) {
          projection[i] = KERNEL::window(x * x);
        } else if (n_dims == 2) {
          projection[i] = 2.0 * integrate(
                                  [&](double y) {
                                    return KERNEL::window(x * x + y * y);
                                  },
                                  0.0, rest);
        } else {
          projection[i] = 2.0 * M_PI *
                          integrate(
                            [&](double rho) {
                              return KERNEL::window(x * x + rho * rho) * rho;
                            },
                            0.0, rest);
        }
      }
      // up to q_max + 1, which the direction average reaches
      int n_ft = (q_max + 1.0) * n_per_unit;
      std::vector<double> window_ft(n_ft + 1);
      for (int iq = 0; iq <= n_ft; iq++) {
        double q = double(iq) / n_per_unit;
        window_ft[iq] =
          2.0 * integrate(
                  [&](double x) {
                    return projection[std::lround(x * n_r)] *
                           std::cos(2.0 * M_PI * q * x);
                  },
                  0.0, 1.0, n_r);
      }
      auto window_ft2 = [&](double q) { return sqr(interpolate(window_ft, q)); };

      // The phase average of |FT|^2 is the mean of |window FT|^2 centered at
      // +-dir, so averaging over directions, |FT|^2 at q e is the mean of
      // |window FT|^2 at |q e - dir|.
      auto& table = tables_[n_dims - 1];
      table.resize(q_max * n_per_unit + 1);
      for (int iq = 0; iq < table.size(); iq++) {
        double q = double(iq) / n_per_unit;
        auto at_cos = [&](double c) {
          return window_ft2(std::sqrt(std::max(0.0, q * q + 1.0 - 2.0 * q * c)));
        };
        if (n_dims == 1) {
          table[iq] = 0.5 * (at_cos(1.0) + at_cos(-1.0));
        } else if (n_dims == 2) {
          table[iq] = integrate(
                        [&](double theta) { return at_cos(std::cos(theta)); },
                        0.0, M_PI, 256) /
                      M_PI;
        } else {
          table[iq] = 0.5 * integrate(at_cos, -1.0, 1.0, 256);
        }
      }
    }
  }

  double operator()(double q, int n_dims) const
  {
    return q < q_max ? interpolate(tables_[n_dims - 1], q) : 0.0;
  }

private:
  static constexpr int n_r = 128;
  static constexpr int n_per_unit = 64;

  // linear, in a table with n_per_unit points per unit, 0 beyond it
  static double interpolate(const std::vector<double>& table, double q)
  {
    double x = q * n_per_unit;
    int i = x;
    if (i + 1 >= table.size()) {
      return 0.0;
    }
    return table[i] + (x - i) * (table[i + 1] - table[i]);
  }

  std::vector<double> tables_[3];
};

/// @brief The levels' expected B spectra, in shells of width dk, and the
/// target spectrum to fit them to.
///
/// A level's spectrum is KernelPowerSpectrum times |D(k)|^2, the discrete
/// curl's, with the kernel's FT summed over aliases at k + 2 pi n / dx, since
/// the noise is sampled on the grid. It's normalized to total 1, as each level
/// is normalized to a total power.
///
/// The target is shell_power per mode, i.e. shell_power(|k|) / (modes per unit
/// k of a complete shell). In complete shells, that sums to shell_power(k) dk,
/// as for AlfvenTurbulence, but it doesn't pile the power of a shell with few
/// modes onto them. It's set only in shells below the lowest axis Nyquist
/// wavenumber; beyond that, shells are incomplete, and the levels' power there
/// is left alone.
///
/// Up to 16 times the coarsest axis dk, modes are summed one by one. Beyond,
/// they're dense enough to integrate over shells, with a fixed number of
/// directions per shell, so the cost doesn't grow with the number of cells.
struct LevelSpectra
{
  double dk;
  int n_shells;
  /// shells [1, n_fit) have a target
  int n_fit;
  /// level i's power in shell j is spectra[i * n_shells + j]
  std::vector<double> spectra;
  std::vector<double> target;

  double& at(int i_level, int j) { return spectra[i_level * n_shells + j]; }
  double at(int i_level, int j) const
  {
    return spectra[i_level * n_shells + j];
  }
};

/// @param s the levels' kernel radii
template <typename KERNEL>
LevelSpectra get_level_spectra(const std::vector<double>& s, Double3 length,
                               Int3 gdims, double correlation_length)
{
  KernelPowerSpectrum<KERNEL> kernel_power;
  int n_levels = s.size();

  Int3 varying;
  Double3 dx = length / Double3(gdims);
  Double3 dk_vec = 2.0 * M_PI / length;
  double dk = std::numeric_limits<double>::max();
  double dk_coarse = 0.0;
  double k_fit = std::numeric_limits<double>::max();
  double k_corner2 = 0.0;
  double dk_volume = 1.0;
  double dx_max = 0.0;
  int n_dims = 0;
  for (int d = 0; d < 3; d++) {
    varying[d] = gdims[d] > 1;
    if (varying[d]) {
      n_dims++;
      dk = std::min(dk, dk_vec[d]);
      dk_coarse = std::max(dk_coarse, dk_vec[d]);
      k_fit = std::min(k_fit, M_PI / dx[d]);
      k_corner2 += sqr(M_PI / dx[d]);
      dk_volume *= dk_vec[d];
      dx_max = std::max(dx_max, dx[d]);
    }
  }
  // like AlfvenTurbulence's, so modes on shell boundaries land in the upper
  auto shell_index = [&](double k) { return int(k / dk + 1e-9); };
  double sphere_area = n_dims == 1 ? 2.0 : n_dims == 2 ? 2.0 * M_PI : 4.0 * M_PI;
  auto get_mode_density = [&](double k) {
    return sphere_area * std::pow(k, n_dims - 1) / dk_volume;
  };

  LevelSpectra res;
  res.dk = dk;
  res.n_shells = shell_index(std::sqrt(k_corner2)) + 1;
  res.n_fit = std::min(res.n_shells, shell_index(k_fit));
  res.spectra.resize(n_levels * res.n_shells);
  res.target.resize(res.n_shells);
  int n_lattice =
    std::min(res.n_shells, int(std::ceil(16.0 * dk_coarse / dk - 1e-9)));

  std::vector<Int3> level_varying(n_levels);
  std::vector<int> level_n_dims(n_levels);
  std::vector<Int3> level_aliases(n_levels);
  for (int i = 0; i < n_levels; i++) {
    Int3 n_cells = get_level_n_cells(length, varying, s[i]);
    for (int d = 0; d < 3; d++) {
      level_varying[i][d] = n_cells[d] > 0;
      level_n_dims[i] += level_varying[i][d];
    }
    // aliases are at least pi / dx away, i.e. q >= s / (2 dx)
    bool aliased = s[i] / (2.0 * dx_max) < kernel_power.q_max;
    level_aliases[i] = aliased ? level_varying[i] : Int3{0, 0, 0};
  }

  // B power of level i at k_vec, which must be 0 along its invariant dims
  auto get_power = [&](int i, Double3 k_vec) {
    double curl2 = 0.0;
    for (int d = 0; d < 3; d++) {
      if (varying[d]) {
        curl2 += sqr(2.0 / dx[d] * std::sin(0.5 * k_vec[d] * dx[d]));
      }
    }
    double sum = 0.0;
    for (Int3 alias : VecRange(-level_aliases[i], level_aliases[i] + 1)) {
      Double3 k_alias = k_vec + Double3(alias) * 2.0 * M_PI / dx;
      sum += kernel_power(k_alias.mag() * s[i] / (2.0 * M_PI), level_n_dims[i]);
    }
    return curl2 * sum;
  };

  // low k: mode by mode
  Int3 n_lo, n_hi;
  for (int d = 0; d < 3; d++) {
    int n_max = varying[d] ? int(n_lattice * dk / dk_vec[d] + 1e-9) : 0;
    n_lo[d] = std::max(-n_max, -gdims[d] / 2);
    n_hi[d] = std::min(n_max, (gdims[d] - 1) / 2);
  }
  for (Int3 n : VecRange(n_lo, n_hi + 1)) {
    Double3 k_vec = Double3(n) * dk_vec;
    double k = k_vec.mag();
    int j = shell_index(k);
    if (j == 0 || j >= n_lattice) {
      continue;
    }
    if (j < res.n_fit) {
      res.target[j] += shell_power(k, correlation_length) / get_mode_density(k);
    }
    for (int i = 0; i < n_levels; i++) {
      bool in_level = true;
      for (int d = 0; d < 3; d++) {
        in_level &= level_varying[i][d] || n[d] == 0;
      }
      if (in_level) {
        res.at(i, j) += get_power(i, k_vec);
      }
    }
  }

  // high k: integrated over shells, within the grid's k-space box. Power is
  // symmetric under k -> -k, so half the directions suffice.
  std::vector<Double3> dirs;
  int d0 = -1, d1 = -1;
  for (int d = 0; d < 3; d++) {
    if (varying[d]) {
      (d0 < 0 ? d0 : d1) = d;
    }
  }
  if (n_dims == 1) {
    dirs.push_back(Double3::unit(d0));
  } else if (n_dims == 2) {
    int n_dirs = 64;
    for (int m = 0; m < n_dirs; m++) {
      double theta = M_PI * (m + 0.5) / n_dirs;
      dirs.push_back(std::cos(theta) * Double3::unit(d0) +
                     std::sin(theta) * Double3::unit(d1));
    }
  } else {
    // Fibonacci lattice on the upper hemisphere
    int n_dirs = 256;
    double golden_angle = M_PI * (3.0 - std::sqrt(5.0));
    for (int m = 0; m < n_dirs; m++) {
      double z = (m + 0.5) / n_dirs;
      double rho = std::sqrt(1.0 - z * z);
      dirs.push_back({rho * std::cos(golden_angle * m),
                      rho * std::sin(golden_angle * m), z});
    }
  }
  for (int j = n_lattice; j < res.n_shells; j++) {
    double k = (j + 0.5) * dk;
    if (j < res.n_fit) {
      res.target[j] = shell_power(k, correlation_length) * dk;
    }
    double modes_per_dir = get_mode_density(k) * dk / dirs.size();
    for (const Double3& dir : dirs) {
      Double3 k_vec = k * dir;
      bool in_grid = true;
      for (int d = 0; d < 3; d++) {
        in_grid &= std::abs(k_vec[d]) <= M_PI / dx[d];
      }
      if (!in_grid) {
        continue;
      }
      for (int i = 0; i < n_levels; i++) {
        // levels with fewer dims have power only on a lower-dimensional
        // subset of modes, which is negligible here
        if (level_varying[i] == varying) {
          res.at(i, j) += modes_per_dir * get_power(i, k_vec);
        }
      }
    }
  }

  for (int i = 0; i < n_levels; i++) {
    double total = 0.0;
    for (int j = 0; j < res.n_shells; j++) {
      total += res.at(i, j);
    }
    for (int j = 0; j < res.n_shells; j++) {
      res.at(i, j) /= total;
    }
  }
  return res;
}

/// @brief Nonnegative level weights, i.e. total B powers, whose sum of
/// spectra best matches the target, by Richardson-Lucy iteration from
/// `weights`, which must be positive.
inline std::vector<double> fit_level_weights(const LevelSpectra& spectra,
                                             std::vector<double> weights)
{
  int n_levels = weights.size();
  std::vector<double> model(spectra.n_shells);
  for (int iter = 0; iter < 10000; iter++) {
    for (int j = 1; j < spectra.n_fit; j++) {
      model[j] = 0.0;
      for (int i = 0; i < n_levels; i++) {
        model[j] += weights[i] * spectra.at(i, j);
      }
    }
    double max_change = 0.0;
    for (int i = 0; i < n_levels; i++) {
      double num = 0.0, den = 0.0;
      for (int j = 1; j < spectra.n_fit; j++) {
        if (model[j] > 0.0) {
          num += spectra.at(i, j) * spectra.target[j] / model[j];
          den += spectra.at(i, j);
        }
      }
      if (den > 0.0) {
        double factor = num / den;
        weights[i] *= factor;
        max_change = std::max(max_change, std::abs(factor - 1.0));
      }
    }
    if (max_change < 1e-6) {
      break;
    }
  }
  return weights;
}

} // namespace detail

/// @brief Turbulence as a sum of sparse-convolution noise levels, each A
/// component drawn independently, then B = curl A.
///
/// Each level's kernels are isotropic, with radius = wavelength = s. They sit
/// in a periodic lattice of cells at least s long, so the noise is periodic on
/// any box. Each level's B is normalized to mean power shell_power(k) dk,
/// where dk is the spacing of levels around it.
///
/// With levels_per_octave > 0, levels are evenly spaced in log k above the
/// wavenumber where that's sparser than even spacing in k (see
/// get_level_wavenumbers), so the cost per cell grows only with the log of the
/// box length, not linearly as with levels_per_octave = 0. The default, 4, is
/// the sparsest spacing tested whose spectrum matches even spacing's to within
/// the noise; 2 shows ripples.
/// @tparam KERNEL copyable and constructible from an Int3 of varying dims,
/// with `randomize(rng::Uniform<double>&)`, `sample(Double3) const` and a
/// `varying` member; see GaborKernel
template <typename MfieldsState, typename KERNEL = GaborKernel>
class NoiseTurbulence : public TurbulenceGenerator<MfieldsState>
{
public:
  using Base = TurbulenceGenerator<MfieldsState>;
  using Mfields = typename Base::Mfields;
  using Kernel = KERNEL;
  using Base::Base;
  using Base::params;

  /// see get_level_wavenumbers
  double levels_per_octave = 4.0;

  void generate_potential(Mfields& vector_potential) override
  {
    using detail::Int5;

    const auto& grid = vector_potential.grid();
    Int3 ibn = vector_potential.ibn();
    Double3 length = grid.domain.length;
    Int3 gdims = grid.domain.gdims;
    Int3 varying{gdims[0] > 1, gdims[1] > 1, gdims[2] > 1};
    Kernel kernel{varying};

    MfieldsState level_fields{grid};

    double max_length = 0.0;
    double min_dx = std::numeric_limits<double>::max();
    for (int d = 0; d < 3; d++) {
      if (varying[d]) {
        max_length = std::max(max_length, length[d]);
        min_dx = std::min(min_dx, double(grid.domain.dx[d]));
      }
    }
    // down to wavelengths of 2 grid cells, in units of 2 pi / max_length
    std::vector<double> xs = detail::get_level_wavenumbers(
      max_length / (2.0 * min_dx), levels_per_octave);
    int n_levels = xs.size();
    int n_dims = 3;
    double impulses_per_volume = 10.0;

    for (int i_level = 0; i_level < n_levels; i_level++) {
      double s = max_length / xs[i_level];
      double k = 2.0 * M_PI / s;
      // the width of k this level stands for, in units of 2 pi / max_length:
      // from midway to the previous level to midway to the next, and from
      // 1/2 for the first level and to 1/2 past it for the last, as with
      // evenly spaced levels (where it's always 1)
      double x_prev = i_level > 0 ? xs[i_level - 1] : 0.0;
      double x_next =
        i_level + 1 < n_levels ? xs[i_level + 1] : xs[i_level] + 1.0;
      double level_dk = 0.5 * (x_next - x_prev);

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

        for (int m : {AX, AY, AZ}) {
          // in units of s
          auto get_pos = [&](Int3 i3) {
            return centering::get_pos(patch, i3, centering::EC, m) / s;
          };

          // draw only the cells this patch's points sample from: their own
          // and their neighbors along varying dims. Positions increase with
          // the index, so the first and last points bound the cells.
          auto get_cell = [&](Int3 i3) {
            return (get_pos(i3) / cell_extent).fint();
          };
          Int3 lo = get_cell(-ibn);
          Int3 hi = get_cell(grid.ldims + ibn - 1);
          detail::NoiseCells<Kernel> cells{{m, i_level},
                                           lo - level_varying,
                                           hi + 1 + level_varying,
                                           cell_extent,
                                           n_impulses,
                                           hyperdims,
                                           params.seed,
                                           kernel};

          for (Int3 i3 : VecRange(-ibn, grid.ldims + ibn)) {
            vp_patch(m, i3[0], i3[1], i3[2]) = cells.sample(get_pos(i3));
          }
        }
      }

      // Give each level's B a mean power of exactly shell_power dk, measured
      // with the same discrete curl as the final field. dk's units are
      // absorbed by the final set_mean_b2.
      set_b_from_potential(level_fields, level_potential);
      double scale = std::sqrt(shell_power(k, params.correlation_length) *
                               level_dk / get_mean_b2(level_fields));

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
  }
};

} // namespace turbulence
} // namespace psc
