#pragma once

#include "dim.hxx"
#include "fields.hxx"
#include "fields3d.hxx"
#include "grid.hxx"
#include "psc_bits.h"

#include <mpi.h>
#include <cmath>

namespace psc
{
namespace turbulence
{

struct TurbulenceParams
{
  /// target mean of |dB|^2 over the domain
  double db2;
  /// L_c in the shell spectrum 1 / (1 + (k L_c)^(5/3))
  double correlation_length;
  /// must be the same on every process, so B is continuous across patches
  int seed;
};

/// @brief Power of the shell at wavenumber `k`, up to normalization:
/// 1/2 <dB^2> ~ int_0^infty shell_power(k) dk. Regardless of dimensionality,
/// it goes as k^(-5/3) for k L_c >> 1, i.e., it includes the Jacobian.
inline double shell_power(double k, double correlation_length)
{
  return 1.0 / (1.0 + std::pow(k * correlation_length, 5.0 / 3.0));
}

// ======================================================================
// helpers shared by the generators

enum
{
  AX = 0,
  AY = 1,
  AZ = 2,
};

/// mean of |H|^2 over the domain's cells, each component averaged to the cell
/// center
template <typename MfieldsState>
double get_mean_b2(MfieldsState& mflds)
{
  const auto& grid = mflds.grid();

  double sum_b2_local = 0.0;

  for (int p = 0; p < mflds.n_patches(); ++p) {
    auto field_patch = make_Fields3d<dim_xyz>(mflds[p]);

    grid.Foreach_3d(0, 0, [&](int jx, int jy, int jz) {
      // note: this is not the most efficient calculation method
      sum_b2_local += 0.5 * (sqr(field_patch(HX, jx, jy, jz)) +
                             sqr(field_patch(HX, jx + 1, jy, jz)) +
                             sqr(field_patch(HY, jx, jy, jz)) +
                             sqr(field_patch(HY, jx, jy + 1, jz)) +
                             sqr(field_patch(HZ, jx, jy, jz)) +
                             sqr(field_patch(HZ, jx, jy, jz + 1)));
    });
  }

  double sum_b2 = 0.0;
  MPI_Allreduce(&sum_b2_local, &sum_b2, 1, MPI_DOUBLE, MPI_SUM, grid.comm());

  return sum_b2 / grid.domain.gdims.prod();
}

/// scales H, ghosts included, so that get_mean_b2 gives `mean_b2`
template <typename MfieldsState>
void set_mean_b2(MfieldsState& mflds, double mean_b2)
{
  const auto& grid = mflds.grid();

  double scale_factor = std::sqrt(mean_b2 / get_mean_b2(mflds));

  for (int p = 0; p < mflds.n_patches(); ++p) {
    auto field_patch = make_Fields3d<dim_xyz>(mflds[p]);

    int n_ghosts = mflds.ibn().max();
    grid.Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
      field_patch(HX, jx, jy, jz) *= scale_factor;
      field_patch(HY, jx, jy, jz) *= scale_factor;
      field_patch(HZ, jx, jy, jz) *= scale_factor;
    });
  }
}

/// sets H = curl A, with the edge-centered A in `vector_potential`'s components
/// AX, AY, AZ
template <typename MfieldsState, typename Mfields>
void set_b_from_potential(MfieldsState& mflds, Mfields& vector_potential)
{
  const auto& grid = mflds.grid();
  auto dx = grid.domain.dx;

  for (int p = 0; p < mflds.n_patches(); ++p) {
    auto field_patch = make_Fields3d<dim_xyz>(mflds[p]);
    auto a = make_Fields3d<dim_xyz>(vector_potential[p]);

    // derivative of A's component m along d, 0 along invariant dims
    auto d_a = [&](int m, int d, int jx, int jy, int jz) {
      if (grid.isInvar(d)) {
        return 0.0;
      }
      Int3 j{jx, jy, jz};
      Int3 jp = j + Int3::unit(d);
      return double(a(m, jp[0], jp[1], jp[2]) - a(m, jx, jy, jz)) / dx[d];
    };

    grid.Foreach_3d(2, 1, [&](int jx, int jy, int jz) {
      field_patch(HX, jx, jy, jz) =
        d_a(AZ, 1, jx, jy, jz) - d_a(AY, 2, jx, jy, jz);
      field_patch(HY, jx, jy, jz) =
        d_a(AX, 2, jx, jy, jz) - d_a(AZ, 0, jx, jy, jz);
      field_patch(HZ, jx, jy, jz) =
        d_a(AY, 0, jx, jy, jz) - d_a(AX, 1, jx, jy, jz);
    });
  }
}

/// @brief Generates magnetic turbulence in an MfieldsState.
///
/// `generate` overwrites H (interior and ghosts) with a divergence-free random
/// field whose mean |H|^2 over the domain is `params.db2`. E and J are left
/// alone. It is deterministic given `params.seed`, independent of the
/// decomposition.
///
/// Subclasses implement `generate_potential`, which sets the edge-centered
/// vector potential A, interior and ghosts, up to normalization; `generate`
/// takes its curl and normalizes it.
template <typename MfieldsState>
class TurbulenceGenerator
{
public:
  using Mfields = ::Mfields<typename MfieldsState::real_t>;

  explicit TurbulenceGenerator(TurbulenceParams params) : params(params) {}
  virtual ~TurbulenceGenerator() = default;

  void generate(MfieldsState& mflds)
  {
    Mfields vector_potential{mflds.grid(), 3, mflds.ibn()};
    generate_potential(vector_potential);
    set_b_from_potential(mflds, vector_potential);
    set_mean_b2(mflds, params.db2);
  }

  /// adds A to `vector_potential`, which must start at 0
  virtual void generate_potential(Mfields& vector_potential) = 0;

  TurbulenceParams params;
};

} // namespace turbulence
} // namespace psc
