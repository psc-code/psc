#pragma once

#include "radiating.hxx"
#include "../axis.hxx"
#include "../bits.hxx"
#include "fields.hxx"
#include "fields3d.hxx"
#include "grid.hxx"
#include "interpolate.hxx"
#include "mpi_dtype_traits.hxx"

#include <algorithm>
#include <memory>
#include <vector>

namespace psc
{
namespace bnd
{
namespace field
{

/**
 * @brief A pulse for `Radiating` that advects a fixed snapshot of H into the
 * domain along one axis at a constant velocity, on top of constant background
 * fields. The snapshot wraps around periodically, so it should be periodic
 * along that axis.
 *
 * Each patch keeps a copy of the snapshot of whichever patch the advected
 * sample positions currently fall in. When the advection crosses a patch
 * boundary, the copies are passed on to the next patch along the axis.
 *
 * Only H is advected: sampled E is just the background. (E is sampled in
 * ghost cells, which a snapshot doesn't hold in general.)
 *
 * @tparam Dim dimension type
 * @tparam Mfields the type that holds the copies. Use the Mfields type load
 * balancing is configured with, so that it moves them with their patches.
 */
template <typename Dim, typename Mfields>
struct AdvectingPulse : PulseBase<typename Mfields::real_t>
{
  using real_t = typename Mfields::real_t;
  using Real3 = Vec3<real_t>;

  /**
   * @param fields the snapshot; only its H components are used, and they
   * must not include the background
   * @param axis the axis along which to advect
   * @param v the velocity along `axis`; positive flows in at the lower
   * boundary, negative at the upper
   * @param background_e constant E
   * @param background_h constant H, added to the snapshot's
   */
  template <typename MF>
  AdvectingPulse(MF& fields, Axis axis, real_t v, Real3 background_e,
                 Real3 background_h)
    : axis{axis},
      v_cell_normalized(v / fields.grid().domain.dx[axis]),
      background_e{background_e},
      background_h{background_h},
      held_fields(new Mfields{fields.grid(), fields.n_comps(), fields.ibn()})
  {
    held_fields->storage().view(_all, _all, _all, _all, _all) =
      fields.storage().view(_all, _all, _all, _all, _all);
  }

  const Grid_t& grid() const { return held_fields->grid(); }

  real_t sample_exterior_field(int m, double t, int p, Real3 x3) override
  {
    if (EX <= m && m <= EZ) {
      return background_e[m - EX];
    }

    Real3 x3_advected = advect(x3, t);
    int shift = wrap_into_patch(x3_advected);
    if (shift != held_shift) {
      LOG_ERROR("advected sample is %d patches away, but held fields are %d\n",
                shift, held_shift);
    }

    ip.set_coeffs(x3_advected);
    auto em = typename decltype(ip)::fields_t(
      held_fields->storage().view(_all, _all, _all, _all, p), -grid().ibn);

    switch (m) {
      case HX: return ip.hx(em) + background_h[0];
      case HY: return ip.hy(em) + background_h[1];
      case HZ: return ip.hz(em) + background_h[2];
      default: return 0.0;
    }
  }

  void tick(double t) override
  {
    // the boundary the fields flow in at
    Real3 x3 = {0.0, 0.0, 0.0};
    if (v_cell_normalized < 0.0) {
      x3[axis] = grid().ldims[axis];
    }
    Real3 x3_advected = advect(x3, t);
    int shift = wrap_into_patch(x3_advected);
    if (shift != held_shift) {
      hold(shift);
    }
  }

  /// @brief Position `x3` (cell-normalized, patch-local) was at time 0
  Real3 advect(Real3 x3, double t) const
  {
    return x3 - Real3::unit(axis) * real_t(t * v_cell_normalized);
  }

  /// @brief Wraps `x3` (cell-normalized, patch-local) into [0, ldims) along
  /// the axis, returning how many patches along the axis it moved from
  int wrap_into_patch(Real3& x3) const
  {
    int ldim = grid().ldims[axis];
    int shift = 0;
    while (x3[axis] < 0.0) {
      x3[axis] += ldim;
      shift -= 1;
    }
    while (x3[axis] >= ldim) {
      x3[axis] -= ldim;
      shift += 1;
    }
    return shift;
  }

  /// @brief Rank and local index of the patch at `idx3`, wrapped periodically
  mrc_patch_info patch_info_at(Int3 idx3) const
  {
    Int3 np = grid().domain.np;
    idx3 = (idx3 % np + np) % np;
    return grid().mrc_domain().levelIdx3PatchInfo(0, idx3);
  }

  /// @brief Exchanges copies so that each patch holds the snapshot of the
  /// patch `shift` patches away from it along the axis
  void hold(int shift)
  {
    const Grid_t& grid = this->grid();
    LOG_INFO("advecting: shifting held fields by %d patches (t=%f)\n",
             shift - held_shift, grid.time());

    // patch p receives from the patch at idx3 + delta and sends to the one at
    // idx3 - delta
    Int3 delta = Int3::unit(axis) * (shift - held_shift);
    auto& storage = held_fields->storage();
    auto received = gt::empty_like(storage);
    int n_local = grid.n_patches();
    int patch_size = storage.size() / std::max(n_local, 1);
    auto mpi_dtype = MpiDtypeTraits<real_t>::value();

    // tag = the receiving patch's local index, which is unique per
    // (sender, receiver) pair since each patch receives exactly one copy
    std::vector<MPI_Request> reqs(2 * n_local);
    for (int p = 0; p < n_local; p++) {
      Int3 idx3 = grid.localPatchInfo(p).idx3;
      auto source = patch_info_at(idx3 + delta);
      MPI_Irecv(received.data() + p * patch_size, patch_size, mpi_dtype,
                source.rank, p, grid.comm(), &reqs[p]);
    }
    for (int p = 0; p < n_local; p++) {
      Int3 idx3 = grid.localPatchInfo(p).idx3;
      auto dest = patch_info_at(idx3 - delta);
      MPI_Isend(storage.data() + p * patch_size, patch_size, mpi_dtype,
                dest.rank, dest.patch, grid.comm(), &reqs[n_local + p]);
    }
    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    storage = std::move(received);

    held_shift = shift;
  }

  Axis axis;
  real_t v_cell_normalized;
  Real3 background_e;
  Real3 background_h;
  // shared, since Radiating copies its pulse
  std::shared_ptr<Mfields> held_fields;
  // each patch holds the snapshot of the patch this many patches away
  int held_shift = 0;
  InterpolateEM<Fields3d<decltype(std::declval<Mfields&>().storage().view(
                           _all, _all, _all, _all, 0)),
                         Dim>,
                opt_ip_1st_ec, Dim>
    ip;
};

} // namespace field
} // namespace bnd
} // namespace psc
