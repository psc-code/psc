#pragma once

#include "psc.h"
#include "../axis.hxx"
#include "kg/Vec3.h"
#include "kg/VecRange.hxx"
#include "field_bc_base.hxx"
#include "field_bc_util.hxx"
#include "../psc_bnd/psc_bnd_util.hxx"

namespace psc
{
namespace bnd
{
namespace field
{

/**
 * @brief "Radiating" open boundary. Prescribe an arbitrary inflowing
 * electromagnetic pulse. See Ruhl 2006 for details.
 * @tparam Dim dimension type
 * @tparam MfieldsState fields type
 * @tparam P pulse type (i.e. a type that implements `PulseBase`)
 */
template <typename Dim, typename MfieldsState, typename P>
struct Radiating : FieldBcBase<MfieldsState>
{
  using dim_t = Dim;
  using Pulse = P;
  using real_t = typename MfieldsState::real_t;
  using Real3 = typename MfieldsState::Real3;

  Radiating(Pulse pulse, Axis d, LoHi lohi) : pulse{pulse}, d{d}, lohi{lohi} {}

  /**
   * @brief Don't do anything to the current.
   * @param mflds fields
   */
  void apply_j_bcs(MfieldsState& mflds) override {}

  /**
   * @brief Set the normal E to 0. They could be self-consistently evolved
   * instead, but ghost corners aren't handled yet. Transverse components are
   * deep enough to not affect 1st-order particles.
   * @param mflds fields
   */
  void apply_e_bcs(MfieldsState& mflds) override
  {
    const Grid_t& grid = mflds.grid();
    pulse.tick(grid.time());

    int d0 = d, d1 = d.next(), d2 = d.prev();
    int E0 = EX + d0, E1 = EX + d1, E2 = EX + d2;

    Int3 d0hat = Int3::unit(d0);
    Int3 d1hat = Int3::unit(d1);
    Int3 d2hat = Int3::unit(d2);

    for (int p = 0; p < mflds.n_patches(); p++) {
      if (lohi == Lo && grid.atBoundaryLo(p, d)) {
        auto F = make_Fields3d<dim_t>(mflds[p]);

        Int3 start = mflds.ib();
        Int3 stop = mflds.ib() + mflds.im();
        start[d0] = -1;
        stop[d0] = start[d0] + 1;

        for (Int3 i3 : VecRange(start, stop)) {
          Real3 x_e1 =
            (Real3(i3) + Real3(d1hat) * real_t(0.5)) * Real3(grid.domain.dx);
          Real3 x_e2 =
            (Real3(i3) + Real3(d2hat) * real_t(0.5)) * Real3(grid.domain.dx);

          F(E0, i3 - d0hat) = 0.0;
          F(E1, i3) = pulse.sample_exterior_field(E1, grid.time(), p, x_e1);
          F(E2, i3) = pulse.sample_exterior_field(E2, grid.time(), p, x_e2);
        }
      }
      if (lohi == Hi && grid.atBoundaryHi(p, d)) {
        auto F = make_Fields3d<dim_t>(mflds[p]);

        Int3 start = mflds.ib();
        Int3 stop = mflds.ib() + mflds.im();
        start[d0] = grid.ldims[d0] + 1;
        stop[d0] = start[d0] + 1;

        for (Int3 i3 : VecRange(start, stop)) {
          Real3 x_e1 =
            (Real3(i3) + Real3(d1hat) * real_t(0.5)) * Real3(grid.domain.dx);
          Real3 x_e2 =
            (Real3(i3) + Real3(d2hat) * real_t(0.5)) * Real3(grid.domain.dx);

          F(E0, i3) = 0.0;
          F(E1, i3) = pulse.sample_exterior_field(E1, grid.time(), p, x_e1);
          F(E2, i3) = pulse.sample_exterior_field(E2, grid.time(), p, x_e2);
        }
      }
    }
  }

  /**
   * @brief Set the first layer of transverse H ghosts such that the inflowing S
   * and P waves are prescribed at the boundary. The definitions of S and P
   * differ from Ruhl 2006 by a factor of 2, and Ruhl's definitions aren't
   * invariant under an x->y-z->x rotation.
   * @param mflds fields
   */
  void apply_h_bcs(MfieldsState& mflds) override
  {
    const Grid_t& grid = mflds.grid();
    real_t dt = grid.dt;
    Real3 dtdx = dt * Real3(grid.domain.dx_inv);

    pulse.tick(grid.time());

    int d0 = d, d1 = d.next(), d2 = d.prev();
    int H0 = HX + d0, H1 = HX + d1, H2 = HX + d2;
    int E1 = EX + d1, E2 = EX + d2;
    int J1 = JXI + d1, J2 = JXI + d2;

    Int3 d0hat = Int3::unit(d0);
    Int3 d1hat = Int3::unit(d1);
    Int3 d2hat = Int3::unit(d2);

    for (int p = 0; p < mflds.n_patches(); p++) {
      if (lohi == Lo && grid.atBoundaryLo(p, d)) {
        psc::bnd::field::detail::set_lower_ghosts_to_nan<dim_t>(mflds, p, d, HX,
                                                                false);

        auto F = make_Fields3d<dim_t>(mflds[p]);

        Int3 start = mflds.ib() + dim_t::get_noninvariant_mask();
        Int3 stop = mflds.ib() + mflds.im();
        start[d0] = 0;
        stop[d0] = start[d0] + 1;

        for (Int3 i3 : VecRange(start, stop)) {
          Real3 x3_s =
            (Real3(i3) + Real3(d1hat) * real_t(0.5)) * Real3(grid.domain.dx);
          Real3 x3_p =
            (Real3(i3) + Real3(d2hat) * real_t(0.5)) * Real3(grid.domain.dx);

          real_t pulse_s =
            pulse.sample_exterior_field(E1, grid.time(), p, x3_s) +
            pulse.sample_exterior_field(H2, grid.time(), p, x3_s);
          real_t pulse_p =
            pulse.sample_exterior_field(E2, grid.time(), p, x3_p) -
            pulse.sample_exterior_field(H1, grid.time(), p, x3_p);

          F(H2, i3 - d0hat) = (2.f * pulse_s - 2.f * F(E1, i3) -
                               dtdx[d2] * (F(H0, i3) - F(H0, i3 - d2hat)) -
                               (1.f - dtdx[d0]) * F(H2, i3) + dt * F(J1, i3)) /
                              (1.f + dtdx[d0]);
          F(H1, i3 - d0hat) = (-2.f * pulse_p + 2.f * F(E2, i3) -
                               dtdx[d1] * (F(H0, i3) - F(H0, i3 - d1hat)) -
                               (1.f - dtdx[d0]) * F(H1, i3) - dt * F(J2, i3)) /
                              (1.f + dtdx[d0]);
        }
      }

      if (lohi == Hi && grid.atBoundaryHi(p, d)) {
        psc::bnd::field::detail::set_upper_ghosts_to_nan<dim_t>(mflds, p, d, HX,
                                                                false);

        auto F = make_Fields3d<dim_t>(mflds[p]);

        Int3 start = mflds.ib() + dim_t::get_noninvariant_mask();
        Int3 stop = mflds.ib() + mflds.im();
        start[d0] = grid.ldims[d0];
        stop[d0] = start[d0] + 1;

        for (Int3 i3 : VecRange(start, stop)) {
          Real3 x3_s =
            (Real3(i3) + Real3(d1hat) * real_t(0.5)) * Real3(grid.domain.dx);
          Real3 x3_p =
            (Real3(i3) + Real3(d2hat) * real_t(0.5)) * Real3(grid.domain.dx);

          real_t pulse_s =
            pulse.sample_exterior_field(E1, grid.time(), p, x3_s) -
            pulse.sample_exterior_field(H2, grid.time(), p, x3_s);
          real_t pulse_p =
            pulse.sample_exterior_field(E2, grid.time(), p, x3_p) +
            pulse.sample_exterior_field(H1, grid.time(), p, x3_p);

          F(H2, i3) = (-2.f * pulse_s + 2.f * F(E1, i3) +
                       dtdx[d2] * (F(H0, i3) - F(H0, i3 - d2hat)) -
                       (1.f - dtdx[d0]) * F(H2, i3 - d0hat) - dt * F(J1, i3)) /
                      (1.f + dtdx[d0]);
          F(H1, i3) = (2.f * pulse_p - 2.f * F(E2, i3) +
                       dtdx[d1] * (F(H0, i3) - F(H0, i3 - d1hat)) -
                       (1.f - dtdx[d0]) * F(H1, i3 - d0hat) + dt * F(J2, i3)) /
                      (1.f + dtdx[d0]);
        }
      }
    }
  }

  Axis d;
  LoHi lohi;

private:
  Pulse pulse;
};

/**
 * @brief The archetypal Pulse type used by the `Radiating` boundary condition.
 * This type isn't used polymorphically, so extending it isn't strictly
 * required.
 * @tparam real_t real type
 */
template <typename real_t>
struct PulseBase
{
  using Real3 = Vec3<real_t>;

  /**
   * @brief Sample a component of an out-of-domain field at a location and time.
   * The field value will be used as part of a boundary condition calculation.
   * @param m field component (e.g. `EX`)
   * @param t time
   * @param p patch index
   * @param x3 cell-normalized location within the patch
   * @return the field value
   */
  virtual real_t sample_exterior_field(int m, double t, int p, Real3 x3) = 0;

  /**
   * @brief Perform any operations that occur once per time step.
   * @param t time
   */
  virtual void tick(double t) {}
};

/**
 * @brief A constant pulse. Use this for open boundaries that have constant
 * external fields.
 * @tparam real_t
 */
template <typename real_t>
struct ConstantPulse : PulseBase<real_t>
{
  using Real3 = Vec3<real_t>;

  ConstantPulse(Real3 e, Real3 h) : e{e}, h{h} {}

  real_t sample_exterior_field(int m, double t, int p, Real3 x3) override
  {
    switch (m) {
      case EX: return e[0];
      case EY: return e[1];
      case EZ: return e[2];
      case HX: return h[0];
      case HY: return h[1];
      case HZ: return h[2];
      default: return 0.0;
    }
  }

  Real3 e;
  Real3 h;
};

} // namespace field
} // namespace bnd
} // namespace psc
