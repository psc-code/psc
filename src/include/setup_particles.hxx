
#pragma once

#include <algorithm>
#include <functional>
#include <type_traits>

#include <mrc_profile.h>

#include "kg/VecRange.hxx"
#include "centering.hxx"
#include "particle.h"
#include "rng.hxx"

struct psc_particle_npt
{
  int kind;  ///< particle kind
  double n;  ///< density
  Double3 p; ///< momentum
  Double3 T; ///< temperature
  psc::particle::Tag tag;
};

struct psc_particle_np
{
  int kind;                   ///< particle kind
  double n;                   ///< density
  std::function<Double3()> p; ///< returns a random momentum
  psc::particle::Tag tag;
};

struct InitNptFunc
{
  // Initialize particles according to a Maxwellian.
  // (int kind, Double3 pos, psc_particle_npt& npt) -> void
  template <
    typename INIT_NPT,
    std::enable_if_t<
      std::is_convertible<
        INIT_NPT, std::function<void(int, Double3, psc_particle_npt&)>>::value,
      bool> = true>
  InitNptFunc(INIT_NPT init_npt_simple)
    : func{[&](int kind, Double3 pos, int p, Int3 idx, psc_particle_npt& npt) {
        init_npt_simple(kind, pos, npt);
      }}
  {}

  // Initialize particles according to a Maxwellian.
  // (int kind, Double3 pos, int patch, Int3 idx, psc_particle_npt& npt) -> void
  template <
    typename INIT_NPT,
    std::enable_if_t<std::is_convertible<
                       INIT_NPT, std::function<void(int, Double3, int, Int3,
                                                    psc_particle_npt&)>>::value,
                     bool> = true>
  InitNptFunc(INIT_NPT init_npt) : func{init_npt}
  {}

  void operator()(int kind, Double3 pos, int patch, Int3 index,
                  psc_particle_npt& npt)
  {
    func(kind, pos, patch, index, npt);
  }

private:
  std::function<void(int, Double3, int, Int3, psc_particle_npt&)> func;
};

struct InitNpFunc
{
  // Initialize particles according to an arbitrary momentum disribution.
  // (int kind, Double3 pos, int patch, Int3 idx, psc_particle_np& np) -> void
  template <
    typename INIT_NP,
    std::enable_if_t<std::is_convertible<
                       INIT_NP, std::function<void(int, Double3, int, Int3,
                                                   psc_particle_np&)>>::value,
                     bool> = true>
  InitNpFunc(INIT_NP init_np) : func{init_np}
  {}

  void operator()(int kind, Double3 pos, int patch, Int3 index,
                  psc_particle_np& np)
  {
    func(kind, pos, patch, index, np);
  }

private:
  std::function<void(int, Double3, int, Int3, psc_particle_np&)> func;
};

// ======================================================================
// get_n_in_cell
//
// helper function for partition / particle setup

template <typename real_t>
int get_n_in_cell(real_t density, real_t prts_per_unit_density,
                  bool fractional_n_particles_per_cell)
{
  static rng::Uniform<float> dist{0, 1};
  if (density == 0.0) {
    return 0;
  }
  if (fractional_n_particles_per_cell) {
    return density * prts_per_unit_density + dist.get();
  }
  return std::max(1, int(density * prts_per_unit_density + .5));
}

/**
 * @brief Boosts velocities using cached intermediate values for a particular
 * Lorentz frame.
 */
struct VelocityBooster
{
  VelocityBooster(Double3 frame_v)
    : frame_gamma(1.0 / std::sqrt(1.0 - frame_v.mag2())),
      frame_u(frame_v * frame_gamma),
      frame_dir(frame_v / frame_v.mag())
  {
    // FIXME kind of hacky
    if (frame_v.mag2() == 0.0) {
      frame_dir = {1, 0, 0};
    }
  }

  /**
   * @param prt_v a particle's "unprimed" proper velocity
   * @return its "primed" proper velocity
   */
  Double3 boost(Double3 prt_u)
  {
    double prt_gamma = std::sqrt(1.0 + prt_u.mag2());
    return prt_u + (frame_gamma - 1.0) * prt_u.dot(frame_dir) * frame_dir -
           frame_u * prt_gamma;
  }

  /**
   * @param prt_v a particle's "unprimed" non-proper velocity
   * @return its "primed" proper velocity
   */
  Double3 boost_and_make_proper(Double3 prt_v)
  {
    double prt_gamma = 1.0 / std::sqrt(1.0 - prt_v.mag2());
    Double3 prt_u = prt_v * prt_gamma;
    return prt_u + (frame_gamma - 1.0) * prt_u.dot(frame_dir) * frame_dir -
           frame_u * prt_gamma;
  }

  double frame_gamma;
  Double3 frame_u;
  Double3 frame_dir;
};

// ======================================================================
// SetupParticles

template <typename MP>
struct SetupParticles
{
  using Mparticles = MP;
  using real_t = typename MP::real_t;

  SetupParticles(const Grid_t& grid, int n_populations = 0)
    : kinds_{grid.kinds},
      norm_{grid.norm},
      n_populations_{n_populations},
      centerer(centering::CC)
  {
    if (n_populations_ == 0) {
      n_populations_ = kinds_.size();
    }
  }

  // ----------------------------------------------------------------------
  // get_n_in_cell

  int get_n_in_cell(real_t density)
  {
    return ::get_n_in_cell(density, real_t(norm_.prts_per_unit_density),
                           fractional_n_particles_per_cell);
  }

  // ----------------------------------------------------------------------
  // op_cellwise
  // Performs a given operation in each cell.
  // op signature: (int n_in_cell, np, Double3 pos) -> void

  template <typename OpFunc>
  void op_cellwise(const Grid_t& grid, int patch, InitNpFunc init_np,
                   OpFunc&& op)
  {
    for (Int3 index : VecRange(Int3{}, grid.ldims)) {
      Double3 pos = centerer.get_pos(grid.patches[patch], index);
      // FIXME, the issue really is that (2nd order) particle pushers
      // don't handle the invariant dim right
      for (int d = 0; d < 3; ++d) {
        if (grid.isInvar(d)) {
          pos[d] = grid.patches[patch].get_nc(index[d], d);
        }
      }

      int n_q_in_cell = 0;
      for (int pop = 0; pop < n_populations_; pop++) {
        psc_particle_np np{};
        if (pop < kinds_.size()) {
          np.kind = pop;
        }
        init_np(pop, pos, patch, index, np);

        int n_in_cell;
        if (pop != neutralizing_population) {
          n_in_cell = get_n_in_cell(np.n);
          n_q_in_cell += kinds_[np.kind].q * n_in_cell;
        } else {
          // FIXME, should handle the case where not the last population
          // is neutralizing
          assert(neutralizing_population == n_populations_ - 1);
          n_in_cell = -n_q_in_cell / kinds_[np.kind].q;
        }
        op(n_in_cell, np, pos);
      }
    }
  }

  // ----------------------------------------------------------------------
  // setupParticle

  psc::particle::Inject setupParticle(const psc_particle_np& np, Double3 pos,
                                      double weight)
  {
    return psc::particle::Inject{pos, np.p(), weight, np.kind, np.tag};
  }

  // ----------------------------------------------------------------------
  // createMaxwellian

  std::function<Double3()> createMaxwellian(const psc_particle_npt& npt)
  {
    assert(npt.kind >= 0 && npt.kind < kinds_.size());
    double beta = norm_.beta;
    double m = kinds_[npt.kind].m;

    return [=]() {
      static rng::Normal<double> dist;

      if (initial_momentum_gamma_correction) {
        // FIXME cache this (static doesn't work)
        VelocityBooster booster{-npt.p};

        Double3 prt_v;
        for (int d = 0; d < 3; d++) {
          // sample velocity in plasma frame
          prt_v[d] = dist.get(0.0, std::sqrt(npt.T[d] / m));
        }

        // boost to lab frame
        // FIXME should really sample from Maxwell-Juttner
        // this hack interprests v as u to handle rare case when v>1
        // v<<1 => v~= u anyways
        return booster.boost(prt_v);
      }

      Double3 p;
      for (int i = 0; i < 3; i++)
        p[i] = dist.get(npt.p[i], beta * std::sqrt(npt.T[i] / m));

      return p;
    };
  }

  // ----------------------------------------------------------------------
  // initNpt_to_initNp

  InitNpFunc initNpt_to_initNp(InitNptFunc& init_npt)
  {
    return InitNpFunc(
      [&](int kind, Double3 pos, int patch, Int3 idx, psc_particle_np& np) {
        psc_particle_npt npt{};
        npt.kind = np.kind;
        init_npt(kind, pos, patch, idx, npt);
        np.n = npt.n;
        np.p = createMaxwellian(npt);
        np.tag = npt.tag;
      });
  }

  // ----------------------------------------------------------------------
  // getWeight

  real_t getWeight(real_t density, int n_in_cell)
  {
    if (fractional_n_particles_per_cell) {
      return 1.;
    } else {
      return density * norm_.prts_per_unit_density / n_in_cell;
    }
  }

  // ----------------------------------------------------------------------
  // setupParticles

  void setupParticles(Mparticles& mprts, InitNptFunc init_npt)
  {
    setupParticles(mprts, initNpt_to_initNp(init_npt));
  }

  void setupParticles(Mparticles& mprts, InitNpFunc init_np)
  {
    static int pr;
    if (!pr) {
      pr = prof_register("setupp", 1., 0, 0);
    }

    prof_start(pr);
    const Grid_t& grid = mprts.grid();

    // mprts.reserve_all(n_prts_by_patch); FIXME

    auto inj = mprts.injector();

    std::vector<rng::Uniform<real_t>> offset_rngs;
    if (random_offsets) {
      int rank;
      MPI_Comm_rank(grid.comm(), &rank);
      if (rank == 0) {
        LOG_WARN(
          "SetupParticles: enabled random offsets. Initial particle positions "
          "will be randomized in each cell, instead of all at cell centers. "
          "Each species uses the same rng seed, so if there are the same "
          "number of particles of each species in each cell, each species will "
          "have the exact same initial position distribution. This results in "
          "a charge density of 0 if there are two species with opposite "
          "charges, but the resulting charge density is nonzero in general. In "
          "the latter case, take special care to ensure Gauss' law isn't "
          "violated.\n");
      }

      int seed = rng::detail::get_process_seed();
      for (int species_count = 0; species_count < grid.kinds.size();
           species_count++) {
        if (centerer.c == centering::Centering::CC) {
          offset_rngs.push_back({-0.5, 0.5, seed});
        } else {
          offset_rngs.push_back({0.0, 1.0, seed});
        }
      }
    }

    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto injector = inj[p];

      op_cellwise(
        grid, p, init_np,
        [&](int n_in_cell, psc_particle_np& np, Double3& cell_pos_cc) {
          for (int cnt = 0; cnt < n_in_cell; cnt++) {
            real_t weight = getWeight(np.n, n_in_cell);

            Double3 prt_pos = cell_pos_cc;
            if (random_offsets) {
              auto rng = offset_rngs[np.kind];
              for (int d = 0; d < 3; d++) {
                if (!grid.isInvar(d)) {
                  prt_pos[d] += rng.get() * grid.domain.dx[d];
                }
              }
            }

            auto prt = setupParticle(np, prt_pos, weight);
            injector(prt);
          }
        });
    }

    prof_stop(pr);
  }

  // ----------------------------------------------------------------------
  // partition

  std::vector<uint> partition(const Grid_t& grid, InitNptFunc init_npt)
  {
    return partition(grid, initNpt_to_initNp(init_npt));
  }

  std::vector<uint> partition(const Grid_t& grid, InitNpFunc init_np)
  {
    std::vector<uint> n_prts_by_patch(grid.n_patches());

    for (int p = 0; p < grid.n_patches(); ++p) {
      op_cellwise(grid, p, init_np,
                  [&](int n_in_cell, psc_particle_np&, Double3&) {
                    n_prts_by_patch[p] += n_in_cell;
                  });
    }

    return n_prts_by_patch;
  }

  // the initial number of particles in a cell for this population will be st so
  // that it achieves neutrality
  int neutralizing_population = {-1};
  bool fractional_n_particles_per_cell = {false};
  bool initial_momentum_gamma_correction = {false};
  bool random_offsets = false;

  centering::Centerer centerer;

private:
  const Grid_t::Kinds kinds_;
  const Grid_t::Normalization norm_;
  int n_populations_;
};
