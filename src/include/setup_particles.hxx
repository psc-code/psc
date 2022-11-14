
#pragma once

#include <functional>
#include <type_traits>
#include "centering.hxx"
#include "distribution.hxx"

struct psc_particle_npt
{
  int kind;    ///< particle kind
  double n;    ///< density
  double p[3]; ///< momentum
  double T[3]; ///< temperature
  psc::particle::Tag tag;
};

struct psc_particle_np
{
  int kind;                   ///< particle kind
  double n;                   ///< density
  std::function<Double3()> p; ///< returns a random momentum
  psc::particle::Tag tag;

  psc_particle_np(const psc_particle_npt& npt, double m, double beta,
                  bool initial_momentum_gamma_correction)
    : kind{npt.kind}, n{npt.n}, tag{npt.tag}
  {
    p = [=]() {
      static distribution::Normal<double> dist;

      Double3 p;
      for (int i = 0; i < 3; i++)
        p[i] = dist.get(npt.p[i], beta * std::sqrt(npt.T[i] / m));

      if (initial_momentum_gamma_correction) {
        double p_squared = sqr(p[0]) + sqr(p[1]) + sqr(p[2]);
        if (p_squared < 1.) {
          double gamma = 1. / sqrt(1. - p_squared);
          p *= gamma;
        }
      }
      return p;
    };
  };
};

// ======================================================================
// SetupParticles

namespace
{
const Centering::Centerer defaultCenterer(Centering::CC);
}

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

template <typename MP>
struct SetupParticles
{
  using Mparticles = MP;
  using real_t = typename MP::real_t;

  SetupParticles(const Grid_t& grid, int n_populations = 0)
    : kinds_{grid.kinds},
      norm_{grid.norm},
      n_populations_{n_populations},
      centerer(defaultCenterer)
  {
    if (n_populations_ == 0) {
      n_populations_ = kinds_.size();
    }
  }

  // ----------------------------------------------------------------------
  // get_n_in_cell
  //
  // helper function for partition / particle setup

  int get_n_in_cell(const psc_particle_npt& npt)
  {
    static distribution::Uniform<float> dist{0, 1};
    if (npt.n == 0) {
      return 0;
    }
    if (fractional_n_particles_per_cell) {
      return npt.n / norm_.cori + dist.get();
    }
    return npt.n / norm_.cori + .5;
  }

  // ----------------------------------------------------------------------
  // op_cellwise
  // Performs a given operation in each cell.
  // op signature: (int n_in_cell, np, Double3 pos) -> void

  template <typename OpFunc>
  void op_cellwise(const Grid_t& grid, int patch, InitNptFunc init_npt,
                   OpFunc&& op)
  {
    Int3 ilo = Int3{}, ihi = grid.ldims;

    for (int jz = ilo[2]; jz < ihi[2]; jz++) {
      for (int jy = ilo[1]; jy < ihi[1]; jy++) {
        for (int jx = ilo[0]; jx < ihi[0]; jx++) {
          Int3 index{jx, jy, jz};

          Double3 pos = centerer.getPos(grid.patches[patch], index);
          // FIXME, the issue really is that (2nd order) particle pushers
          // don't handle the invariant dim right
          for (int a = 0; a < 3; ++a) {
            if (grid.isInvar(a) == 1) {
              pos[a] = grid.patches[patch].get_nc(index[a], a);
            }
          }

          int n_q_in_cell = 0;
          for (int pop = 0; pop < n_populations_; pop++) {
            psc_particle_npt npt{};
            if (pop < kinds_.size()) {
              npt.kind = pop;
            }
            init_npt(pop, pos, patch, index, npt);

            int n_in_cell;
            if (pop != neutralizing_population) {
              n_in_cell = get_n_in_cell(npt);
              n_q_in_cell += kinds_[npt.kind].q * n_in_cell;
            } else {
              // FIXME, should handle the case where not the last population
              // is neutralizing
              assert(neutralizing_population == n_populations_ - 1);
              n_in_cell = -n_q_in_cell / kinds_[npt.kind].q;
            }
            auto np = npt_to_np(npt);
            op(n_in_cell, np, pos);
          }
        }
      }
    }
  }

  // ----------------------------------------------------------------------
  // setupParticle

  psc::particle::Inject setupParticle(const psc_particle_np& np, Double3 pos,
                                      double wni)
  {
    return psc::particle::Inject{pos, np.p(), wni, np.kind, np.tag};
  }

  // ----------------------------------------------------------------------
  // npt_to_np

  psc_particle_np npt_to_np(const psc_particle_npt& npt)
  {
    assert(npt.kind >= 0 && npt.kind < kinds_.size());
    return psc_particle_np(npt, kinds_[npt.kind].m, norm_.beta,
                           initial_momentum_gamma_correction);
  }

  // ----------------------------------------------------------------------
  // setup_particles

  void operator()(Mparticles& mprts, InitNptFunc init_npt)
  {
    setupParticles(mprts, init_npt);
  }

  // ----------------------------------------------------------------------
  // setupParticles

  void setupParticles(Mparticles& mprts, InitNptFunc init_npt)
  {
    static int pr;
    if (!pr) {
      pr = prof_register("setupp", 1., 0, 0);
    }

    prof_start(pr);
    const auto& grid = mprts.grid();

    // mprts.reserve_all(n_prts_by_patch); FIXME

    auto inj = mprts.injector();

    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto injector = inj[p];

      op_cellwise(grid, p, init_npt,
                  [&](int n_in_cell, psc_particle_np& np, Double3& pos) {
                    for (int cnt = 0; cnt < n_in_cell; cnt++) {
                      real_t wni;
                      if (fractional_n_particles_per_cell) {
                        wni = 1.;
                      } else {
                        wni = np.n / (n_in_cell * norm_.cori);
                      }
                      auto prt = setupParticle(np, pos, wni);
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
    std::vector<uint> n_prts_by_patch(grid.n_patches());

    for (int p = 0; p < grid.n_patches(); ++p) {
      op_cellwise(grid, p, init_npt,
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

  Centering::Centerer centerer;

private:
  const Grid_t::Kinds kinds_;
  const Grid_t::Normalization norm_;
  int n_populations_;
};
