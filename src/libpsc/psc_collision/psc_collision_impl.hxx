
#pragma once

#include "collision.hxx"
#include "binary_collision.hxx"
#include "fields.hxx"
#include "fields3d.hxx"

#include <cmath>
#include <numeric>
#include <random>

extern void* global_collision; // FIXME

// ======================================================================
// CollisionHost

template <typename _Mparticles, typename _MfieldsState, typename _Mfields,
          typename Rng>
struct CollisionHost
{
  using Mparticles = _Mparticles;
  using real_t = typename Mparticles::real_t;
  using MfieldsState = _MfieldsState;
  using Mfields = _Mfields;
  using Particles = typename Mparticles::Patch;
  using Particle = typename Mparticles::Particle;

  enum
  {
    STATS_MIN,
    STATS_MED,
    STATS_MAX,
    STATS_NLARGE,
    STATS_NCOLL,
    NR_STATS,
  };

  struct psc_collision_stats
  {
    real_t s[NR_STATS];
  };

  CollisionHost(const Grid_t& grid, int interval, double nu)
    : interval_{interval},
      nu_{nu},
      mflds_stats_{grid, NR_STATS, grid.ibn},
      mflds_rei_{grid, NR_STATS, grid.ibn}
  {
    assert(nu_ > 0.);
    global_collision = this;
  }

  // ----------------------------------------------------------------------
  // collide

  void operator()(Mparticles& mprts)
  {
    auto& grid = mprts.grid();

    for (int p = 0; p < mprts.n_patches(); p++) {
      auto prts = mprts[p];

      const int* ldims = grid.ldims;
      int nr_cells = ldims[0] * ldims[1] * ldims[2];
      int* offsets = (int*)calloc(nr_cells + 1, sizeof(*offsets));
      struct psc_collision_stats stats_total = {};

      find_cell_offsets(mprts, p, offsets);

      auto F = make_Fields3d<dim_xyz>(mflds_stats_[p]);
      grid.Foreach_3d(0, 0, [&](int ix, int iy, int iz) {
        int c = (iz * ldims[1] + iy) * ldims[0] + ix;

        update_rei_before(mprts, p, offsets[c], offsets[c + 1], ix, iy, iz);

        struct psc_collision_stats stats = {};
        // mprintf("p %d ijk %d:%d:%d # %d\n", p, ix, iy, iz, offsets[c+1] -
        // offsets[c]);
        auto permute = randomize_in_cell(offsets[c], offsets[c + 1]);
        collide_in_cell(mprts, p, permute, &stats);

        update_rei_after(mprts, p, offsets[c], offsets[c + 1], ix, iy, iz);

        for (int s = 0; s < NR_STATS; s++) {
          F(s, ix, iy, iz) = stats.s[s];
          stats_total.s[s] += stats.s[s];
        }
      });

#if 0
      mprintf("p%d: min %g med %g max %g nlarge %g ncoll %g\n", p,
	      stats_total.s[STATS_MIN] / nr_cells,
	      stats_total.s[STATS_MED] / nr_cells,
	      stats_total.s[STATS_MAX] / nr_cells,
	      stats_total.s[STATS_NLARGE] / nr_cells,
	      stats_total.s[STATS_NCOLL] / nr_cells);
#endif

      free(offsets);
    }
  }

  // ----------------------------------------------------------------------
  // calc_stats

  static int compare(const void* _a, const void* _b)
  {
    const real_t* a = (const real_t*)_a;
    const real_t* b = (const real_t*)_b;

    if (*a < *b) {
      return -1;
    } else if (*a > *b) {
      return 1;
    } else {
      return 0;
    }
  }

  static void calc_stats(struct psc_collision_stats* stats, real_t* nudts,
                         int cnt)
  {
    qsort(nudts, cnt, sizeof(*nudts), compare);
    stats->s[STATS_NLARGE] = 0;
    for (int n = cnt - 1; n >= 0; n--) {
      if (nudts[n] < 1.) {
        break;
      }
      stats->s[STATS_NLARGE]++;
    }
    stats->s[STATS_MIN] = nudts[0];
    stats->s[STATS_MAX] = nudts[cnt - 1];
    stats->s[STATS_MED] = nudts[cnt / 2];
    stats->s[STATS_NCOLL] = cnt;
    /* mprintf("min %g med %g max %g nlarge %g ncoll %g\n", */
    /* 	  stats->s[STATS_MIN], */
    /* 	  stats->s[STATS_MED], */
    /* 	  stats->s[STATS_MAX], */
    /* 	  stats->s[STATS_NLARGE], */
    /* 	  stats->s[STATS_NCOLL]); */
  }

  // ----------------------------------------------------------------------
  // find_cell_offsets

  static void find_cell_offsets(const Mparticles& mprts, int p, int offsets[])
  {
    const int* ldims = mprts.grid().ldims;
    int last = 0;
    offsets[last] = 0;
    int n_prts = mprts.size(p);
    for (int n = 0; n < n_prts; n++) {
      int cell_index = mprts.validCellIndex(mprts[p][n].x);
      assert(cell_index >= last);
      while (last < cell_index) {
        offsets[++last] = n;
      }
    }
    while (last < ldims[0] * ldims[1] * ldims[2]) {
      offsets[++last] = n_prts;
    }
  }

  // ----------------------------------------------------------------------
  // randomize_in_cell

  std::vector<int> randomize_in_cell(int n_start, int n_end)
  {
    std::vector<int> permute(n_end - n_start);
    std::iota(permute.begin(), permute.end(), n_start);
    std::shuffle(permute.begin(), permute.end(), rng_);
    return permute;
  }

  // ----------------------------------------------------------------------
  // update_rei_before

  void update_rei_before(const Mparticles& mprts, int p, int n_start, int n_end,
                         int i, int j, int k)
  {
    real_t fnqs = mprts.grid().norm.fnqs;
    auto F = make_Fields3d<dim_xyz>(mflds_rei_[p]);
    F(0, i, j, k) = 0.;
    F(1, i, j, k) = 0.;
    F(2, i, j, k) = 0.;
    for (int n = n_start; n < n_end; n++) {
      const auto& prt = mprts[p][n];
      auto fac = mprts.prt_m(prt) * mprts.prt_w(prt) * fnqs;
      F(0, i, j, k) -= prt.u[0] * fac;
      F(1, i, j, k) -= prt.u[1] * fac;
      F(2, i, j, k) -= prt.u[2] * fac;
    }
  }

  // ----------------------------------------------------------------------
  // update_rei_after

  void update_rei_after(const Mparticles& mprts, int p, int n_start, int n_end,
                        int i, int j, int k)
  {
    real_t fnqs = mprts.grid().norm.fnqs;
    real_t dt = mprts.grid().dt;
    auto F = make_Fields3d<dim_xyz>(mflds_rei_[p]);
    for (int n = n_start; n < n_end; n++) {
      const auto& prt = mprts[p][n];
      auto fac = mprts.prt_m(prt) * mprts.prt_w(prt) * fnqs;
      F(0, i, j, k) += prt.u[0] * fac;
      F(1, i, j, k) += prt.u[1] * fac;
      F(2, i, j, k) += prt.u[2] * fac;
    }
    F(0, i, j, k) /= (this->interval_ * dt);
    F(1, i, j, k) /= (this->interval_ * dt);
    F(2, i, j, k) /= (this->interval_ * dt);
  }

  // ----------------------------------------------------------------------
  // collide_in_cell

  void collide_in_cell(Mparticles& mprts, int p,
                       const std::vector<int>& permute,
                       struct psc_collision_stats* stats)
  {
    const auto& grid = mprts.grid();
    int nn = permute.size();

    if (nn < 2) { // can't collide only one (or zero) particles
      return;
    }

    // all particles need to have same weight!
    real_t wni = mprts.prt_w(mprts[p][permute[0]]);
    real_t nudt1 = wni * grid.norm.cori * nn * this->interval_ * grid.dt * nu_;

    real_t* nudts = (real_t*)malloc((nn / 2 + 2) * sizeof(*nudts));
    int cnt = 0;

    int n = 0;
    if (nn % 2 == 1) { // odd # of particles: do 3-collision
      nudts[cnt++] = do_bc(mprts, p, permute[0], permute[1], .5 * nudt1);
      nudts[cnt++] = do_bc(mprts, p, permute[0], permute[2], .5 * nudt1);
      nudts[cnt++] = do_bc(mprts, p, permute[1], permute[2], .5 * nudt1);
      n = 3;
    }
    for (; n < nn; n += 2) { // do remaining particles as pair
      nudts[cnt++] = do_bc(mprts, p, permute[n], permute[n + 1], nudt1);
    }

    calc_stats(stats, nudts, cnt);
    free(nudts);
  }

  real_t do_bc(Mparticles& mprts, int p, int n1, int n2, real_t nudt1)
  {
    Rng rng;
    BinaryCollision<Mparticles, Particle> bc(mprts);
    return bc(mprts[p][n1], mprts[p][n2], nudt1, rng);
  }

  int interval() const { return interval_; }

private:
  // parameters
  double nu_;
  int interval_;

  std::mt19937 rng_;

public: // FIXME
  // for output
  Mfields mflds_stats_;
  Mfields mflds_rei_;
};

template <typename _Mparticles, typename _MfieldsState, typename _Mfields>
using Collision_ = CollisionHost<_Mparticles, _MfieldsState, _Mfields,
                                 RngC<typename _Mparticles::real_t>>;
