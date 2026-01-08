
#pragma once

#include "sort.hxx"

#include <psc_particles.h>

#include <mrc_profile.h>
#include <cassert>

// ======================================================================
// SortCountsort

template <typename MP>
struct SortCountsort
{
  using Mparticles = MP;
  using Particle = typename Mparticles::Particle;

  void operator()(Mparticles& mprts)
  {
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto& prts = mprts[p];
      unsigned int n_prts = prts.size();

      unsigned int n_cells = prts.pi_.n_cells_;
      unsigned int* cnts = new unsigned int[n_cells]{};

      // count
      for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
        int cni = prts.validCellIndex(*prt_iter);
        cnts[cni]++;
      }

      // calc offsets
      int cur = 0;
      for (int i = 0; i < n_cells; i++) {
        int n = cnts[i];
        cnts[i] = cur;
        cur += n;
      }
      assert(cur == n_prts);

      // move into new position
      auto particles2 = new Particle[n_prts];
      for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
        unsigned int cni = prts.validCellIndex(*prt_iter);
        particles2[cnts[cni]] = *prt_iter;
        cnts[cni]++;
      }

      // back to in-place
      memcpy(&*prts.begin(), particles2, n_prts * sizeof(*particles2));

      delete[] particles2;
      delete[] cnts;
    }
  }
};

// ======================================================================
// SortCountsort2
// use a separate array of cell indices

template <typename MP>
struct SortCountsort2
{
  using Mparticles = MP;
  using Particle = typename Mparticles::Particle;

  void operator()(Mparticles& mprts)
  {
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto&& prts = mprts[p];
      unsigned int n_prts = mprts.size(p);

      unsigned int n_cells = mprts.pi_.n_cells_;
      unsigned int* cnis = new unsigned int[n_prts];
      // FIXME, might as well merge counting here, too
      int i = 0;
      for (auto prt_iter = mprts.begin(p); prt_iter != mprts.end(p);
           ++prt_iter, ++i) {
        cnis[i] = mprts.validCellIndex(prt_iter->x);
      }

      unsigned int* cnts = new unsigned int[n_cells]{};

      // count
      for (int i = 0; i < n_prts; i++) {
        unsigned int cni = cnis[i];
        cnts[cni]++;
      }

      // calc offsets
      int cur = 0;
      for (int i = 0; i < n_cells; i++) {
        int n = cnts[i];
        cnts[i] = cur;
        cur += n;
      }
      assert(cur == n_prts);

      // move into new position
      auto particles2 = new Particle[n_prts];
      for (int i = 0; i < n_prts; i++) {
        unsigned int cni = cnis[i];
        int n = 1;
        while (i + n < n_prts && cnis[i + n] == cni) {
          n++;
        }
        memcpy(&particles2[cnts[cni]], &prts[i], n * sizeof(*particles2));
        cnts[cni] += n;
        i += n - 1;
      }

      // back to in-place
      memcpy(&*mprts.begin(p), particles2, n_prts * sizeof(*particles2));

      delete[] particles2;
      delete[] cnis;
      delete[] cnts;
    }
  }
};

// ======================================================================
// SortNone

template <typename MP>
struct SortNone
{
  using mparticles_t = MP;

  void sort(mparticles_t mprts) {}
};
