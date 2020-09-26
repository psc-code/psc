
#pragma once

#include "vpic_iface.h"

#ifdef USE_VPIC

// ======================================================================
// SortVpicWrap
//
// wraps the actual vpic sort; only works with Mparticles which themselves wrap
// the vpic particles

template <typename Mparticles>
struct SortVpicWrap
{
  void operator()(Mparticles& mprts)
  {
    auto step = mprts.grid().timestep();
    // Sort the particles for performance if desired.

    for (auto& sp : mprts[0]) {
      if (sp.sort_interval > 0 && (step % sp.sort_interval) == 0) {
        mpi_printf(MPI_COMM_WORLD, "Performance sorting \"%s\"\n", sp.name);
        TIC ::sort_p(&sp);
        TOC(sort_p, 1);
      }
    }
  }
};

#endif

// ======================================================================
// SortVpic
//
// sorts vpic-style particles

template <typename Mparticles>
struct SortVpic
{
  using Grid = typename Mparticles::Grid;
  using Species = typename Mparticles::Species;
  using Particle = typename Mparticles::Particle;

  // ----------------------------------------------------------------------
  // operator()

  void operator()(Mparticles& mprts)
  {
    auto step = mprts.grid().timestep();
    // Sort the particles for performance if desired.

    for (auto& sp : mprts[0]) {
      if (sp.sort_interval > 0 && (step % sp.sort_interval) == 0) {
        mpi_printf(MPI_COMM_WORLD, "Performance sorting \"%s\"\n", sp.name);
        TIC sort_p(sp);
        TOC(sort_p, 1);
      }
    }
  }

private:
  // ----------------------------------------------------------------------
  // sort_p

  static void sort_p(Species& sp)
  {
    const auto& g = sp.vgrid();
    sp.last_sorted = g.step;

    int n_prts = sp.np;
    int vl = VOXEL(1, 1, 1, g.nx, g.ny, g.nz);
    int vh = VOXEL(g.nx, g.ny, g.nz, g.nx, g.ny, g.nz) + 1;

    static int* RESTRICT ALIGNED(128) next;
    if (!next) {
      next = new int[g.nv];
    }
    int* RESTRICT ALIGNED(128) partition = sp.partition;

    static Particle* RESTRICT ALIGNED(128) p_aux;
    static size_t n_alloced;
    if (n_prts > n_alloced) {
      delete[] p_aux;
      p_aux = new Particle[n_prts];
      n_alloced = n_prts;
    }
    Particle* RESTRICT ALIGNED(128) p = sp.p;

    // zero counts
    for (int v = vl; v < vh; v++) {
      next[v] = 0;
    }

    // find counts
    for (int i = 0; i < n_prts; i++) {
      next[p[i].i]++;
    }

    // prefix sum
    int sum = 0;
    for (int v = vl; v < vh; v++) {
      int count = next[v];
      next[v] = sum;
      partition[v] = sum;
      sum += count;
    }
    partition[vh] = sum;

    // reorder
    for (int i = 0; i < n_prts; i++) {
      int v = p[i].i;
      int j = next[v]++;
      p_aux[j] = p[i];
    }

    // fix up unused part of partition
    for (int i = 0; i < vl; i++) {
      partition[i] = 0;
    }
    for (int i = vh; i < g.nv; i++) {
      partition[i] = n_prts;
    }

    // OPT: just swap pointer?
    for (int i = 0; i < n_prts; i++) {
      p[i] = p_aux[i];
    }
  }
};
