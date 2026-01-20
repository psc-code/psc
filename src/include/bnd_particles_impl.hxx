
#pragma once

#include <mrc_profile.h>

#include <cstring>

#include "balance.hxx"
#include "ddc_particles.hxx"
#include "particles.hxx"
#include "particles_simple.hxx"

extern int pr_time_step_no_comm;
extern double* psc_balance_comp_time_by_patch;

// ======================================================================
// BndParticlesCommon

template <typename MP>
struct BndParticlesCommon
{
  using Mparticles = MP;
  using real_t = typename Mparticles::real_t;
  using Real3 = Vec3<real_t>;
  using ddcp_t = ddc_particles<Mparticles>;
  using BndBuffers = typename Mparticles::BndBuffers;

  // ----------------------------------------------------------------------
  // ctor

  BndParticlesCommon(const Grid_t& grid) : ddcp{}, balance_generation_cnt_{-1}
  {
    reset(grid);
  }

  // ----------------------------------------------------------------------
  // dtor

  ~BndParticlesCommon() { delete ddcp; }

  // ----------------------------------------------------------------------
  // reset

  void reset(const Grid_t& grid)
  {
    delete ddcp;
    ddcp = new ddcp_t{grid};
    balance_generation_cnt_ = psc_balance_generation_cnt;
  }

  // ----------------------------------------------------------------------
  // process_and_exchange

  void process_and_exchange(Mparticles& mprts, BndBuffers& bufs)
  {
    static int pr_B, pr_C;
    if (!pr_B) {
      pr_B = prof_register("xchg_prep", 1., 0, 0);
      pr_C = prof_register("xchg_comm", 1., 0, 0);
    }

    // prof_restart(pr_time_step_no_comm);
    prof_start(pr_B);
#pragma omp parallel for
    for (int p = 0; p < ddcp->nr_patches; p++) {
      if (psc_balance_comp_time_by_patch)
        psc_balance_comp_time_by_patch[p] -= MPI_Wtime();
      process_patch(mprts.grid(), mprts.particleIndexer(), bufs, p);
      if (psc_balance_comp_time_by_patch)
        psc_balance_comp_time_by_patch[p] += MPI_Wtime();
    }
    prof_stop(pr_B);
    // prof_stop(pr_time_step_no_comm);

    prof_start(pr_C);
    ddcp->comm(bufs);
    prof_stop(pr_C);
  }

protected:
  void process_patch(const Grid_t& grid, const ParticleIndexer<real_t>& pi,
                     BndBuffers& buf, int p);

protected:
  ddcp_t* ddcp;
  int balance_generation_cnt_;
};

// ----------------------------------------------------------------------
// BndParticlesCommon::process_patch

template <typename MP>
void BndParticlesCommon<MP>::process_patch(const Grid_t& grid,
                                           const ParticleIndexer<real_t>& pi,
                                           BndBuffers& bufs, int p)
{
  // New-style boundary requirements.
  // These will need revisiting when it comes to non-periodic domains.

  const Int3& ldims = grid.ldims;
  const auto& gpatch = grid.patches[p];
  Real3 patch_size = Real3(gpatch.xe - gpatch.xb);

  auto* dpatch = &ddcp->patches_[p];
  for (int dir1 = 0; dir1 < N_DIR; dir1++) {
    dpatch->nei[dir1].send_buf.resize(0);
  }

  using BndBuffer = typename Mparticles::BndBuffer;
  BndBuffer& buf = bufs[p];
  unsigned int n_begin = 0;
  unsigned int n_end = buf.size();
  unsigned int head = n_begin;

  for (int n = n_begin; n < n_end; n++) {
    auto* prt = &buf[n];
    real_t* xi = prt->x;
    real_t* pxi = prt->u;

    Int3 pos = pi.cellPosition(xi);

    if (pi.isValidCellPosition(pos)) {
      // fast path
      // particle is still inside patch: move into right position
      buf[head++] = *prt;
      continue;
    }

    // slow path
    // handle particles which (seemingly) left the patch
    // (may end up in the same patch, anyway, though)
    bool drop = false;
    int dir[3];
    for (int d = 0; d < 3; d++) {
      if (pos[d] < 0) {
        if (!grid.atBoundaryLo(p, d) || grid.bc.prt_lo[d] == BND_PRT_PERIODIC) {
          xi[d] += patch_size[d];
          dir[d] = -1;
          int ci = pi.cellPosition(xi[d], d);
          if (ci >= ldims[d]) {
            xi[d] = 0.;
            dir[d] = 0;
          }
        } else {
          switch (grid.bc.prt_lo[d]) {
            case BND_PRT_REFLECTING: {
              xi[d] = -xi[d];
              pxi[d] = -pxi[d];
              dir[d] = 0;
              break;
            }
            case BND_PRT_OPEN: {
              drop = true;
              break;
            }
            case BND_PRT_ABSORBING: {
              drop = true;
              break;
            }
            default: assert(0);
          }
        }
      } else if (pos[d] >= ldims[d]) {
        if (!grid.atBoundaryHi(p, d) || grid.bc.prt_hi[d] == BND_PRT_PERIODIC) {
          xi[d] -= patch_size[d];
          dir[d] = +1;
          int ci = pi.cellPosition(xi[d], d);
          if (ci < 0) {
            xi[d] = 0.;
          }
        } else {
          switch (grid.bc.prt_hi[d]) {
            case BND_PRT_REFLECTING: {
              xi[d] = 2.f * patch_size[d] - xi[d];
              pxi[d] = -pxi[d];
              dir[d] = 0;
              int ci = pi.cellPosition(xi[d], d);
              if (ci >= ldims[d]) {
                xi[d] *= (1. - 1e-6);
              }
              break;
            }
            case BND_PRT_OPEN: {
              drop = true;
              break;
            }
            case BND_PRT_ABSORBING: {
              drop = true;
              break;
            }
            default: assert(0);
          }
        }
      } else {
        // computational bnd
        dir[d] = 0;
      }
      if (!drop) {
        if (xi[d] < 0.f && xi[d] > -1e-6f) {
          mprintf("d %d xi %g\n", d, xi[d]);
          xi[d] = 0.f;
        }
        assert(xi[d] >= 0.f);
        assert(xi[d] <= patch_size[d]);
      }
    }
    if (!drop) {
      if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
        pi.validCellIndex(xi);
        buf[head++] = *prt;
      } else {
        auto* nei = &dpatch->nei[mrc_ddc_dir2idx(dir)];
        nei->send_buf.push_back(*prt);
      }
    }
  }
  buf.resize(head);
}

// ======================================================================
// BndParticles_

template <typename MP>
struct BndParticles_ : BndParticlesCommon<MP>
{
  using Base = BndParticlesCommon<MP>;
  using Mparticles = MP;

  using Base::Base;

  // ----------------------------------------------------------------------
  // operator() (exchange on correct particle type)

  void operator()(Mparticles& mprts)
  {
    if (psc_balance_generation_cnt > this->balance_generation_cnt_) {
      this->balance_generation_cnt_ = psc_balance_generation_cnt;
      this->reset(mprts.grid());
    }

    auto&& bufs = mprts.bndBuffers();
    this->process_and_exchange(mprts, bufs);

    // struct psc_mfields *mflds = psc_mfields_get_as(psc->flds, "c", JXI, JXI +
    // 3); psc_bnd_particles_open_boundary(bnd, particles, mflds);
    // psc_mfields_put_as(mflds, psc->flds, JXI, JXI + 3);
  }
};
