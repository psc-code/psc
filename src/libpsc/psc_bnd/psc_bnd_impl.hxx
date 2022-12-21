
#pragma once

#include "psc.h"
#include "fields.hxx"
#include "bnd.hxx"
#include "balance.hxx"

#include <mrc_profile.h>
#include <mrc_ddc.h>

template <typename S>
struct BndContext
{
  using storage_type = S;

  storage_type& mflds_gt;
  const Int3& ib;
};

template <typename MF>
struct Bnd_ : BndBase
{
  using Mfields = MF;
  using MfieldsHost = hostMirror_t<Mfields>;
  using real_t = typename Mfields::real_t;
  using storage_type = typename Mfields::Storage;
  using storage_host_type = typename MfieldsHost::Storage;
  using BndCtx = BndContext<storage_host_type>;

  // ----------------------------------------------------------------------
  // ctor

  Bnd_(const Grid_t& grid, const int ibn[3])
  {
    balance_generation_cnt_ = psc_balance_generation_cnt;
  }

  // ----------------------------------------------------------------------
  // reset

  void reset(const Grid_t& grid)
  {
    // FIXME, not really a pretty way of doing this
    this->~Bnd_();
    new (this) Bnd_(grid, grid.ibn);
  }

  // ----------------------------------------------------------------------
  // add_ghosts

  void add_ghosts(const Grid_t& grid, storage_type& mflds_gt, const Int3& ib,
                  int mb, int me)
  {
    if (psc_balance_generation_cnt != balance_generation_cnt_) {
      balance_generation_cnt_ = psc_balance_generation_cnt;
      reset(grid);
    }

    // FIXME
    // I don't think we need as many points, and only stencil star
    // rather then box
    assert(Int3(mflds_gt.shape(0), mflds_gt.shape(1), mflds_gt.shape(2)) ==
           grid.ldims + 2 * grid.ibn);

    auto&& h_mflds_gt = gt::host_mirror(mflds_gt);
    gt::copy(mflds_gt, h_mflds_gt);
    BndCtx ctx{h_mflds_gt, ib};
    mrc_ddc_set_param_int(grid.ddc(), "size_of_type", sizeof(real_t));
    mrc_ddc_set_funcs(grid.ddc(), const_cast<mrc_ddc_funcs*>(&ddc_funcs));
    mrc_ddc_add_ghosts(grid.ddc(), mb, me, &ctx);
    gt::copy(h_mflds_gt, mflds_gt);
  }

  void add_ghosts(Mfields& mflds, int mb, int me)
  {
    add_ghosts(mflds.grid(), mflds.storage(), mflds.ib(), mb, me);
  }

  // ----------------------------------------------------------------------
  // fill_ghosts

  void fill_ghosts(const Grid_t& grid, storage_type& mflds_gt, const Int3& ib,
                   int mb, int me)
  {
    if (psc_balance_generation_cnt != balance_generation_cnt_) {
      balance_generation_cnt_ = psc_balance_generation_cnt;
      reset(grid);
    }
    // FIXME
    // I don't think we need as many points, and only stencil star
    // rather then box
    assert(Int3(mflds_gt.shape(0), mflds_gt.shape(1), mflds_gt.shape(2)) ==
           grid.ldims + 2 * grid.ibn);

    auto&& h_mflds_gt = gt::host_mirror(mflds_gt);
    gt::copy(mflds_gt, h_mflds_gt);
    BndCtx ctx{h_mflds_gt, ib};
    mrc_ddc_set_param_int(grid.ddc(), "size_of_type", sizeof(real_t));
    mrc_ddc_set_funcs(grid.ddc(), const_cast<mrc_ddc_funcs*>(&ddc_funcs));
    mrc_ddc_fill_ghosts(grid.ddc(), mb, me, &ctx);
    gt::copy(h_mflds_gt, mflds_gt);
  }

  void fill_ghosts(Mfields& mflds, int mb, int me)
  {
    fill_ghosts(mflds.grid(), mflds.storage(), mflds.ib(), mb, me);
  }

  // ----------------------------------------------------------------------
  // copy_to_buf

  static void copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3],
                          void* _buf, void* _ctx)
  {
    BndCtx* ctx = static_cast<BndCtx*>(_ctx);
    real_t* buf = static_cast<real_t*>(_buf);
    const Int3& ib = ctx->ib;

    for (int m = mb; m < me; m++) {
      for (int iz = ilo[2]; iz < ihi[2]; iz++) {
        for (int iy = ilo[1]; iy < ihi[1]; iy++) {
          for (int ix = ilo[0]; ix < ihi[0]; ix++) {
            MRC_DDC_BUF3(buf, m - mb, ix, iy, iz) =
              ctx->mflds_gt(ix - ib[0], iy - ib[1], iz - ib[2], m, p);
          }
        }
      }
    }
  }

  static void add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3],
                           void* _buf, void* _ctx)
  {
    BndCtx* ctx = static_cast<BndCtx*>(_ctx);
    real_t* buf = static_cast<real_t*>(_buf);
    const Int3& ib = ctx->ib;

    for (int m = mb; m < me; m++) {
      for (int iz = ilo[2]; iz < ihi[2]; iz++) {
        for (int iy = ilo[1]; iy < ihi[1]; iy++) {
          for (int ix = ilo[0]; ix < ihi[0]; ix++) {
            ctx->mflds_gt(ix - ib[0], iy - ib[1], iz - ib[2], m, p) +=
              MRC_DDC_BUF3(buf, m - mb, ix, iy, iz);
          }
        }
      }
    }
  }

  static void copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3],
                            void* _buf, void* _ctx)
  {
    BndCtx* ctx = static_cast<BndCtx*>(_ctx);
    real_t* buf = static_cast<real_t*>(_buf);
    const Int3& ib = ctx->ib;

    for (int m = mb; m < me; m++) {
      for (int iz = ilo[2]; iz < ihi[2]; iz++) {
        for (int iy = ilo[1]; iy < ihi[1]; iy++) {
          for (int ix = ilo[0]; ix < ihi[0]; ix++) {
            ctx->mflds_gt(ix - ib[0], iy - ib[1], iz - ib[2], m, p) =
              MRC_DDC_BUF3(buf, m - mb, ix, iy, iz);
          }
        }
      }
    }
  }

  constexpr static mrc_ddc_funcs ddc_funcs = {
    .copy_to_buf = copy_to_buf,
    .copy_from_buf = copy_from_buf,
    .add_from_buf = add_from_buf,
  };

private:
  int balance_generation_cnt_;
};

template <typename MF>
constexpr mrc_ddc_funcs Bnd_<MF>::ddc_funcs;
