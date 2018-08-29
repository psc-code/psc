
#pragma once

#include "psc.h"
#include "fields.hxx"
#include "bnd.hxx"
#include "balance.hxx"

#include <mrc_profile.h>
#include <mrc_ddc.h>

template<typename MF>
struct Bnd_ : BndBase
{
  using Mfields = MF;
  using fields_t = typename Mfields::fields_t;
  using real_t = typename Mfields::real_t;
  using Fields = Fields3d<fields_t>;

  // ----------------------------------------------------------------------
  // ctor

  Bnd_(const Grid_t& grid, const int ibn[3])
  {
    static struct mrc_ddc_funcs ddc_funcs = {
      .copy_to_buf   = copy_to_buf,
      .copy_from_buf = copy_from_buf,
      .add_from_buf  = add_from_buf,
    };

    ddc_ = grid.create_ddc();
    mrc_ddc_set_funcs(ddc_, &ddc_funcs);
    mrc_ddc_set_param_int3(ddc_, "ibn", ibn);
    mrc_ddc_set_param_int(ddc_, "max_n_fields", 24);
    mrc_ddc_set_param_int(ddc_, "size_of_type", sizeof(real_t));
    assert(ibn[0] > 0 || ibn[1] > 0 || ibn[2] > 0);
    mrc_ddc_setup(ddc_);
    balance_generation_cnt_ = psc_balance_generation_cnt;
  }

  // ----------------------------------------------------------------------
  // dtor

  ~Bnd_()
  {
    mrc_ddc_destroy(ddc_);
  }
  
  // ----------------------------------------------------------------------
  // reset
  
  void reset(const Grid_t& grid)
  {
    // FIXME, not really a pretty way of doing this
    this->~Bnd_();
    new(this) Bnd_(grid, grid.ibn);
  }
  
  // ----------------------------------------------------------------------
  // add_ghosts
  
  void add_ghosts(Mfields& mflds, int mb, int me)
  {
    if (psc_balance_generation_cnt != balance_generation_cnt_) {
      reset(mflds.grid());
    }
    mrc_ddc_add_ghosts(ddc_, mb, me, &mflds);
  }

  void add_ghosts(MfieldsBase& mflds_base, int mb, int me) override
  {
    assert(0);
#if 0
    auto& mf = mflds_base.get_as<Mfields>(mb, me);
    add_ghosts(mf, mb, me);
    mflds_base.put_as(mf, mb, me);
#endif
  }

  // ----------------------------------------------------------------------
  // fill_ghosts

  void fill_ghosts(Mfields& mflds, int mb, int me)
  {
    if (psc_balance_generation_cnt != balance_generation_cnt_) {
      reset(mflds.grid());
    }
    // FIXME
    // I don't think we need as many points, and only stencil star
    // rather then box
    mrc_ddc_fill_ghosts(ddc_, mb, me, &mflds);
  }

  void fill_ghosts(MfieldsBase& mflds_base, int mb, int me) override
  {
    assert(0);
#if 0
    auto& mf = mflds_base.get_as<Mfields>(mb, me);
    fill_ghosts(mf, mb, me);
    mflds_base.put_as(mf, mb, me);
#endif
  }

  // ----------------------------------------------------------------------
  // copy_to_buf

  static void copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3],
			  void *_buf, void *ctx)
  {
    auto& mf = *static_cast<Mfields*>(ctx);
    Fields F(mf[p]);
    real_t *buf = static_cast<real_t*>(_buf);
    
    for (int m = mb; m < me; m++) {
      for (int iz = ilo[2]; iz < ihi[2]; iz++) {
	for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	  for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	    MRC_DDC_BUF3(buf, m - mb, ix,iy,iz) = F(m, ix,iy,iz);
	  }
	}
      }
    }
  }

  static void add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3],
			   void *_buf, void *ctx)
  {
    auto& mf = *static_cast<Mfields*>(ctx);
    Fields F(mf[p]);
    real_t *buf = static_cast<real_t*>(_buf);
    
    for (int m = mb; m < me; m++) {
      for (int iz = ilo[2]; iz < ihi[2]; iz++) {
	for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	  for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	    F(m, ix,iy,iz) += MRC_DDC_BUF3(buf, m - mb, ix,iy,iz);
	  }
	}
      }
    }
  }
  
  static void copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3],
			    void *_buf, void *ctx)
  {
    auto& mf = *static_cast<Mfields*>(ctx);
    Fields F(mf[p]);
    real_t *buf = static_cast<real_t*>(_buf);
    
    for (int m = mb; m < me; m++) {
      for (int iz = ilo[2]; iz < ihi[2]; iz++) {
	for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	  for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	    F(m, ix,iy,iz) = MRC_DDC_BUF3(buf, m - mb, ix,iy,iz);
	  }
	}
      }
    }
  }

private:
  mrc_ddc *ddc_;
  int balance_generation_cnt_;
};
