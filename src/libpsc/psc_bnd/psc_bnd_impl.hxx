
#include "psc.h"
#include "fields.hxx"
#include "bnd.hxx"

#include <mrc_profile.h>
#include <mrc_ddc.h>

template<typename MF>
struct Bnd_ : BndBase
{
  using mfields_t = MF;
  using Mfields = typename MF::sub_t;
  using fields_t = typename mfields_t::fields_t;
  using real_t = typename mfields_t::real_t;
  using Fields = Fields3d<fields_t>;

  // ----------------------------------------------------------------------
  // ctor

  Bnd_(const Grid_t& grid, mrc_domain* domain, int ibn[3])
  {
    static struct mrc_ddc_funcs ddc_funcs = {
      .copy_to_buf   = copy_to_buf,
      .copy_from_buf = copy_from_buf,
      .add_from_buf  = add_from_buf,
    };

    ddc_ = mrc_domain_create_ddc(domain);
    mrc_ddc_set_funcs(ddc_, &ddc_funcs);
    mrc_ddc_set_param_int3(ddc_, "ibn", ibn);
    mrc_ddc_set_param_int(ddc_, "max_n_fields", 24);
    mrc_ddc_set_param_int(ddc_, "size_of_type", sizeof(real_t));
    mrc_ddc_setup(ddc_);
  }

  // ----------------------------------------------------------------------
  // dtor

  ~Bnd_()
  {
    mrc_ddc_destroy(ddc_);
  }
  
  // ----------------------------------------------------------------------
  // reset
  
  void reset() override
  {
    this->~Bnd_();
    new(this) Bnd_(ppsc->grid(), ppsc->mrc_domain, ppsc->ibn);
  }
  
  // ----------------------------------------------------------------------
  // add_ghosts
  
  void add_ghosts(mfields_t& mf, int mb, int me)
  {
    mrc_ddc_add_ghosts(ddc_, mb, me, &mf);
  }

  void add_ghosts(PscMfieldsBase mflds_base, int mb, int me) override
  {
    auto mf = mflds_base.get_as<mfields_t>(mb, me);
    add_ghosts(mf, mb, me);
    mf.put_as(mflds_base, mb, me);
  }

  // ----------------------------------------------------------------------
  // fill_ghosts

  void fill_ghosts(mfields_t& mf, int mb, int me)
  {
    // FIXME
    // I don't think we need as many points, and only stencil star
    // rather then box
    mrc_ddc_fill_ghosts(ddc_, mb, me, &mf);
  }

  void fill_ghosts(PscMfieldsBase mflds_base, int mb, int me) override
  {
    auto mf = mflds_base.get_as<mfields_t>(mb, me);
    fill_ghosts(mf, mb, me);
    mf.put_as(mflds_base, mb, me);
  }

  // ----------------------------------------------------------------------
  // copy_to_buf

  static void copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3],
			  void *_buf, void *ctx)
  {
    mfields_t _mf = *static_cast<mfields_t*>(ctx);
    Mfields& mf = *_mf.sub();
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
    mfields_t _mf = *static_cast<mfields_t*>(ctx);
    Mfields& mf = *_mf.sub();
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
    mfields_t _mf = *static_cast<mfields_t*>(ctx);
    Mfields& mf = *_mf.sub();
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
};
