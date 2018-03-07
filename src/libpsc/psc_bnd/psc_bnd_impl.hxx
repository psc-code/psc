
#include "psc.h"
#include "fields.hxx"
#include "bnd.hxx"

#include <mrc_profile.h>
#include <mrc_ddc.h>

template<typename MF>
struct Bnd_ : BndBase
{
  using mfields_t = MF;
  using fields_t = typename mfields_t::fields_t;
  using real_t = typename mfields_t::real_t;
  using Fields = Fields3d<fields_t>;

  // ----------------------------------------------------------------------
  // create
  
  static void create(struct psc_bnd *bnd)
  {
    static struct mrc_ddc_funcs ddc_funcs = {
      .copy_to_buf   = copy_to_buf,
      .copy_from_buf = copy_from_buf,
      .add_from_buf  = add_from_buf,
    };

    struct mrc_ddc *ddc = mrc_domain_create_ddc(bnd->psc->mrc_domain);
    mrc_ddc_set_funcs(ddc, &ddc_funcs);
    mrc_ddc_set_param_int3(ddc, "ibn", bnd->psc->ibn);
    mrc_ddc_set_param_int(ddc, "max_n_fields", 24);
    mrc_ddc_set_param_int(ddc, "size_of_type", sizeof(real_t));
    mrc_ddc_setup(ddc);
    bnd->ddc = ddc;
  }

  // ----------------------------------------------------------------------
  // reset
  
  void reset()
  {
    assert(0);
    // mrc_ddc_destroy(bnd_->ddc);
    // ops->create_ddc(bnd_);
  }
  
  // ----------------------------------------------------------------------
  // add_ghosts
  
  static void add_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds_base,
			 int mb, int me)
  {
    mfields_t mf = mflds_base->get_as<MF>(mb, me);
    mrc_ddc_add_ghosts(bnd->ddc, mb, me, &mf);
    mf.put_as(mflds_base, mb, me);
  }

  // ----------------------------------------------------------------------
  // fill_ghosts

  static void fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds_base,
			  int mb, int me)
  {
    mfields_t mf = mflds_base->get_as<MF>(mb, me);
    // FIXME
    // I don't think we need as many points, and only stencil star
    // rather then box
    mrc_ddc_fill_ghosts(bnd->ddc, mb, me, &mf);
    mf.put_as(mflds_base, mb, me);
  }

  // ----------------------------------------------------------------------
  // copy_to_buf

  static void copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3],
			  void *_buf, void *ctx)
  {
    mfields_t mf = *static_cast<mfields_t*>(ctx);
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
    mfields_t mf = *static_cast<mfields_t*>(ctx);
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
    mfields_t mf = *static_cast<mfields_t*>(ctx);
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

};
