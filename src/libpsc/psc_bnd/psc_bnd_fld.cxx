
#include "psc.h"
#include "fields.hxx"

#include <mrc_profile.h>
#include <mrc_ddc.h>

using Fields = Fields3d<fields_t>;

// ======================================================================
// ddc funcs

void
psc_bnd_fld_sub_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  struct psc_mfields *mflds = reinterpret_cast<struct psc_mfields*>(ctx);
  mfields_t mf(mflds);
  Fields F(mf[p]);
  fields_t::real_t *buf = reinterpret_cast<fields_t::real_t*>(_buf);

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

void
psc_bnd_fld_sub_add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  struct psc_mfields *mflds = reinterpret_cast<struct psc_mfields*>(ctx);
  mfields_t mf(mflds);
  Fields F(mf[p]);
  fields_t::real_t *buf = reinterpret_cast<fields_t::real_t*>(_buf);

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

void
psc_bnd_fld_sub_copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  struct psc_mfields *mflds = reinterpret_cast<struct psc_mfields*>(ctx);
  mfields_t mf(mflds);
  Fields F(mf[p]);
  fields_t::real_t *buf = reinterpret_cast<fields_t::real_t*>(_buf);

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

static struct mrc_ddc_funcs ddc_funcs = {
  .copy_to_buf   = psc_bnd_fld_sub_copy_to_buf,
  .copy_from_buf = psc_bnd_fld_sub_copy_from_buf,
  .add_from_buf  = psc_bnd_fld_sub_add_from_buf,
};

// ----------------------------------------------------------------------
// psc_bnd_fld_sub_create

void
psc_bnd_fld_sub_create(struct psc_bnd *bnd)
{
  struct mrc_ddc *ddc = mrc_domain_create_ddc(bnd->psc->mrc_domain);
  mrc_ddc_set_funcs(ddc, &ddc_funcs);
  mrc_ddc_set_param_int3(ddc, "ibn", bnd->psc->ibn);
  mrc_ddc_set_param_int(ddc, "max_n_fields", 24);
  mrc_ddc_set_param_int(ddc, "size_of_type", sizeof(fields_t::real_t));
  mrc_ddc_setup(ddc);
  bnd->ddc = ddc;
}

// ----------------------------------------------------------------------
// psc_bnd_fld_sub_add_ghosts

void
psc_bnd_fld_sub_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds_base, int mb, int me)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, mb, me);
  mrc_ddc_add_ghosts(bnd->ddc, mb, me, mflds);
  psc_mfields_put_as(mflds, mflds_base, mb, me);
}

// ----------------------------------------------------------------------
// psc_bnd_fld_sub_fill_ghosts

void
psc_bnd_fld_sub_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds_base, int mb, int me)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, mb, me);
  // FIXME
  // I don't think we need as many points, and only stencil star
  // rather then box
  mrc_ddc_fill_ghosts(bnd->ddc, mb, me, mflds);
  psc_mfields_put_as(mflds, mflds_base, mb, me);
}
