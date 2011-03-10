
#include "psc_bnd_private.h"

#include "psc.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>

#define to_psc_bnd_c(bnd) ((struct psc_bnd_c *)((bnd)->obj.subctx))

void create_bnd(void);

// ======================================================================
// ddc funcs

static void
copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  mfields_base_t *flds = ctx;
  fields_base_t *pf = &flds->f[p];
  fields_base_real_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_DDC_BUF3(buf, m - mb, ix,iy,iz) = F3_BASE(pf, m, ix,iy,iz);
	}
      }
    }
  }
}

static void
add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  mfields_base_t *flds = ctx;
  fields_base_t *pf = &flds->f[p];
  fields_base_real_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  F3_BASE(pf, m, ix,iy,iz) += MRC_DDC_BUF3(buf, m - mb, ix,iy,iz);
	}
      }
    }
  }
}

static void
copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *ctx)
{
  mfields_base_t *flds = ctx;
  fields_base_t *pf = &flds->f[p];
  fields_base_real_t *buf = _buf;

  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  F3_BASE(pf, m, ix,iy,iz) = MRC_DDC_BUF3(buf, m - mb, ix,iy,iz);
	}
      }
    }
  }
}

struct mrc_ddc_funcs ddc_funcs = {
  .copy_to_buf   = copy_to_buf,
  .copy_from_buf = copy_from_buf,
  .add_from_buf  = add_from_buf,
};

// ----------------------------------------------------------------------
// psc_bnd_c_setup

static void
psc_bnd_c_setup(struct mrc_obj *obj)
{
  struct psc_bnd *bnd = to_psc_bnd(obj);
  struct psc_bnd_c *bnd_c = to_psc_bnd_c(bnd);

  bnd_c->ddc = mrc_domain_create_ddc(psc.mrc_domain);
  mrc_ddc_set_funcs(bnd_c->ddc, &ddc_funcs);
  mrc_ddc_set_param_int3(bnd_c->ddc, "ibn", psc.ibn);
  mrc_ddc_set_param_int(bnd_c->ddc, "max_n_fields", 6);
  mrc_ddc_set_param_int(bnd_c->ddc, "size_of_type", sizeof(fields_base_real_t));
  mrc_ddc_setup(bnd_c->ddc);

  create_bnd();
}

// ======================================================================
// psc_bnd: subclass "c"

struct psc_bnd_ops psc_bnd_c_ops = {
  .name                  = "c",
  .size                  = sizeof(struct psc_bnd_c),
  .setup                 = psc_bnd_c_setup,
};
