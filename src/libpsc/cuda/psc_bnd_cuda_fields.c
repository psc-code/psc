
#include "psc_bnd_cuda_fields.h"

#include "psc_cuda.h"
#include "psc_bnd_private.h"
#include "../psc_bnd/psc_bnd_fld.h"

#include <mrc_ddc.h>
#include <mrc_profile.h>

static void
psc_fields_cuda_bnd_prep(struct psc_fields *pf, int nr_fields)
{
  // FIXME, could be done once and cached
  struct psc_fields_cuda_bnd *cf = &psc_fields_cuda(pf)->bnd;
  int sz = 1;
  for (int d = 0; d < 3; d++) {
    if (pf->im[d] == 1 - 2 * pf->ib[d]) { // only 1 non-ghost point
      cf->im[d] = 1;
      cf->ib[d] = 0;
    } else {
      cf->im[d] = pf->im[d];
      cf->ib[d] = pf->ib[d];
    }
    sz *= cf->im[d];
  }
  cf->arr = malloc(nr_fields * sz * sizeof(*cf->arr));
  cf->arr_off = cf->arr 
    - ((cf->ib[2] * cf->im[1] + cf->ib[1]) * cf->im[0] + cf->ib[0]);
}

static void
psc_fields_cuda_bnd_post(struct psc_fields *pf)
{
  struct psc_fields_cuda_bnd *cf = &psc_fields_cuda(pf)->bnd;
  free(cf->arr);
}

// ======================================================================
// ddc funcs

void
psc_bnd_fld_cuda_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *_ctx)
{
  struct psc_mfields *mflds = _ctx;
  struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
  struct psc_fields_cuda_bnd *cf = &psc_fields_cuda(flds)->bnd;
  fields_cuda_real_t *buf = _buf;

  me -= mb; // FIXME, the "mix" bnd needs this adjustment
  mb = 0;
  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_DDC_BUF3(buf, m - mb, ix,iy,iz) = F3_CF(cf, m, ix,iy,iz);
	}
      }
    }
  }
}

void
psc_bnd_fld_cuda_add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *_ctx)
{
  struct psc_mfields *mflds = _ctx;
  struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
  struct psc_fields_cuda_bnd *cf = &psc_fields_cuda(flds)->bnd;
  fields_cuda_real_t *buf = _buf;

  me -= mb;
  mb = 0;
  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  F3_CF(cf, m, ix,iy,iz) += MRC_DDC_BUF3(buf, m - mb, ix,iy,iz);
	}
      }
    }
  }
}

void
psc_bnd_fld_cuda_copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *_ctx)
{
  struct psc_mfields *mflds = _ctx;
  struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
  struct psc_fields_cuda_bnd *cf = &psc_fields_cuda(flds)->bnd;
  fields_cuda_real_t *buf = _buf;

  me -= mb;
  mb = 0;
  for (int m = mb; m < me; m++) {
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  F3_CF(cf, m, ix,iy,iz) = MRC_DDC_BUF3(buf, m - mb, ix,iy,iz);
	}
      }
    }
  }
}

static struct mrc_ddc_funcs ddc_funcs = {
  .copy_to_buf   = psc_bnd_fld_cuda_copy_to_buf,
  .copy_from_buf = psc_bnd_fld_cuda_copy_from_buf,
  .add_from_buf  = psc_bnd_fld_cuda_add_from_buf,
};

// ----------------------------------------------------------------------
// psc_bnd_fld_cuda_create

void
psc_bnd_fld_cuda_create(struct psc_bnd *bnd)
{
  struct psc *psc = bnd->psc;

  bnd->ddc = mrc_domain_create_ddc(psc->mrc_domain);
  mrc_ddc_set_funcs(bnd->ddc, &ddc_funcs);
  mrc_ddc_set_param_int3(bnd->ddc, "ibn", psc->ibn);
  mrc_ddc_set_param_int(bnd->ddc, "max_n_fields", 6);
  mrc_ddc_set_param_int(bnd->ddc, "size_of_type", sizeof(fields_cuda_real_t));
  mrc_ddc_setup(bnd->ddc);
}

// ----------------------------------------------------------------------

void
psc_bnd_fld_cuda_add_ghosts_prep(struct psc_bnd *bnd, struct psc_fields *pf, int mb, int me)
{
  psc_fields_cuda_bnd_prep(pf, me - mb);
  __fields_cuda_from_device_inside(pf, mb, me);
}

void
psc_bnd_fld_cuda_add_ghosts_post(struct psc_bnd *bnd, struct psc_fields *pf, int mb, int me)
{
  __fields_cuda_to_device_inside(pf, mb, me);
  psc_fields_cuda_bnd_post(pf);
}

// ----------------------------------------------------------------------

void
psc_bnd_fld_cuda_fill_ghosts_prep(struct psc_bnd *bnd, struct psc_fields *pf, int mb, int me)
{
  psc_fields_cuda_bnd_prep(pf, me - mb);
  __fields_cuda_from_device_inside(pf, mb, me);
}

void
psc_bnd_fld_cuda_fill_ghosts_post(struct psc_bnd *bnd, struct psc_fields *pf, int mb, int me)
{
  __fields_cuda_to_device_outside(pf, mb, me);
  psc_fields_cuda_bnd_post(pf);
}

// ----------------------------------------------------------------------
// psc_bnd_fld_cuda_add_ghosts

void
psc_bnd_fld_cuda_add_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base, int mb, int me)
{
  int size;
  MPI_Comm_size(psc_bnd_comm(bnd), &size);
  if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // double periodic single patch
    mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, mb, me);
    cuda_add_ghosts_periodic_yz(0, psc_mfields_get_patch(flds, 0), mb, me);
    psc_mfields_put_cuda(flds, flds_base, mb, me);
  } else if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] != BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // z-periodic single patch
    mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, mb, me);
    cuda_add_ghosts_periodic_z(0, psc_mfields_get_patch(flds, 0), mb, me);
    psc_mfields_put_cuda(flds, flds_base, mb, me);
  } else {
    mfields_cuda_t *flds_cuda = psc_mfields_get_cuda(flds_base, mb, me);

    psc_foreach_patch(ppsc, p) {
      psc_bnd_fld_cuda_add_ghosts_prep(bnd, psc_mfields_get_patch(flds_cuda, p), mb, me);
    }

    mrc_ddc_add_ghosts(bnd->ddc, 0, me - mb, flds_cuda);

    psc_foreach_patch(ppsc, p) {
      psc_bnd_fld_cuda_add_ghosts_post(bnd, psc_mfields_get_patch(flds_cuda, p), mb, me);
    }

    psc_mfields_put_cuda(flds_cuda, flds_base, mb, me);
  }
}

// ----------------------------------------------------------------------
// psc_bnd_fld_cuda_fill_ghosts

void
psc_bnd_fld_cuda_fill_ghosts(struct psc_bnd *bnd, mfields_base_t *flds_base, int mb, int me)
{
  static int pr1, pr3, pr5;
  if (!pr1) {
    pr1 = prof_register("cuda_fill_ghosts_1", 1., 0, 0);
    pr3 = prof_register("cuda_fill_ghosts_3", 1., 0, 0);
    pr5 = prof_register("cuda_fill_ghosts_5", 1., 0, 0);
  }
  int size;
  MPI_Comm_size(psc_bnd_comm(bnd), &size);
  if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // double periodic single patch
    mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, mb, me);
    cuda_fill_ghosts_periodic_yz(0, psc_mfields_get_patch(flds, 0), mb, me);
    psc_mfields_put_cuda(flds, flds_base, mb, me);
  } else if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] != BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // z-periodic single patch
    mfields_cuda_t *flds = psc_mfields_get_cuda(flds_base, mb, me);
    cuda_fill_ghosts_periodic_z(0, psc_mfields_get_patch(flds, 0), mb, me);
    psc_mfields_put_cuda(flds, flds_base, mb, me);
  } else {
    mfields_cuda_t *flds_cuda = psc_mfields_get_cuda(flds_base, mb, me);

    prof_start(pr1);
    psc_foreach_patch(ppsc, p) {
      psc_bnd_fld_cuda_fill_ghosts_prep(bnd, psc_mfields_get_patch(flds_cuda, p), mb, me);
    }
    prof_stop(pr1);

    prof_start(pr3);
    mrc_ddc_fill_ghosts(bnd->ddc, 0, me - mb, flds_cuda);
    prof_stop(pr3);

    prof_start(pr5);
    psc_foreach_patch(ppsc, p) {
      psc_bnd_fld_cuda_fill_ghosts_post(bnd, psc_mfields_get_patch(flds_cuda, p), mb, me);
    }
    prof_stop(pr5);

    psc_mfields_put_cuda(flds_cuda, flds_base, mb, me);
  }
}

