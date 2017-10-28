
#include "psc_bnd_cuda.h"
#include "psc_cuda.h"
#include "psc_bnd_cuda_fields.h"
#include "cuda_mfields.h"

#include <mrc_ddc.h>
#include <mrc_profile.h>

// ======================================================================
// ddc funcs

void
psc_bnd_fld_cuda_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *_ctx)
{
  struct psc_mfields *mflds = _ctx;
  struct cuda_mfields_bnd *cbnd = psc_mfields_cuda(mflds)->cbnd;
  struct cuda_mfields_bnd_patch *cf = &cbnd->bnd_by_patch[p];
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
  struct cuda_mfields_bnd *cbnd = psc_mfields_cuda(mflds)->cbnd;
  struct cuda_mfields_bnd_patch *cf = &cbnd->bnd_by_patch[p];
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
  struct cuda_mfields_bnd *cbnd = psc_mfields_cuda(mflds)->cbnd;
  struct cuda_mfields_bnd_patch *cf = &cbnd->bnd_by_patch[p];
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
// psc_bnd_fld_cuda_add_ghosts

void
psc_bnd_fld_cuda_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds_base, int mb, int me)
{
  int size;
  MPI_Comm_size(psc_bnd_comm(bnd), &size);
  if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // double periodic single patch
    struct psc_mfields *flds = psc_mfields_get_as(flds_base, "cuda", mb, me);
    cuda_add_ghosts_periodic_yz(flds, 0, mb, me);
    psc_mfields_put_as(flds, flds_base, mb, me);
  } else if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] != BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // z-periodic single patch
    struct psc_mfields *flds = psc_mfields_get_as(flds_base, "cuda", mb, me);
    cuda_add_ghosts_periodic_z(flds, 0, mb, me);
    psc_mfields_put_as(flds, flds_base, mb, me);
  } else {
    struct psc_mfields *flds_cuda = psc_mfields_get_as(flds_base, "cuda", mb, me);

    __fields_cuda_from_device_inside(flds_cuda, mb, me);
    mrc_ddc_add_ghosts(bnd->ddc, 0, me - mb, flds_cuda);
    __fields_cuda_to_device_inside(flds_cuda, mb, me);

    psc_mfields_put_as(flds_cuda, flds_base, mb, me);
  }
}

// ----------------------------------------------------------------------
// psc_bnd_fld_cuda_fill_ghosts

void
psc_bnd_fld_cuda_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds_base, int mb, int me)
{
  static int pr1, pr2, pr3, pr4, pr5;
  if (!pr1) {
    pr1 = prof_register("cuda_fill_ghosts_1", 1., 0, 0);
    pr2 = prof_register("cuda_fill_ghosts_2", 1., 0, 0);
    pr3 = prof_register("cuda_fill_ghosts_3", 1., 0, 0);
    pr4 = prof_register("cuda_fill_ghosts_4", 1., 0, 0);
    pr5 = prof_register("cuda_fill_ghosts_5", 1., 0, 0);
  }
  int size;
  MPI_Comm_size(psc_bnd_comm(bnd), &size);
  if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // double periodic single patch
    struct psc_mfields *flds = psc_mfields_get_as(flds_base, "cuda", mb, me);
    cuda_fill_ghosts_periodic_yz(flds, 0, mb, me);
    psc_mfields_put_as(flds, flds_base, mb, me);
  } else if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] != BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // z-periodic single patch
    struct psc_mfields *flds = psc_mfields_get_as(flds_base, "cuda", mb, me);
    cuda_fill_ghosts_periodic_z(flds, 0, mb, me);
    psc_mfields_put_as(flds, flds_base, mb, me);
  } else {
    struct psc_mfields *flds_cuda = psc_mfields_get_as(flds_base, "cuda", mb, me);

    EXTERN_C void __fields_cuda_fill_ghosts_setup(struct psc_mfields *mflds, struct mrc_ddc *ddc);
    __fields_cuda_fill_ghosts_setup(flds_cuda, bnd->ddc);

    EXTERN_C void __fields_cuda_from_device_inside_only(struct psc_mfields *mflds, int mb, int me);
    prof_start(pr1);
    __fields_cuda_from_device_inside(flds_cuda, mb, me); // FIXME _only
    prof_stop(pr1);

    prof_start(pr2);
    mrc_ddc_fill_ghosts_begin(bnd->ddc, 0, me - mb, flds_cuda);
    prof_stop(pr2);
    prof_start(pr3);
#if 0
    mrc_ddc_fill_ghosts_local(bnd->ddc, 0, me - mb, flds_cuda);
#endif
#if 1
    EXTERN_C void __fields_cuda_fill_ghosts_local(struct psc_mfields *mflds, int mb, int me);
    __fields_cuda_fill_ghosts_local(flds_cuda, mb, me);
#endif
    prof_stop(pr3);
    prof_start(pr4);
    mrc_ddc_fill_ghosts_end(bnd->ddc, 0, me - mb, flds_cuda);
    prof_stop(pr4);

    prof_start(pr5);
    __fields_cuda_to_device_outside(flds_cuda, mb, me);
    prof_stop(pr5);

    psc_mfields_put_as(flds_cuda, flds_base, mb, me);
  }
}

// ======================================================================
// psc_bnd: subclass "cuda"

struct psc_bnd_ops psc_bnd_cuda_ops = {
  .name                    = "cuda",
  .create_ddc              = psc_bnd_fld_cuda_create,
  .add_ghosts              = psc_bnd_fld_cuda_add_ghosts,
  .fill_ghosts             = psc_bnd_fld_cuda_fill_ghosts,
};

