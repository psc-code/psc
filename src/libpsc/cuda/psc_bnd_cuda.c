
#include "psc_bnd_cuda.h"
#include "psc_fields_cuda.h"
#include "psc_bnd_cuda_fields.h"
#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "cuda_mfields.h"

#include <mrc_ddc_private.h>
#include <mrc_profile.h>

struct psc_bnd_cuda {
  struct cuda_mfields_bnd *cbnd;
};

#define psc_bnd_cuda(bnd) mrc_to_subobj(bnd, struct psc_bnd_cuda)

static inline struct cuda_mfields_bnd *
psc_bnd_cuda_cbnd(struct psc_bnd *bnd)
{
  struct psc_bnd_cuda *sub = psc_bnd_cuda(bnd);
  assert(sub->cbnd);
  return sub->cbnd;
}

// ======================================================================
// ddc funcs

void
psc_bnd_fld_cuda_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *_ctx)
{
  struct psc_bnd *bnd = _ctx;
  struct cuda_mfields_bnd *cbnd = psc_bnd_cuda_cbnd(bnd);
  struct cuda_mfields_bnd_patch *cf = cuda_mfields_bnd_get_patch(cbnd, p);
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
  struct psc_bnd *bnd = _ctx;
  struct cuda_mfields_bnd *cbnd = psc_bnd_cuda_cbnd(bnd);
  struct cuda_mfields_bnd_patch *cf = cuda_mfields_bnd_get_patch(cbnd, p);
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
  struct psc_bnd *bnd = _ctx;
  struct cuda_mfields_bnd *cbnd = psc_bnd_cuda_cbnd(bnd);
  struct cuda_mfields_bnd_patch *cf = cuda_mfields_bnd_get_patch(cbnd, p);
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
// psc_bnd_cuda_create_ddc

void
psc_bnd_cuda_create_ddc(struct psc_bnd *bnd)
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
// psc_bnd_cuda_fill_ghosts_setup

struct cuda_mfields_bnd_entry {
  int patch;
  int nei_patch;
  int dir1;
};

static void
psc_bnd_cuda_fill_ghosts_setup(struct psc_bnd *bnd)
{
  struct cuda_mfields_bnd *cbnd = psc_bnd_cuda_cbnd(bnd);

  int n_entries;
  struct cuda_mfields_bnd_entry *entries;
  {
    struct mrc_ddc_multi *multi = mrc_ddc_multi(bnd->ddc);
    struct mrc_ddc_pattern2 *patt2 = &multi->fill_ghosts2;
    struct mrc_ddc_rank_info *ri = patt2->ri;
    
    n_entries = ri[multi->mpi_rank].n_recv_entries;
    entries = calloc(n_entries, sizeof(*entries));
    for (int i = 0; i < n_entries; i++) {
      struct mrc_ddc_sendrecv_entry *re = &ri[multi->mpi_rank].recv_entry[i];
      entries[i].patch = re->patch;
      entries[i].nei_patch = re->nei_patch;
      entries[i].dir1 = re->dir1;
      mprintf("i %d patch %d dir1 %d nei_patch %d\n", i, re->patch, re->dir1, re->nei_patch);
    }
  }
  
  cbnd->h_nei_patch = calloc(9 * cbnd->n_patches, sizeof(*cbnd->h_nei_patch));

  for (int p = 0; p < cbnd->n_patches; p++) {
    for (int dir1 = 0; dir1 < 9; dir1++) {
      cbnd->h_nei_patch[p * 9 + dir1] = -1;
    }
  }
  for (int i = 0; i < n_entries; i++) {
    cbnd->h_nei_patch[entries[i].patch * 9 + entries[i].dir1 / 3] = entries[i].nei_patch;
  }

  free(entries);

  cuda_mfields_bnd_setup_d_nei_patch(cbnd);
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_setup

static void
psc_bnd_cuda_setup(struct psc_bnd *bnd)
{
  psc_bnd_setup_super(bnd); // this calls back to set up ::ddc
  // FIXME then, we should to the above cbnd setup, but we don't have cbnd/cmflds
  // at this point... yet
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_add_ghosts

void
psc_bnd_cuda_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds_base, int mb, int me)
{
  struct cuda_mfields_bnd *cbnd = psc_bnd_cuda_cbnd(bnd);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "cuda", mb, me);
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds)->cmflds;

  int size;
  MPI_Comm_size(psc_bnd_comm(bnd), &size);
  if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // double periodic single patch
    cuda_add_ghosts_periodic_yz(cmflds, 0, mb, me);
  } else if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] != BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // z-periodic single patch
    cuda_add_ghosts_periodic_z(cmflds, 0, mb, me);
  } else {
    __fields_cuda_from_device_inside(cbnd, cmflds, mb, me);
    mrc_ddc_add_ghosts(bnd->ddc, 0, me - mb, bnd);
    __fields_cuda_to_device_inside(cbnd, cmflds, mb, me);
  }

  psc_mfields_put_as(mflds, mflds_base, mb, me);
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_fill_ghosts

void
psc_bnd_cuda_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds_base, int mb, int me)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "cuda", mb, me);
  struct cuda_mfields *cmflds = psc_mfields_cuda(mflds)->cmflds;

  struct psc_bnd_cuda *sub = psc_bnd_cuda(bnd);
  if (!sub->cbnd) {
    sub->cbnd = psc_mfields_cuda(mflds)->cbnd;
    psc_bnd_cuda_fill_ghosts_setup(bnd);
  }

  struct cuda_mfields_bnd *cbnd = psc_bnd_cuda_cbnd(bnd);

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
    cuda_fill_ghosts_periodic_yz(cmflds, 0, mb, me);
  } else if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] != BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // z-periodic single patch
    cuda_fill_ghosts_periodic_z(cmflds, 0, mb, me);
  } else {
    prof_start(pr1);
    __fields_cuda_from_device_inside(cbnd, cmflds, mb, me); // FIXME _only
    prof_stop(pr1);

    prof_start(pr2);
    mrc_ddc_fill_ghosts_begin(bnd->ddc, 0, me - mb, bnd);
    prof_stop(pr2);
    prof_start(pr3);
#if 0
    mrc_ddc_fill_ghosts_local(bnd->ddc, 0, me - mb, mflds);
#endif
#if 1
    __fields_cuda_fill_ghosts_local(cbnd, cmflds, mb, me);
#endif
    prof_stop(pr3);
    prof_start(pr4);
    mrc_ddc_fill_ghosts_end(bnd->ddc, 0, me - mb, mflds);
    prof_stop(pr4);

    prof_start(pr5);
    __fields_cuda_to_device_outside(cbnd, cmflds, mb, me);
    prof_stop(pr5);
  }

  psc_mfields_put_as(mflds, mflds_base, mb, me);
}

// ======================================================================
// psc_bnd: subclass "cuda"

struct psc_bnd_ops psc_bnd_cuda_ops = {
  .name                    = "cuda",
  .size                    = sizeof(struct psc_bnd_cuda),
  .setup                   = psc_bnd_cuda_setup,
  .create_ddc              = psc_bnd_cuda_create_ddc,
  .add_ghosts              = psc_bnd_cuda_add_ghosts,
  .fill_ghosts             = psc_bnd_cuda_fill_ghosts,
};

