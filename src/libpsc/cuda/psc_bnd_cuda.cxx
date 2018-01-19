
#include "psc_bnd_private.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"
#include "cuda_iface_bnd.h"

#include <mrc_ddc_private.h>
#include <mrc_profile.h>

struct psc_bnd_cuda {
  struct cuda_mfields_bnd *cbnd;
};

#define psc_bnd_cuda(bnd) mrc_to_subobj(bnd, struct psc_bnd_cuda)

// ======================================================================
// ddc funcs

void
psc_bnd_fld_cuda_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *_ctx)
{
  struct psc_bnd *bnd = (struct psc_bnd *) _ctx;
  struct cuda_mfields_bnd *cbnd = psc_bnd_cuda(bnd)->cbnd;
  struct cuda_mfields_bnd_patch *cf = cuda_mfields_bnd_get_patch(cbnd, p);
  fields_cuda_real_t *buf = (fields_cuda_real_t *) _buf;

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
  struct psc_bnd *bnd = (struct psc_bnd *) _ctx;
  struct cuda_mfields_bnd *cbnd = psc_bnd_cuda(bnd)->cbnd;
  struct cuda_mfields_bnd_patch *cf = cuda_mfields_bnd_get_patch(cbnd, p);
  fields_cuda_real_t *buf = (fields_cuda_real_t *) _buf;

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
  struct psc_bnd *bnd = (struct psc_bnd *) _ctx;
  struct cuda_mfields_bnd *cbnd = psc_bnd_cuda(bnd)->cbnd;
  struct cuda_mfields_bnd_patch *cf = cuda_mfields_bnd_get_patch(cbnd, p);
  fields_cuda_real_t *buf = (fields_cuda_real_t *) _buf;

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
// psc_bnd_cuda_setup

static void
psc_bnd_cuda_setup(struct psc_bnd *bnd)
{
  struct psc_bnd_cuda *sub = psc_bnd_cuda(bnd);

  psc_bnd_setup_super(bnd); // this calls back ::create_ddc() to set up ::ddc

  struct cuda_mfields_bnd_params prm;

  struct psc *psc = bnd->psc;
  prm.n_patches = psc->nr_patches;
  assert(prm.n_patches > 0);
  // FIXME, it'd be nicer if the interface was ldims / ibn based
  for (int d = 0; d < 3; d++) {
    prm.ib[d] = -psc->ibn[d];
    prm.im[d] = psc->patch[0].ldims[d] + 2 * psc->ibn[d];
  }
  
  struct mrc_ddc_multi *multi = mrc_ddc_multi(bnd->ddc);
  struct mrc_ddc_pattern2 *patt2 = &multi->fill_ghosts2;
  struct mrc_ddc_rank_info *ri = patt2->ri;
    
  prm.n_recv_entries = ri[multi->mpi_rank].n_recv_entries;
  prm.recv_entry = (struct cuda_mfields_bnd_entry *) calloc(prm.n_recv_entries, sizeof(*prm.recv_entry));

  for (int i = 0; i < prm.n_recv_entries; i++) {
    struct mrc_ddc_sendrecv_entry *re = &ri[multi->mpi_rank].recv_entry[i];
    prm.recv_entry[i].patch     = re->patch;
    prm.recv_entry[i].nei_patch = re->nei_patch;
    prm.recv_entry[i].dir1      = re->dir1;
    //mprintf("i %d patch %d dir1 %d nei_patch %d\n", i, re->patch, re->dir1, re->nei_patch);
  }

  sub->cbnd = cuda_mfields_bnd_create();
  cuda_mfields_bnd_ctor(sub->cbnd, &prm);

  free(prm.recv_entry);
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_destroy

static void
psc_bnd_cuda_destroy(struct psc_bnd *bnd)
{
  struct psc_bnd_cuda *sub = psc_bnd_cuda(bnd);

  cuda_mfields_bnd_dtor(sub->cbnd);
  cuda_mfields_bnd_destroy(sub->cbnd);
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_add_ghosts

void
psc_bnd_cuda_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds_base, int mb, int me)
{
  struct cuda_mfields_bnd *cbnd = psc_bnd_cuda(bnd)->cbnd;
  mfields_cuda_t mf = mflds_base->get_as<mfields_cuda_t>(mb, me);

  int size;
  MPI_Comm_size(psc_bnd_comm(bnd), &size);
  if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // double periodic single patch
    cuda_add_ghosts_periodic_yz(mf->cmflds, 0, mb, me);
  } else if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] != BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // z-periodic single patch
    cuda_add_ghosts_periodic_z(mf->cmflds, 0, mb, me);
  } else {
    cuda_mfields_bnd_from_device_inside(cbnd, mf->cmflds, mb, me);
    mrc_ddc_add_ghosts(bnd->ddc, 0, me - mb, bnd);
    cuda_mfields_bnd_to_device_inside(cbnd, mf->cmflds, mb, me);
  }

  mf.put_as(mflds_base, mb, me);
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_fill_ghosts

void
psc_bnd_cuda_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds_base, int mb, int me)
{
  struct cuda_mfields_bnd *cbnd = psc_bnd_cuda(bnd)->cbnd;

  static int pr1, pr2, pr3, pr4, pr5;
  if (!pr1) {
    pr1 = prof_register("cuda_fill_ghosts_1", 1., 0, 0);
    pr2 = prof_register("cuda_fill_ghosts_2", 1., 0, 0);
    pr3 = prof_register("cuda_fill_ghosts_3", 1., 0, 0);
    pr4 = prof_register("cuda_fill_ghosts_4", 1., 0, 0);
    pr5 = prof_register("cuda_fill_ghosts_5", 1., 0, 0);
  }

  mfields_cuda_t mf = mflds_base->get_as<mfields_cuda_t>(mb, me);

  int size;
  MPI_Comm_size(psc_bnd_comm(bnd), &size);

  if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // double periodic single patch
    cuda_fill_ghosts_periodic_yz(mf->cmflds, 0, mb, me);
  } else if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] != BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // z-periodic single patch
    cuda_fill_ghosts_periodic_z(mf->cmflds, 0, mb, me);
  } else {
    prof_start(pr1);
    cuda_mfields_bnd_from_device_inside(cbnd, mf->cmflds, mb, me); // FIXME _only
    prof_stop(pr1);

    prof_start(pr2);
    mrc_ddc_fill_ghosts_begin(bnd->ddc, 0, me - mb, bnd);
    prof_stop(pr2);
    prof_start(pr3);
#if 0
    mrc_ddc_fill_ghosts_local(bnd->ddc, 0, me - mb, mflds);
#endif
#if 1
    cuda_mfields_bnd_fill_ghosts_local(cbnd, mf->cmflds, mb, me);
#endif
    prof_stop(pr3);
    prof_start(pr4);
    mrc_ddc_fill_ghosts_end(bnd->ddc, 0, me - mb, mf.mflds());
    prof_stop(pr4);

    prof_start(pr5);
    cuda_mfields_bnd_to_device_outside(cbnd, mf->cmflds, mb, me);
    prof_stop(pr5);
  }

  mf.put_as(mflds_base, mb, me);
}

// ======================================================================
// psc_bnd: subclass "cuda"

struct psc_bnd_ops_cuda : psc_bnd_ops {
  psc_bnd_ops_cuda() {
    name                    = "cuda";
    size                    = sizeof(struct psc_bnd_cuda);
    setup                   = psc_bnd_cuda_setup;
    destroy                 = psc_bnd_cuda_destroy;
    create_ddc              = psc_bnd_cuda_create_ddc;
    add_ghosts              = psc_bnd_cuda_add_ghosts;
    fill_ghosts             = psc_bnd_cuda_fill_ghosts;
  }
} psc_bnd_cuda_ops;

