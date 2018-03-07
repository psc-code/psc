
#include "psc_bnd_private.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "bnd.hxx"

#include <mrc_ddc_private.h>
#include <mrc_profile.h>

void psc_bnd_fld_cuda_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *_ctx);
void psc_bnd_fld_cuda_copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *_ctx);
void psc_bnd_fld_cuda_add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *_ctx);

static struct mrc_ddc_funcs ddc_funcs = {
  .copy_to_buf   = psc_bnd_fld_cuda_copy_to_buf,
  .copy_from_buf = psc_bnd_fld_cuda_copy_from_buf,
  .add_from_buf  = psc_bnd_fld_cuda_add_from_buf,
};

struct BndCuda : BndBase {
  struct cuda_mfields_bnd *cbnd;

  // ----------------------------------------------------------------------
  // ctor

  BndCuda(const Grid_t& grid, mrc_domain *domain, int ibn[3])
  {
    ddc_ = mrc_domain_create_ddc(domain);
    mrc_ddc_set_funcs(ddc_, &ddc_funcs);
    mrc_ddc_set_param_int3(ddc_, "ibn", ibn);
    mrc_ddc_set_param_int(ddc_, "max_n_fields", 6);
    mrc_ddc_set_param_int(ddc_, "size_of_type", sizeof(fields_cuda_real_t));
    mrc_ddc_setup(ddc_);

    struct cuda_mfields_bnd_params prm;

    prm.n_patches = grid.n_patches();
    assert(prm.n_patches > 0);
    // FIXME, it'd be nicer if the interface was ldims / ibn based
    for (int d = 0; d < 3; d++) {
      prm.ib[d] = -ibn[d];
      prm.im[d] = grid.ldims[d] + 2 * ibn[d];
    }
    
    struct mrc_ddc_multi *multi = mrc_ddc_multi(ddc_);
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
    
    cbnd = cuda_mfields_bnd_create();
    cuda_mfields_bnd_ctor(cbnd, &prm);
    
    free(prm.recv_entry);
  }

  // ----------------------------------------------------------------------
  // dtor

  ~BndCuda()
  {
    cuda_mfields_bnd_dtor(cbnd);
    cuda_mfields_bnd_destroy(cbnd);

    mrc_ddc_destroy(ddc_);
  }

  // ----------------------------------------------------------------------
  // reset
  
  void reset()
  {
    assert(0);
    // mrc_ddc_destroy(bnd_->ddc);
    // ops->create_ddc(bnd_);
  }

  //private:
  mrc_ddc* ddc_;
};

// ======================================================================
// ddc funcs

void
psc_bnd_fld_cuda_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *_buf, void *_ctx)
{
  struct psc_bnd *_bnd = (struct psc_bnd *) _ctx;
  PscBnd<BndCuda> bnd(_bnd);
  struct cuda_mfields_bnd *cbnd = bnd->cbnd;
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
  struct psc_bnd *_bnd = (struct psc_bnd *) _ctx;
  PscBnd<BndCuda> bnd(_bnd);
  struct cuda_mfields_bnd *cbnd = bnd->cbnd;
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
  struct psc_bnd *_bnd = (struct psc_bnd *) _ctx;
  PscBnd<BndCuda> bnd(_bnd);
  struct cuda_mfields_bnd *cbnd = bnd->cbnd;
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

// ----------------------------------------------------------------------
// psc_bnd_cuda_add_ghosts

void
psc_bnd_cuda_add_ghosts(struct psc_bnd *_bnd, struct psc_mfields *mflds_base, int mb, int me)
{
  PscBnd<BndCuda> bnd(_bnd);
  PscMfieldsCuda mf = mflds_base->get_as<PscMfieldsCuda>(mb, me);

  int size;
  MPI_Comm_size(psc_bnd_comm(_bnd), &size);
  if (size == 1 && ppsc->n_patches() == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // double periodic single patch
    cuda_add_ghosts_periodic_yz(mf->cmflds, 0, mb, me);
  } else if (size == 1 && ppsc->n_patches() == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] != BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // z-periodic single patch
    cuda_add_ghosts_periodic_z(mf->cmflds, 0, mb, me);
  } else {
    cuda_mfields_bnd_from_device_inside(bnd->cbnd, mf->cmflds, mb, me);
    mrc_ddc_add_ghosts(bnd->ddc_, 0, me - mb, _bnd);
    cuda_mfields_bnd_to_device_inside(bnd->cbnd, mf->cmflds, mb, me);
  }

  mf.put_as(mflds_base, mb, me);
}

// ----------------------------------------------------------------------
// psc_bnd_cuda_fill_ghosts

void
psc_bnd_cuda_fill_ghosts(struct psc_bnd *_bnd, struct psc_mfields *mflds_base, int mb, int me)
{
  PscBnd<BndCuda> bnd(_bnd);

  static int pr1, pr2, pr3, pr4, pr5;
  if (!pr1) {
    pr1 = prof_register("cuda_fill_ghosts_1", 1., 0, 0);
    pr2 = prof_register("cuda_fill_ghosts_2", 1., 0, 0);
    pr3 = prof_register("cuda_fill_ghosts_3", 1., 0, 0);
    pr4 = prof_register("cuda_fill_ghosts_4", 1., 0, 0);
    pr5 = prof_register("cuda_fill_ghosts_5", 1., 0, 0);
  }

  PscMfieldsCuda mf = mflds_base->get_as<PscMfieldsCuda>(mb, me);

  int size;
  MPI_Comm_size(psc_bnd_comm(_bnd), &size);

  if (size == 1 && ppsc->n_patches() == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // double periodic single patch
    cuda_fill_ghosts_periodic_yz(mf->cmflds, 0, mb, me);
  } else if (size == 1 && ppsc->n_patches() == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] != BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    // z-periodic single patch
    cuda_fill_ghosts_periodic_z(mf->cmflds, 0, mb, me);
  } else {
    prof_start(pr1);
    cuda_mfields_bnd_from_device_inside(bnd->cbnd, mf->cmflds, mb, me); // FIXME _only
    prof_stop(pr1);

    prof_start(pr2);
    mrc_ddc_fill_ghosts_begin(bnd->ddc_, 0, me - mb, _bnd);
    prof_stop(pr2);
    prof_start(pr3);
#if 0
    mrc_ddc_fill_ghosts_local(bnd->ddc_, 0, me - mb, mflds);
#endif
#if 1
    cuda_mfields_bnd_fill_ghosts_local(bnd->cbnd, mf->cmflds, mb, me);
#endif
    prof_stop(pr3);
    prof_start(pr4);
    mrc_ddc_fill_ghosts_end(bnd->ddc_, 0, me - mb, mf.mflds());
    prof_stop(pr4);

    prof_start(pr5);
    cuda_mfields_bnd_to_device_outside(bnd->cbnd, mf->cmflds, mb, me);
    prof_stop(pr5);
  }

  mf.put_as(mflds_base, mb, me);
}

// ======================================================================
// psc_bnd: subclass "cuda"

struct psc_bnd_ops_cuda : psc_bnd_ops {
  using PscBnd_t = PscBndWrapper<BndCuda>;
  psc_bnd_ops_cuda() {
    name                    = "cuda";
    size                    = PscBnd_t::size;
    reset                   = PscBnd_t::reset;
    setup                   = PscBnd_t::setup;
    destroy                 = PscBnd_t::destroy;
    add_ghosts              = psc_bnd_cuda_add_ghosts;
    fill_ghosts             = psc_bnd_cuda_fill_ghosts;
  }
} psc_bnd_cuda_ops;

