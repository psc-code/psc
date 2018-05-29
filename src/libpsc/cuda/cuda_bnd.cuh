
#pragma once

#include "psc_fields_cuda.h"
#include "psc_fields_single.h"
#include "fields.hxx"

#include "mrc_ddc_private.h"

#include <thrust/gather.h>

#define mrc_ddc_multi(ddc) mrc_to_subobj(ddc, struct mrc_ddc_multi)

static void
mrc_ddc_multi_set_mpi_type(struct mrc_ddc *ddc)
{
  if (ddc->size_of_type == sizeof(float)) {
    ddc->mpi_type = MPI_FLOAT;
  } else if (ddc->size_of_type == sizeof(double)) {
    ddc->mpi_type = MPI_DOUBLE;
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_free_buffers

static void
mrc_ddc_multi_free_buffers(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2)
{
  free(patt2->send_buf);
  free(patt2->recv_buf);
  free(patt2->local_buf);
}

// ----------------------------------------------------------------------
// mrc_ddc_multi_alloc_buffers

static void
mrc_ddc_multi_alloc_buffers(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2,
			    int n_fields)
{
  if (ddc->size_of_type > patt2->max_size_of_type ||
      n_fields > patt2->max_n_fields) {

    if (ddc->size_of_type > patt2->max_size_of_type) {
      patt2->max_size_of_type = ddc->size_of_type;
    }
    if (n_fields > patt2->max_n_fields) {
      patt2->max_n_fields = n_fields;
    }

    mrc_ddc_multi_free_buffers(ddc, patt2);

    patt2->recv_buf = malloc(patt2->n_recv * patt2->max_n_fields * patt2->max_size_of_type);
    patt2->send_buf = malloc(patt2->n_send * patt2->max_n_fields * patt2->max_size_of_type);
    patt2->local_buf = malloc(patt2->local_buf_size * patt2->max_n_fields * patt2->max_size_of_type);
  }
}

// ----------------------------------------------------------------------
// ddc_run_begin

static void
ddc_run_begin(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2,
	      int mb, int me, void *ctx,
	      void (*to_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx))
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);
  struct mrc_ddc_rank_info *ri = patt2->ri;

  // communicate aggregated buffers
  // post receives
  patt2->recv_cnt = 0;
  char* p = (char*) patt2->recv_buf;
  for (int r = 0; r < sub->mpi_size; r++) {
    if (r != sub->mpi_rank && ri[r].n_recv_entries) {
      MPI_Irecv(p, ri[r].n_recv * (me - mb), ddc->mpi_type,
		r, 0, ddc->obj.comm, &patt2->recv_req[patt2->recv_cnt++]);
      p += ri[r].n_recv * (me - mb) * ddc->size_of_type;
    }
  }  
  assert(p == (char*)patt2->recv_buf + patt2->n_recv * (me - mb) * ddc->size_of_type);

  // post sends
  patt2->send_cnt = 0;
  p = (char*) patt2->send_buf;
  for (int r = 0; r < sub->mpi_size; r++) {
    if (r != sub->mpi_rank && ri[r].n_send_entries) {
      void *p0 = p;
      for (int i = 0; i < ri[r].n_send_entries; i++) {
	struct mrc_ddc_sendrecv_entry *se = &ri[r].send_entry[i];
	to_buf(mb, me, se->patch, se->ilo, se->ihi, p, ctx);
	p += se->len * (me - mb) * ddc->size_of_type;
      }
      MPI_Isend(p0, ri[r].n_send * (me - mb), ddc->mpi_type,
		r, 0, ddc->obj.comm, &patt2->send_req[patt2->send_cnt++]);
    }
  }  
  assert(p == (char*) patt2->send_buf + patt2->n_send * (me - mb) * ddc->size_of_type);
}

// ----------------------------------------------------------------------
// ddc_run_end

static void
ddc_run_end(struct mrc_ddc *ddc, struct mrc_ddc_pattern2 *patt2,
	    int mb, int me, void *ctx,
	    void (*from_buf)(int mb, int me, int p, int ilo[3], int ihi[3], void *buf, void *ctx))
{
  struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);
  struct mrc_ddc_rank_info *ri = patt2->ri;

  MPI_Waitall(patt2->recv_cnt, patt2->recv_req, MPI_STATUSES_IGNORE);

  char* p = (char*) patt2->recv_buf;
  for (int r = 0; r < sub->mpi_size; r++) {
    if (r != sub->mpi_rank) {
      for (int i = 0; i < ri[r].n_recv_entries; i++) {
	struct mrc_ddc_sendrecv_entry *re = &ri[r].recv_entry[i];
	from_buf(mb, me, re->patch, re->ilo, re->ihi, p, ctx);
	p += re->len * (me - mb) * ddc->size_of_type;
      }
    }
  }

  MPI_Waitall(patt2->send_cnt, patt2->send_req, MPI_STATUSES_IGNORE);
}

// ======================================================================
// CudaBnd

struct CudaBnd
{
  using Mfields = MfieldsCuda;
  using fields_t = typename Mfields::fields_t;
  using real_t = typename Mfields::real_t;
  using Fields = Fields3d<fields_t>;

  // ----------------------------------------------------------------------
  // ctor
  
  CudaBnd(const Grid_t& grid, mrc_domain* domain, int ibn[3])
  {
    static struct mrc_ddc_funcs ddc_funcs;
    ddc_funcs.copy_to_buf   = copy_to_buf;
    ddc_funcs.copy_from_buf = copy_from_buf;
    ddc_funcs.add_from_buf  = add_from_buf;

    ddc_ = mrc_domain_create_ddc(domain);
    mrc_ddc_set_funcs(ddc_, &ddc_funcs);
    mrc_ddc_set_param_int3(ddc_, "ibn", ibn);
    mrc_ddc_set_param_int(ddc_, "max_n_fields", 24);
    mrc_ddc_set_param_int(ddc_, "size_of_type", sizeof(real_t));
    mrc_ddc_setup(ddc_);
  }

  // ----------------------------------------------------------------------
  // dtor
  
  ~CudaBnd()
  {
    mrc_ddc_destroy(ddc_);
  }

  // ----------------------------------------------------------------------
  // add_ghosts
  
  void add_ghosts(Mfields& mflds, int mb, int me)
  {
    auto& mflds_single = mflds.get_as<MfieldsSingle>(mb, me);
    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc_);
    
    mrc_ddc_multi_set_mpi_type(ddc_);
    mrc_ddc_multi_alloc_buffers(ddc_, &sub->add_ghosts2, me - mb);
    ddc_run_begin(ddc_, &sub->add_ghosts2, mb, me, &mflds_single, copy_to_buf);
    add_local(&sub->add_ghosts2, mb, me, mflds_single);
    ddc_run_end(ddc_, &sub->add_ghosts2, mb, me, &mflds_single, add_from_buf);
    mflds.put_as(mflds_single, mb, me);
  }
  
  void add_local(struct mrc_ddc_pattern2 *patt2,
		 int mb, int me, MfieldsSingle& mflds)
  {
    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc_);
    struct mrc_ddc_rank_info *ri = patt2->ri;
    
    for (int i = 0; i < ri[sub->mpi_rank].n_send_entries; i++) {
      struct mrc_ddc_sendrecv_entry *se = &ri[sub->mpi_rank].send_entry[i];
      struct mrc_ddc_sendrecv_entry *re = &ri[sub->mpi_rank].recv_entry[i];
      if (se->ilo[0] == se->ihi[0] ||
	  se->ilo[1] == se->ihi[1] ||
	  se->ilo[2] == se->ihi[2]) { // FIXME, we shouldn't even create these
	continue;
      }
      copy_to_buf(mb, me, se->patch, se->ilo, se->ihi, (real_t*) patt2->local_buf, mflds);
      add_from_buf(mb, me, se->nei_patch, re->ilo, re->ihi, patt2->local_buf, &mflds);
    }
  }

  // ----------------------------------------------------------------------
  // fill_ghosts

  void fill_ghosts(Mfields& mflds, int mb, int me)
  {
    // FIXME
    // I don't think we need as many points, and only stencil star
    // rather then box
    auto& mflds_single = mflds.get_as<MfieldsSingle>(mb, me);

    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc_);
    
    mrc_ddc_multi_set_mpi_type(ddc_);
    mrc_ddc_multi_alloc_buffers(ddc_, &sub->fill_ghosts2, me - mb);
    ddc_run_begin(ddc_, &sub->fill_ghosts2, mb, me, &mflds_single, copy_to_buf);
    fill_local(&sub->fill_ghosts2, mb, me, mflds_single);
    ddc_run_end(ddc_, &sub->fill_ghosts2, mb, me, &mflds_single, copy_from_buf);

    mflds.put_as(mflds_single, mb, me);
  }

  void fill_local(struct mrc_ddc_pattern2 *patt2, int mb, int me, MfieldsSingle& mflds)
  {
    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc_);
    struct mrc_ddc_rank_info *ri = patt2->ri;
    
    // overlap: local exchange
    for (int i = 0; i < ri[sub->mpi_rank].n_send_entries; i++) {
      struct mrc_ddc_sendrecv_entry *se = &ri[sub->mpi_rank].send_entry[i];
      struct mrc_ddc_sendrecv_entry *re = &ri[sub->mpi_rank].recv_entry[i];
      if (se->ilo[0] == se->ihi[0] ||
	  se->ilo[1] == se->ihi[1] ||
	  se->ilo[2] == se->ihi[2]) { // FIXME, we shouldn't even create these
	continue;
      }
      copy_to_buf(mb, me, se->patch, se->ilo, se->ihi, (real_t*) patt2->local_buf, mflds);
      copy_from_buf(mb, me, se->nei_patch, re->ilo, re->ihi, patt2->local_buf, &mflds);
    }
  }

  // ----------------------------------------------------------------------
  // copy_to_buf

  static void copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3],
			  real_t *buf, MfieldsSingle& mflds)
  {
    auto F = mflds[p];
    auto F0 = &F(mb, ilo[0], ilo[1], ilo[2]);
 
    uint size = (me - mb) * (ihi[0] - ilo[0]) * (ihi[1] - ilo[1]) * (ihi[2] - ilo[2]);
    std::vector<uint> map(size);
    auto cur = map.begin();
    for (int m = mb; m < me; m++) {
      for (int iz = ilo[2]; iz < ihi[2]; iz++) {
	for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	  for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	    *cur++ = &F(m, ix,iy,iz) - F0;
	  }
	}
      }
    }
    thrust::gather(map.begin(), map.end(), F0, buf);
    // for (auto cur : map) {
    //   *buf++ = F0[cur];
    // }
  }

  static void copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3],
			  void *_buf, void *ctx)
  {
    auto& mf = *static_cast<MfieldsSingle*>(ctx);
    auto F = mf[p];
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

  // ----------------------------------------------------------------------
  // add_from_buf

  static void add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3],
			   void *_buf, void *ctx)
  {
    auto& mf = *static_cast<MfieldsSingle*>(ctx);
    auto F = mf[p];
    real_t *buf = static_cast<real_t*>(_buf);
    
    for (int m = mb; m < me; m++) {
      for (int iz = ilo[2]; iz < ihi[2]; iz++) {
	for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	  for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	    real_t val = F(m, ix,iy,iz) + MRC_DDC_BUF3(buf, m - mb, ix,iy,iz);
	    F(m, ix,iy,iz) = val;
	  }
	}
      }
    }
  }
  
  // ----------------------------------------------------------------------
  // copy_from_buf

  static void copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3],
			    void *_buf, void *ctx)
  {
    auto& mf = *static_cast<MfieldsSingle*>(ctx);
    auto F = mf[p];
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
  mrc_ddc* ddc_;
};

