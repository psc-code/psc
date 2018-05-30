
#pragma once

#include "psc_fields_cuda.h"
#include "cuda_mfields.h"
#include "fields.hxx"

#include "mrc_ddc_private.h"

#include <thrust/gather.h>
#include <thrust/scatter.h>

using real_t = float;

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
  }
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
    cuda_mfields& cmflds = *mflds.cmflds;
    thrust::device_ptr<real_t> d_flds{cmflds.data()};
    thrust::host_vector<real_t> h_flds{d_flds, d_flds + cmflds.n_fields * cmflds.n_cells};

    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc_);
    
    mrc_ddc_multi_set_mpi_type(ddc_);
    mrc_ddc_multi_alloc_buffers(ddc_, &sub->add_ghosts2, me - mb);
    ddc_run_begin(&sub->add_ghosts2, mb, me, cmflds, h_flds);
    add_local(&sub->add_ghosts2, mb, me, h_flds, cmflds);
    ddc_run_end(&sub->add_ghosts2, mb, me, cmflds, h_flds, add_from_buf);
    thrust::copy(h_flds.begin(), h_flds.end(), d_flds);
  }
  
  // ----------------------------------------------------------------------
  // fill_ghosts

  void fill_ghosts(Mfields& mflds, int mb, int me)
  {
    cuda_mfields& cmflds = *mflds.cmflds;
    thrust::device_ptr<real_t> d_flds{cmflds.data()};
    thrust::host_vector<real_t> h_flds{d_flds, d_flds + cmflds.n_fields * cmflds.n_cells};
    // FIXME
    // I don't think we need as many points, and only stencil star
    // rather then box
    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc_);
    
    mrc_ddc_multi_set_mpi_type(ddc_);
    mrc_ddc_multi_alloc_buffers(ddc_, &sub->fill_ghosts2, me - mb);
    ddc_run_begin(&sub->fill_ghosts2, mb, me, cmflds, h_flds);
    fill_local(&sub->fill_ghosts2, mb, me, h_flds, cmflds);
    ddc_run_end(&sub->fill_ghosts2, mb, me, cmflds, h_flds, copy_from_buf);

    thrust::copy(h_flds.begin(), h_flds.end(), d_flds);
  }

  // ----------------------------------------------------------------------
  // ddc_run_begin

  void ddc_run_begin(struct mrc_ddc_pattern2 *patt2, int mb, int me,
		     cuda_mfields& cmflds, thrust::host_vector<real_t>& h_flds)
  {
    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc_);
    struct mrc_ddc_rank_info *ri = patt2->ri;

    // communicate aggregated buffers
    // post receives
    patt2->recv_cnt = 0;
    real_t* p = (real_t*) patt2->recv_buf;
    for (int r = 0; r < sub->mpi_size; r++) {
      if (r != sub->mpi_rank && ri[r].n_recv_entries) {
	MPI_Irecv(p, ri[r].n_recv * (me - mb), ddc_->mpi_type,
		  r, 0, ddc_->obj.comm, &patt2->recv_req[patt2->recv_cnt++]);
	p += ri[r].n_recv * (me - mb);
      }
    }  
    assert(p == (real_t*) patt2->recv_buf + patt2->n_recv * (me - mb));

    // post sends
    patt2->send_cnt = 0;
    p = (real_t*) patt2->send_buf;
    for (int r = 0; r < sub->mpi_size; r++) {
      if (r != sub->mpi_rank && ri[r].n_send_entries) {
	real_t *p0 = p;
	for (int i = 0; i < ri[r].n_send_entries; i++) {
	  struct mrc_ddc_sendrecv_entry *se = &ri[r].send_entry[i];
	  thrust::host_vector<uint> map_send(se->len * (me - mb));
	  map_setup(map_send, 0, mb, me, se->patch, se->ilo, se->ihi, cmflds);
	  real_t *buf = p;
	  for (auto cur : map_send) {
	    *buf++ = h_flds[cur];
	  }
	  p += map_send.size();
	}
	MPI_Isend(p0, ri[r].n_send * (me - mb), ddc_->mpi_type,
		  r, 0, ddc_->obj.comm, &patt2->send_req[patt2->send_cnt++]);
      }
    }  
    assert(p == (real_t*) patt2->send_buf + patt2->n_send * (me - mb));
  }

  // ----------------------------------------------------------------------
  // ddc_run_end

  void ddc_run_end(struct mrc_ddc_pattern2 *patt2, int mb, int me,
		   cuda_mfields& cmflds, thrust::host_vector<real_t>& h_flds,
		   void (*from_buf)(const thrust::host_vector<uint>& map,
				    real_t *buf, thrust::host_vector<real_t>& h_flds))
  {
    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc_);
    struct mrc_ddc_rank_info *ri = patt2->ri;

    MPI_Waitall(patt2->recv_cnt, patt2->recv_req, MPI_STATUSES_IGNORE);

    real_t* p = (real_t*) patt2->recv_buf;
    for (int r = 0; r < sub->mpi_size; r++) {
      if (r != sub->mpi_rank) {
	for (int i = 0; i < ri[r].n_recv_entries; i++) {
	  struct mrc_ddc_sendrecv_entry *re = &ri[r].recv_entry[i];
	  thrust::host_vector<uint> map_recv(re->len * (me - mb));
	  map_setup(map_recv, 0, mb, me, re->patch, re->ilo, re->ihi, cmflds);
	  from_buf(map_recv, p, h_flds);
	  p += map_recv.size();
	}
      }
    }

    MPI_Waitall(patt2->send_cnt, patt2->send_req, MPI_STATUSES_IGNORE);
  }

  void add_local(struct mrc_ddc_pattern2 *patt2, int mb, int me,
		 thrust::host_vector<real_t>& h_flds, cuda_mfields& cmflds)
  {
    thrust::host_vector<uint> map_send, map_recv;
    setup_local_maps(map_send, map_recv, patt2, mb, me, cmflds);

    thrust::host_vector<real_t> buf(map_send.size());
#if 0
    thrust::gather(map_send.begin(), map_send.end(), h_flds.begin(), buf.begin());
#else
    for (int i = 0; i < map_send.size(); i++) {
      buf[i] = h_flds[map_send[i]];
    }
#endif
    for (int i = 0; i < map_send.size(); i++) {
      h_flds[map_recv[i]] += buf[i];
    }
  }

  void fill_local(struct mrc_ddc_pattern2 *patt2, int mb, int me,
		  thrust::host_vector<real_t>& h_flds, cuda_mfields& cmflds)
  {
    thrust::host_vector<uint> map_send, map_recv;
    setup_local_maps(map_send, map_recv, patt2, mb, me, cmflds);

    thrust::host_vector<real_t> buf(map_send.size());
#if 0
    thrust::gather(map_send.begin(), map_send.end(), h_flds.begin(), buf.begin());
    thrust::scatter(buf.begin(), buf.end(), map_recv.begin(), h_flds.begin());
#else
    for (int i = 0; i < map_send.size(); i++) {
      buf[i] = h_flds[map_send[i]];
    }
    for (int i = 0; i < map_send.size(); i++) {
      h_flds[map_recv[i]] = buf[i];
    }
#endif
  }
  
  void setup_local_maps(thrust::host_vector<uint>& map_send, thrust::host_vector<uint>& map_recv,
			struct mrc_ddc_pattern2* patt2, int mb, int me, cuda_mfields& cmflds)
  {
    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc_);
    struct mrc_ddc_rank_info *ri = patt2->ri;

    uint buf_size = 0;
    for (int i = 0; i < ri[sub->mpi_rank].n_send_entries; i++) {
      struct mrc_ddc_sendrecv_entry *se = &ri[sub->mpi_rank].send_entry[i];
      if (se->ilo[0] == se->ihi[0] ||
	  se->ilo[1] == se->ihi[1] ||
	  se->ilo[2] == se->ihi[2]) { // FIXME, we shouldn't even create these
	continue;
      }
      uint size = (me - mb) * (se->ihi[0] - se->ilo[0]) * (se->ihi[1] - se->ilo[1]) * (se->ihi[2] - se->ilo[2]);
      buf_size += size;
    }

    map_send.resize(buf_size);
    map_recv.resize(buf_size);

    uint off = 0;
    for (int i = 0; i < ri[sub->mpi_rank].n_send_entries; i++) {
      struct mrc_ddc_sendrecv_entry *se = &ri[sub->mpi_rank].send_entry[i];
      struct mrc_ddc_sendrecv_entry *re = &ri[sub->mpi_rank].recv_entry[i];
      if (se->ilo[0] == se->ihi[0] ||
	  se->ilo[1] == se->ihi[1] ||
	  se->ilo[2] == se->ihi[2]) { // FIXME, we shouldn't even create these
	continue;
      }
      uint size = (me - mb) * (se->ihi[0] - se->ilo[0]) * (se->ihi[1] - se->ilo[1]) * (se->ihi[2] - se->ilo[2]);
      map_setup(map_send, off, mb, me, se->patch, se->ilo, se->ihi, cmflds);
      map_setup(map_recv, off, mb, me, re->patch, re->ilo, re->ihi, cmflds);
      off += size;
    }
  }

  static void map_setup(thrust::host_vector<uint>& map, uint off, int mb, int me, int p, int ilo[3], int ihi[3],
			cuda_mfields& cmflds)
  {
    auto cur = &map[off];
    for (int m = mb; m < me; m++) {
      for (int iz = ilo[2]; iz < ihi[2]; iz++) {
	for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	  for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	    *cur++ = cmflds.index(m, ix,iy,iz, p);
	  }
	}
      }
    }
  }

  // ----------------------------------------------------------------------
  // add_from_buf

  static void add_from_buf(const thrust::host_vector<uint>& map,
			   real_t* buf, thrust::host_vector<real_t>& h_flds)
  {
    for (auto cur : map) {
      h_flds[cur] += *buf++;
    }
  }
  
  // ----------------------------------------------------------------------
  // copy_from_buf

  static void copy_from_buf(const thrust::host_vector<uint>& map,
			    real_t* buf, thrust::host_vector<real_t>& h_flds)
  {
    for (auto cur : map) {
      h_flds[cur] = *buf++;
    }
  }

private:
  mrc_ddc* ddc_;
};

