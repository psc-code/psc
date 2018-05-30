
#pragma once

#include "psc_fields_cuda.h"
#include "cuda_mfields.h"
#include "fields.hxx"

#include "mrc_ddc_private.h"

#include <thrust/gather.h>
#include <thrust/scatter.h>

#define mrc_ddc_multi(ddc) mrc_to_subobj(ddc, struct mrc_ddc_multi)

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

  // ======================================================================
  // Maps

  struct Maps
  {
    Maps(mrc_ddc* ddc, mrc_ddc_pattern2* patt2, int mb, int me, cuda_mfields& cmflds)
      : patt{patt2}, mb{mb}, me{me}
    {
      setup_remote_maps(send, recv, ddc, patt2, mb, me, cmflds);
      setup_local_maps(local_send, local_recv, ddc, patt2, mb, me, cmflds);
      local_buf.resize(local_send.size());
      send_buf.resize(send.size());
      recv_buf.resize(recv.size());
    }
    
    thrust::host_vector<uint> send, recv;
    thrust::host_vector<uint> local_send, local_recv;
    thrust::host_vector<real_t> local_buf;
    thrust::host_vector<real_t> send_buf;
    thrust::host_vector<real_t> recv_buf;
    mrc_ddc_pattern2* patt;
    int mb, me;
  };
  
  // ----------------------------------------------------------------------
  // add_ghosts

  void add_ghosts(Mfields& mflds, int mb, int me)
  {
    cuda_mfields& cmflds = *mflds.cmflds;
    mrc_ddc_multi* sub = mrc_ddc_multi(ddc_);
    
    Maps maps(ddc_, &sub->add_ghosts2, mb, me, cmflds);

    ddc_run(maps, &sub->add_ghosts2, mb, me, cmflds, ScatterAdd{});
  }
  
  // ----------------------------------------------------------------------
  // fill_ghosts

  void fill_ghosts(Mfields& mflds, int mb, int me)
  {
    // FIXME
    // I don't think we need as many points, and only stencil star
    // rather then box
    cuda_mfields& cmflds = *mflds.cmflds;
    mrc_ddc_multi* sub = mrc_ddc_multi(ddc_);
    
    Maps maps(ddc_, &sub->fill_ghosts2, mb, me, cmflds);

    ddc_run(maps, &sub->fill_ghosts2, mb, me, cmflds, Scatter{});
  }

  // ----------------------------------------------------------------------
  // ddc_run

  template<typename S>
  void ddc_run(Maps& maps, mrc_ddc_pattern2* patt2, int mb, int me, cuda_mfields& cmflds,
	       S scatter)
  {
    thrust::device_ptr<real_t> d_flds{cmflds.data()};
    thrust::host_vector<real_t> h_flds{d_flds, d_flds + cmflds.n_fields * cmflds.n_cells};

    ddc_run_begin(maps, h_flds);

    // local part
    thrust::gather(maps.local_send.begin(), maps.local_send.end(), h_flds.begin(),
		   maps.local_buf.begin());
    scatter(maps.local_recv, maps.local_buf, h_flds);

    ddc_run_end(maps, h_flds, scatter);

    thrust::copy(h_flds.begin(), h_flds.end(), d_flds);
  }
  
  // ----------------------------------------------------------------------
  // ddc_run_begin

  void ddc_run_begin(Maps& maps, thrust::host_vector<real_t>& h_flds)
  {
    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc_);
    auto patt2 = maps.patt;
    struct mrc_ddc_rank_info *ri = patt2->ri;
    MPI_Datatype mpi_dtype = MPI_FLOAT; // FIXME Mfields_traits<Mfields>::mpi_dtype();
    int mm = maps.me - maps.mb;

    // communicate aggregated buffers
    // post receives
    patt2->recv_cnt = 0;
    auto p_recv = maps.recv_buf.begin();
    for (int r = 0; r < sub->mpi_size; r++) {
      if (r != sub->mpi_rank && ri[r].n_recv_entries) {
	MPI_Irecv(&*p_recv, ri[r].n_recv * mm, mpi_dtype,
		  r, 0, ddc_->obj.comm, &patt2->recv_req[patt2->recv_cnt++]);
	p_recv += ri[r].n_recv * mm;
      }
    }
    assert(p_recv == maps.recv_buf.end());

    // gather what's to be sent
    thrust::gather(maps.send.begin(), maps.send.end(), h_flds.begin(), maps.send_buf.begin());

    // post sends
    patt2->send_cnt = 0;
    auto p_send = maps.send_buf.begin();
    for (int r = 0; r < sub->mpi_size; r++) {
      if (r != sub->mpi_rank && ri[r].n_send_entries) {
	MPI_Isend(&*p_send, ri[r].n_send * mm, mpi_dtype,
		  r, 0, ddc_->obj.comm, &patt2->send_req[patt2->send_cnt++]);
	p_send += ri[r].n_send * mm;
      }
    }  
    assert(p_send == maps.send_buf.end());
  }

  // ----------------------------------------------------------------------
  // ddc_run_end

  template<typename S>
  void ddc_run_end(Maps& maps,
		   thrust::host_vector<real_t>& h_flds,
		   S scatter)
  {
    MPI_Waitall(maps.patt->recv_cnt, maps.patt->recv_req, MPI_STATUSES_IGNORE);

    scatter(maps.recv, maps.recv_buf, h_flds);

    MPI_Waitall(maps.patt->send_cnt, maps.patt->send_req, MPI_STATUSES_IGNORE);
  }

  // ----------------------------------------------------------------------
  // setup_remote_maps

  static void setup_remote_maps(thrust::host_vector<uint>& map_send, thrust::host_vector<uint>& map_recv,
				mrc_ddc* ddc, struct mrc_ddc_pattern2* patt2, int mb, int me, cuda_mfields& cmflds)
  {
    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);
    struct mrc_ddc_rank_info *ri = patt2->ri;

    map_send.resize(patt2->n_send * (me - mb));
    map_recv.resize(patt2->n_recv * (me - mb));

    uint off_send = 0, off_recv = 0;
    for (int r = 0; r < sub->mpi_size; r++) {
      if (r == sub->mpi_rank) {
	continue;
      }

      for (int i = 0; i < ri[r].n_send_entries; i++) {
	struct mrc_ddc_sendrecv_entry *se = &ri[r].send_entry[i];
	map_setup(map_send, off_send, mb, me, se->patch, se->ilo, se->ihi, cmflds);
	off_send += se->len * (me - mb);
      }
      for (int i = 0; i < ri[r].n_recv_entries; i++) {
	struct mrc_ddc_sendrecv_entry *re = &ri[r].recv_entry[i];
	map_setup(map_recv, off_recv, mb, me, re->patch, re->ilo, re->ihi, cmflds);
	off_recv += re->len * (me - mb);
      }
    }
  }

  // ----------------------------------------------------------------------
  // setup_local_maps
  
  static void setup_local_maps(thrust::host_vector<uint>& map_send, thrust::host_vector<uint>& map_recv,
			       mrc_ddc* ddc, struct mrc_ddc_pattern2* patt2, int mb, int me, cuda_mfields& cmflds)
  {
    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc);
    struct mrc_ddc_rank_info *ri = patt2->ri;

    uint buf_size = 0;
    for (int i = 0; i < ri[sub->mpi_rank].n_send_entries; i++) {
      struct mrc_ddc_sendrecv_entry *se = &ri[sub->mpi_rank].send_entry[i];
      if (se->ilo[0] == se->ihi[0] ||
	  se->ilo[1] == se->ihi[1] ||
	  se->ilo[2] == se->ihi[2]) { // FIXME, we shouldn't even create these
	continue;
      }
      buf_size += se->len * (me - mb);
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
      uint size = se->len * (me - mb);
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

  struct ScatterAdd
  {
    void operator()(const thrust::host_vector<uint>& map,
		    const thrust::host_vector<real_t>& buf, thrust::host_vector<real_t>& h_flds)
    {
      auto p = buf.begin();
      for (auto cur : map) {
	h_flds[cur] += *p++;
      }
    }
  };

  struct Scatter
  {
    void operator()(const thrust::host_vector<uint>& map,
		    const thrust::host_vector<real_t>& buf, thrust::host_vector<real_t>& h_flds)
    {
#if 1
      thrust::scatter(buf.begin(), buf.end(), map.begin(), h_flds.begin());
#else
      auto p = buf.begin();
      for (auto cur : map) {
	h_flds[cur] = *p++;
      }
#endif
    }
  };

private:
  mrc_ddc* ddc_;
};

