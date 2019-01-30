
#pragma once

#include "psc_fields_cuda.h"
#include "cuda_mfields.h"
#include "fields.hxx"
#include "cuda_bits.h"

#include "mrc_ddc_private.h"

#include <unordered_map>

#include <thrust/gather.h>
#include <thrust/scatter.h>
#include <thrust/sort.h>

#define mrc_ddc_multi(ddc) mrc_to_subobj(ddc, struct mrc_ddc_multi)

template<typename real_t>
__global__
static void k_scatter(const real_t* buf, const uint* map, real_t* flds, unsigned int size)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < size) {
    flds[map[i]] = buf[i];
  }
}

template<typename real_t>
__global__
static void k_scatter_add(const real_t* buf, const uint* map, real_t* flds, unsigned int size)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < size) {
    atomicAdd(&flds[map[i]], buf[i]);
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

  // ======================================================================
  // Scatter
  
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

    void operator()(const thrust::device_vector<uint>& map,
		    const thrust::device_vector<real_t>& buf, thrust::device_ptr<real_t> d_flds)
    {
      if (buf.empty()) return;
      
      const int THREADS_PER_BLOCK = 256;
      dim3 dimGrid((buf.size() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
      k_scatter_add<<<dimGrid, THREADS_PER_BLOCK>>>(buf.data().get(), map.data().get(), d_flds.get(), buf.size());
      cuda_sync_if_enabled();
    }
  };

  struct Scatter
  {
    void operator()(const thrust::host_vector<uint>& map,
		    const thrust::host_vector<real_t>& buf, thrust::host_vector<real_t>& h_flds)
    {
      thrust::scatter(buf.begin(), buf.end(), map.begin(), h_flds.begin());
    }

    void operator()(const thrust::device_vector<uint>& map,
		    const thrust::device_vector<real_t>& buf, thrust::device_ptr<real_t> d_flds)
    {
#if 1
      thrust::scatter(buf.begin(), buf.end(), map.begin(), d_flds);
#else
      if (buf.empty()) return;

      const int THREADS_PER_BLOCK = 256;
      dim3 dimGrid((buf.size() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
      k_scatter<<<dimGrid, THREADS_PER_BLOCK>>>(buf.data().get(), map.data().get(), d_flds.get(), buf.size());
      cuda_sync_if_enabled();
#endif
    }
  };

  // ----------------------------------------------------------------------
  // ctor
  
  CudaBnd(const Grid_t& grid, Int3 ibn)
  {
    static struct mrc_ddc_funcs ddc_funcs;

    ddc_ = grid.mrc_domain().create_ddc();
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

      d_send = send;
      d_recv = recv;
      d_local_send = local_send;
      d_local_recv = local_recv;
      d_local_buf.resize(local_buf.size());
      d_send_buf.resize(send_buf.size());
      d_recv_buf.resize(recv_buf.size());
    }
    
    thrust::host_vector<uint> send, recv;
    thrust::host_vector<uint> local_send, local_recv;
    thrust::host_vector<real_t> local_buf;
    thrust::host_vector<real_t> send_buf;
    thrust::host_vector<real_t> recv_buf;

    thrust::device_vector<uint> d_recv, d_send;
    thrust::device_vector<uint> d_local_recv, d_local_send;
    thrust::device_vector<real_t> d_local_buf;
    thrust::device_vector<real_t> d_send_buf;
    thrust::device_vector<real_t> d_recv_buf;

    mrc_ddc_pattern2* patt;
    int mb, me;
  };

  // ----------------------------------------------------------------------
  // run

  template<typename T>
  void run(cuda_mfields& cmflds, int mb, int me, mrc_ddc_pattern2* patt2, std::unordered_map<int, Maps>& maps,
	   T scatter)
  {
    static int pr_ddc_run, pr_ddc_sync1, pr_ddc_sync2;
    if (!pr_ddc_run) {
      pr_ddc_run = prof_register("ddc_run", 1., 0, 0);
      pr_ddc_sync1 = prof_register("ddc_sync1", 1., 0, 0);
      pr_ddc_sync2 = prof_register("ddc_sync2", 1., 0, 0);
    }

#if 0
    prof_start(pr_ddc_sync1);
    MPI_Barrier(MPI_COMM_WORLD);
    prof_stop(pr_ddc_sync1);
#endif

    int key = mb + 100*me;
    auto map = maps.find(key);
    if (map == maps.cend()) {
      auto pair = maps.emplace(std::make_pair(key, Maps{ddc_, patt2, mb, me, cmflds}));
      map = pair.first;
    }

    prof_start(pr_ddc_run);
    ddc_run(map->second, patt2, mb, me, cmflds, scatter);
    prof_stop(pr_ddc_run);

#if 0
    prof_start(pr_ddc_sync2);
    MPI_Barrier(MPI_COMM_WORLD);
    prof_stop(pr_ddc_sync2);
#endif
  }
  
  // ----------------------------------------------------------------------
  // add_ghosts

  void add_ghosts(cuda_mfields& cmflds, int mb, int me)
  {
    mrc_ddc_multi* sub = mrc_ddc_multi(ddc_);

    run(cmflds, mb, me, &sub->add_ghosts2, maps_add_, ScatterAdd{});
  }
  
  // ----------------------------------------------------------------------
  // fill_ghosts

  void fill_ghosts(cuda_mfields& cmflds, int mb, int me)
  {
    // FIXME
    // I don't think we need as many points, and only stencil star
    // rather then box
    mrc_ddc_multi* sub = mrc_ddc_multi(ddc_);
    
    run(cmflds, mb, me, &sub->fill_ghosts2, maps_fill_, Scatter{});
  }

  // ----------------------------------------------------------------------
  // ddc_run

  template<typename S>
  void ddc_run(Maps& maps, mrc_ddc_pattern2* patt2, int mb, int me, cuda_mfields& cmflds,
	       S scatter)
  {
    static int pr_ddc0, pr_ddc1, pr_ddc2, pr_ddc3, pr_ddc4, pr_ddc5;
    static int pr_ddc6, pr_ddc7, pr_ddc8, pr_ddc9, pr_ddc10;
    if (!pr_ddc1) {
      pr_ddc0 = prof_register("ddc0", 1., 0, 0);
      pr_ddc1 = prof_register("ddc1", 1., 0, 0);
      pr_ddc2 = prof_register("ddc2", 1., 0, 0);
      pr_ddc3 = prof_register("ddc3", 1., 0, 0);
      pr_ddc4 = prof_register("ddc4", 1., 0, 0);
      pr_ddc5 = prof_register("ddc5", 1., 0, 0);
      pr_ddc6 = prof_register("ddc6", 1., 0, 0);
      pr_ddc7 = prof_register("ddc7", 1., 0, 0);
      pr_ddc8 = prof_register("ddc8", 1., 0, 0);
      pr_ddc9 = prof_register("ddc9", 1., 0, 0);
      pr_ddc10 = prof_register("ddc10", 1., 0, 0);
    }

#if 0
    thrust::device_ptr<real_t> d_flds{cmflds.data()};
    thrust::host_vector<real_t> h_flds{d_flds, d_flds + cmflds.n_fields * cmflds.n_cells};

    postReceives(maps);
    thrust::gather(maps.send.begin(), maps.send.end(), h_flds.begin(), maps.send_buf.begin());
    postSends(maps);

    MPI_Waitall(maps.patt->recv_cnt, maps.patt->recv_req, MPI_STATUSES_IGNORE);
    scatter(maps.recv, maps.recv_buf, h_flds);
    MPI_Waitall(maps.patt->send_cnt, maps.patt->send_req, MPI_STATUSES_IGNORE);

    // local part
    thrust::gather(maps.local_send.begin(), maps.local_send.end(), h_flds.begin(),
		   maps.local_buf.begin());
    scatter(maps.local_recv, maps.local_buf, h_flds);
    thrust::copy(h_flds.begin(), h_flds.end(), d_flds);
#else
    thrust::device_ptr<real_t> d_flds{cmflds.data()};
    prof_start(pr_ddc0);
    MPI_Barrier(MPI_COMM_WORLD);
    prof_stop(pr_ddc0);

    prof_start(pr_ddc1);
    postReceives(maps);
    prof_stop(pr_ddc1);

    prof_start(pr_ddc2);
    thrust::gather(maps.d_send.begin(), maps.d_send.end(), d_flds, maps.d_send_buf.begin());
    prof_stop(pr_ddc2);

    prof_start(pr_ddc3);
    thrust::copy(maps.d_send_buf.begin(), maps.d_send_buf.end(), maps.send_buf.begin());
    prof_stop(pr_ddc3);

    prof_start(pr_ddc4);
    postSends(maps);
    prof_stop(pr_ddc4);

    // local part
    prof_start(pr_ddc5);
    thrust::gather(maps.d_local_send.begin(), maps.d_local_send.end(), d_flds,
		   maps.d_local_buf.begin());
    prof_stop(pr_ddc5);

    prof_start(pr_ddc6);
    scatter(maps.d_local_recv, maps.d_local_buf, d_flds);
    prof_stop(pr_ddc6);

    prof_start(pr_ddc7);
    MPI_Waitall(maps.patt->recv_cnt, maps.patt->recv_req, MPI_STATUSES_IGNORE);
    prof_stop(pr_ddc7);

    prof_start(pr_ddc8);
    thrust::copy(maps.recv_buf.begin(), maps.recv_buf.end(), maps.d_recv_buf.begin());
    prof_stop(pr_ddc8);

    prof_start(pr_ddc9);
    scatter(maps.d_recv, maps.d_recv_buf, d_flds);
    prof_stop(pr_ddc9);

    prof_start(pr_ddc10);
    MPI_Waitall(maps.patt->send_cnt, maps.patt->send_req, MPI_STATUSES_IGNORE);
    prof_stop(pr_ddc10);
#endif

  }
  
  // ----------------------------------------------------------------------
  // postReceives
  
  void postReceives(Maps& maps)
  {
    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc_);
    struct mrc_ddc_rank_info *ri = maps.patt->ri;
    MPI_Datatype mpi_dtype = MPI_FLOAT; // FIXME Mfields_traits<Mfields>::mpi_dtype();
    int mm = maps.me - maps.mb;

    maps.patt->recv_cnt = 0;
    auto p_recv = maps.recv_buf.begin();
    for (int r = 0; r < sub->mpi_size; r++) {
      if (r != sub->mpi_rank && ri[r].n_recv_entries) {
	MPI_Irecv(&*p_recv, ri[r].n_recv * mm, mpi_dtype,
		  r, 0, ddc_->obj.comm, &maps.patt->recv_req[maps.patt->recv_cnt++]);
	p_recv += ri[r].n_recv * mm;
      }
    }
    assert(p_recv == maps.recv_buf.end());
  }

  // ----------------------------------------------------------------------
  // postSends
  
  void postSends(Maps& maps)
  {
    struct mrc_ddc_multi *sub = mrc_ddc_multi(ddc_);
    struct mrc_ddc_rank_info *ri = maps.patt->ri;
    MPI_Datatype mpi_dtype = MPI_FLOAT; // FIXME Mfields_traits<Mfields>::mpi_dtype();
    int mm = maps.me - maps.mb;

    maps.patt->send_cnt = 0;
    auto p_send = maps.send_buf.begin();
    for (int r = 0; r < sub->mpi_size; r++) {
      if (r != sub->mpi_rank && ri[r].n_send_entries) {
	MPI_Isend(&*p_send, ri[r].n_send * mm, mpi_dtype,
		  r, 0, ddc_->obj.comm, &maps.patt->send_req[maps.patt->send_cnt++]);
	p_send += ri[r].n_send * mm;
      }
    }  
    assert(p_send == maps.send_buf.end());
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

private:
  mrc_ddc* ddc_;
  std::unordered_map<int, Maps> maps_add_;
  std::unordered_map<int, Maps> maps_fill_;
};

