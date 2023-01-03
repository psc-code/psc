
#pragma once

#include "psc_fields_cuda.h"
#include "fields.hxx"
#include "cuda_bits.h"

#include "mrc_ddc_private.h"

#include <unordered_map>

#include <thrust/gather.h>
#include <thrust/scatter.h>
#include <thrust/sort.h>

extern std::size_t mem_bnd;

#define mrc_ddc_multi(ddc) mrc_to_subobj(ddc, struct mrc_ddc_multi)

template <typename real_t>
__global__ static void k_scatter(const real_t* buf, const uint* map,
                                 real_t* flds, unsigned int size)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < size) {
    flds[map[i]] = buf[i];
  }
}

template <typename real_t>
__global__ static void k_scatter_add(const real_t* buf, const uint* map,
                                     real_t* flds, unsigned int size)
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
  using storage_type = MfieldsCuda::Storage;
  using real_t = storage_type::value_type;

  // ======================================================================
  // Scatter

  struct ScatterAdd
  {
    void operator()(const thrust::host_vector<uint>& map,
                    const thrust::host_vector<real_t>& buf,
                    thrust::host_vector<real_t>& h_flds)
    {
      auto p = buf.begin();
      for (auto cur : map) {
        h_flds[cur] += *p++;
      }
    }

    void operator()(const psc::device_vector<uint>& map,
                    const psc::device_vector<real_t>& buf,
                    thrust::device_ptr<real_t> d_flds)
    {
      if (buf.empty())
        return;

      const int THREADS_PER_BLOCK = 256;
      dim3 dimGrid((buf.size() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
      k_scatter_add<<<dimGrid, THREADS_PER_BLOCK>>>(
        buf.data().get(), map.data().get(), d_flds.get(), buf.size());
      cuda_sync_if_enabled();
    }
  };

  struct Scatter
  {
    void operator()(const thrust::host_vector<uint>& map,
                    const thrust::host_vector<real_t>& buf,
                    thrust::host_vector<real_t>& h_flds)
    {
      thrust::scatter(buf.begin(), buf.end(), map.begin(), h_flds.begin());
    }

    void operator()(const psc::device_vector<uint>& map,
                    const psc::device_vector<real_t>& buf,
                    thrust::device_ptr<real_t> d_flds)
    {
#if 1
      thrust::scatter(buf.begin(), buf.end(), map.begin(), d_flds);
#else
      if (buf.empty())
        return;

      const int THREADS_PER_BLOCK = 256;
      dim3 dimGrid((buf.size() + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
      k_scatter<<<dimGrid, THREADS_PER_BLOCK>>>(
        buf.data().get(), map.data().get(), d_flds.get(), buf.size());
      cuda_sync_if_enabled();
#endif
    }
  };

  // ----------------------------------------------------------------------
  // ctor

  CudaBnd(const Grid_t& grid, Int3 ibn)
    : balance_generation_cnt_{psc_balance_generation_cnt}
  {}

  CudaBnd(const CudaBnd& bnd) = delete;
  CudaBnd& operator=(const CudaBnd& bnd) = delete;

  // ======================================================================
  // Maps

  struct Maps
  {
    Maps(mrc_ddc* ddc, mrc_ddc_pattern2* patt2, int mb, int me,
         storage_type& mflds_gt, const Int3& mflds_ib)
      : patt{patt2}, mb{mb}, me{me}
    {
      setup_remote_maps(send, recv, ddc, patt2, mb, me, mflds_gt, mflds_ib);
      setup_local_maps(local_send, local_recv, ddc, patt2, mb, me, mflds_gt,
                       mflds_ib);
      send_buf.resize(send.size());
      recv_buf.resize(recv.size());

      d_send = send;
      mem_bnd += allocated_bytes(d_send);
      d_recv = recv;
      mem_bnd += allocated_bytes(d_recv);
      d_local_send = local_send;
      mem_bnd += allocated_bytes(d_local_send);
      d_local_recv = local_recv;
      mem_bnd += allocated_bytes(d_local_recv);
    }

    Maps(const Maps&) = delete;
    Maps(Maps&&) = default;

    ~Maps()
    {
      mem_bnd -= allocated_bytes(d_send);
      mem_bnd -= allocated_bytes(d_recv);
      mem_bnd -= allocated_bytes(d_local_send);
      mem_bnd -= allocated_bytes(d_local_recv);
    }

    thrust::host_vector<uint> send, recv;
    thrust::host_vector<uint> local_send, local_recv;
    thrust::host_vector<real_t> send_buf;
    thrust::host_vector<real_t> recv_buf;

    psc::device_vector<uint> d_recv, d_send;
    psc::device_vector<uint> d_local_recv, d_local_send;

    mrc_ddc_pattern2* patt;
    int mb, me;
  };

  // ----------------------------------------------------------------------
  // run

  template <typename T>
  void run(const Grid_t& grid, storage_type& mflds_gt, const Int3& mflds_ib,
           int mb, int me, mrc_ddc_pattern2* patt2,
           std::unordered_map<int, Maps>& maps, T scatter)
  {
    // static int pr_ddc_run, pr_ddc_sync1, pr_ddc_sync2;
    // if (!pr_ddc_run) {
    //   pr_ddc_run = prof_register("ddc_run", 1., 0, 0);
    //   pr_ddc_sync1 = prof_register("ddc_sync1", 1., 0, 0);
    //   pr_ddc_sync2 = prof_register("ddc_sync2", 1., 0, 0);
    // }

#if 0
    prof_start(pr_ddc_sync1);
    MPI_Barrier(MPI_COMM_WORLD);
    prof_stop(pr_ddc_sync1);
#endif

    int key = mb + 100 * me;
    auto map = maps.find(key);
    if (map == maps.cend()) {
      auto pair = maps.emplace(std::make_pair(
        key, Maps{grid.ddc(), patt2, mb, me, mflds_gt, mflds_ib}));
      map = pair.first;
    }

    // prof_start(pr_ddc_run);
    ddc_run(grid, map->second, patt2, mb, me, mflds_gt, mflds_ib, scatter);
    // prof_stop(pr_ddc_run);

#if 0
    prof_start(pr_ddc_sync2);
    MPI_Barrier(MPI_COMM_WORLD);
    prof_stop(pr_ddc_sync2);
#endif
  }

  // ----------------------------------------------------------------------
  // add_ghosts

  void add_ghosts(const Grid_t& grid, storage_type& mflds_gt,
                  const Int3& mflds_ib, int mb, int me)
  {
    if (psc_balance_generation_cnt != balance_generation_cnt_) {
      clear();
      balance_generation_cnt_ = psc_balance_generation_cnt;
    }

    mrc_ddc_multi* sub = mrc_ddc_multi(grid.ddc());

    run(grid, mflds_gt, mflds_ib, mb, me, &sub->add_ghosts2, maps_add_,
        ScatterAdd{});
  }

  // ----------------------------------------------------------------------
  // fill_ghosts

  void fill_ghosts(const Grid_t& grid, storage_type& mflds_gt,
                   const Int3& mflds_ib, int mb, int me)
  {
    if (psc_balance_generation_cnt != balance_generation_cnt_) {
      clear();
      balance_generation_cnt_ = psc_balance_generation_cnt;
    }

    // FIXME
    // I don't think we need as many points, and only stencil star
    // rather then box
    mrc_ddc_multi* sub = mrc_ddc_multi(grid.ddc());

    run(grid, mflds_gt, mflds_ib, mb, me, &sub->fill_ghosts2, maps_fill_,
        Scatter{});
  }

  // ----------------------------------------------------------------------
  // ddc_run

  template <typename S>
  void ddc_run(const Grid_t& grid, Maps& maps, mrc_ddc_pattern2* patt2, int mb,
               int me, storage_type& mflds_gt, const Int3& mflds_ib, S scatter)
  {
    // static int pr_ddc0, pr_ddc1, pr_ddc2, pr_ddc3, pr_ddc4, pr_ddc5;
    // static int pr_ddc6, pr_ddc7, pr_ddc8, pr_ddc9, pr_ddc10;
    // if (!pr_ddc1) {
    //   pr_ddc0 = prof_register("ddc0", 1., 0, 0);
    //   pr_ddc1 = prof_register("ddc1", 1., 0, 0);
    //   pr_ddc2 = prof_register("ddc2", 1., 0, 0);
    //   pr_ddc3 = prof_register("ddc3", 1., 0, 0);
    //   pr_ddc4 = prof_register("ddc4", 1., 0, 0);
    //   pr_ddc5 = prof_register("ddc5", 1., 0, 0);
    //   pr_ddc6 = prof_register("ddc6", 1., 0, 0);
    //   pr_ddc7 = prof_register("ddc7", 1., 0, 0);
    //   pr_ddc8 = prof_register("ddc8", 1., 0, 0);
    //   pr_ddc9 = prof_register("ddc9", 1., 0, 0);
    //   pr_ddc10 = prof_register("ddc10", 1., 0, 0);
    // }

#if 0
    thrust::device_ptr<real_t> d_flds{mflds.gt().data()};
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
    auto d_flds = mflds_gt.data();
    prof_barrier("ddc_run");

    // prof_start(pr_ddc1);
    postReceives(grid, maps);
    // prof_stop(pr_ddc1);

    {
      psc::device_vector<real_t> d_send_buf(maps.send_buf.size());
      // prof_start(pr_ddc2);
      thrust::gather(maps.d_send.begin(), maps.d_send.end(), d_flds,
                     d_send_buf.begin());
      // prof_stop(pr_ddc2);

      // prof_start(pr_ddc3);
      thrust::copy(d_send_buf.begin(), d_send_buf.end(), maps.send_buf.begin());
      // prof_stop(pr_ddc3);
    }

    // prof_start(pr_ddc4);
    postSends(grid, maps);
    // prof_stop(pr_ddc4);

    // local part
    {
      psc::device_vector<real_t> d_local_buf(maps.d_local_send.size());
      // prof_start(pr_ddc5);
      thrust::gather(maps.d_local_send.begin(), maps.d_local_send.end(), d_flds,
                     d_local_buf.begin());
      // prof_stop(pr_ddc5);

      // prof_start(pr_ddc6);
      scatter(maps.d_local_recv, d_local_buf, d_flds);
      // prof_stop(pr_ddc6);
    }

    {
      psc::device_vector<real_t> d_recv_buf(maps.recv_buf.size());
      // prof_start(pr_ddc7);
      MPI_Waitall(maps.patt->recv_cnt, maps.patt->recv_req,
                  MPI_STATUSES_IGNORE);
      // prof_stop(pr_ddc7);

      // prof_start(pr_ddc8);
      thrust::copy(maps.recv_buf.begin(), maps.recv_buf.end(),
                   d_recv_buf.begin());
      // prof_stop(pr_ddc8);

      // prof_start(pr_ddc9);
      scatter(maps.d_recv, d_recv_buf, d_flds);
      // prof_stop(pr_ddc9);
    }

    // prof_start(pr_ddc10);
    MPI_Waitall(maps.patt->send_cnt, maps.patt->send_req, MPI_STATUSES_IGNORE);
    // prof_stop(pr_ddc10);
#endif
  }

  // ----------------------------------------------------------------------
  // postReceives

  void postReceives(const Grid_t& grid, Maps& maps)
  {
    struct mrc_ddc_multi* sub = mrc_ddc_multi(grid.ddc());
    struct mrc_ddc_rank_info* ri = maps.patt->ri;
    MPI_Datatype mpi_dtype =
      MPI_FLOAT; // FIXME Mfields_traits<Mfields>::mpi_dtype();
    int mm = maps.me - maps.mb;

    maps.patt->recv_cnt = 0;
    auto p_recv = maps.recv_buf.begin();
    for (int r = 0; r < sub->mpi_size; r++) {
      if (r != sub->mpi_rank && ri[r].n_recv_entries) {
        MPI_Irecv(&*p_recv, ri[r].n_recv * mm, mpi_dtype, r, 0,
                  grid.ddc()->obj.comm,
                  &maps.patt->recv_req[maps.patt->recv_cnt++]);
        p_recv += ri[r].n_recv * mm;
      }
    }
    assert(p_recv == maps.recv_buf.end());
  }

  // ----------------------------------------------------------------------
  // postSends

  void postSends(const Grid_t& grid, Maps& maps)
  {
    struct mrc_ddc_multi* sub = mrc_ddc_multi(grid.ddc());
    struct mrc_ddc_rank_info* ri = maps.patt->ri;
    MPI_Datatype mpi_dtype =
      MPI_FLOAT; // FIXME Mfields_traits<Mfields>::mpi_dtype();
    int mm = maps.me - maps.mb;

    maps.patt->send_cnt = 0;
    auto p_send = maps.send_buf.begin();
    for (int r = 0; r < sub->mpi_size; r++) {
      if (r != sub->mpi_rank && ri[r].n_send_entries) {
        MPI_Isend(&*p_send, ri[r].n_send * mm, mpi_dtype, r, 0,
                  grid.ddc()->obj.comm,
                  &maps.patt->send_req[maps.patt->send_cnt++]);
        p_send += ri[r].n_send * mm;
      }
    }
    assert(p_send == maps.send_buf.end());
  }

  // ----------------------------------------------------------------------
  // setup_remote_maps

  static void setup_remote_maps(thrust::host_vector<uint>& map_send,
                                thrust::host_vector<uint>& map_recv,
                                mrc_ddc* ddc, struct mrc_ddc_pattern2* patt2,
                                int mb, int me, storage_type& mflds_gt,
                                const Int3& mflds_ib)
  {
    struct mrc_ddc_multi* sub = mrc_ddc_multi(ddc);
    struct mrc_ddc_rank_info* ri = patt2->ri;

    map_send.resize(patt2->n_send * (me - mb));
    map_recv.resize(patt2->n_recv * (me - mb));

    uint off_send = 0, off_recv = 0;
    for (int r = 0; r < sub->mpi_size; r++) {
      if (r == sub->mpi_rank) {
        continue;
      }

      for (int i = 0; i < ri[r].n_send_entries; i++) {
        struct mrc_ddc_sendrecv_entry* se = &ri[r].send_entry[i];
        map_setup(map_send, off_send, mb, me, se->patch, se->ilo, se->ihi,
                  mflds_gt, mflds_ib);
        off_send += se->len * (me - mb);
      }
      for (int i = 0; i < ri[r].n_recv_entries; i++) {
        struct mrc_ddc_sendrecv_entry* re = &ri[r].recv_entry[i];
        map_setup(map_recv, off_recv, mb, me, re->patch, re->ilo, re->ihi,
                  mflds_gt, mflds_ib);
        off_recv += re->len * (me - mb);
      }
    }
  }

  // ----------------------------------------------------------------------
  // setup_local_maps

  static void setup_local_maps(thrust::host_vector<uint>& map_send,
                               thrust::host_vector<uint>& map_recv,
                               mrc_ddc* ddc, struct mrc_ddc_pattern2* patt2,
                               int mb, int me, storage_type& mflds_gt,
                               const Int3& mflds_ib)
  {
    struct mrc_ddc_multi* sub = mrc_ddc_multi(ddc);
    struct mrc_ddc_rank_info* ri = patt2->ri;

    uint buf_size = 0;
    for (int i = 0; i < ri[sub->mpi_rank].n_send_entries; i++) {
      struct mrc_ddc_sendrecv_entry* se = &ri[sub->mpi_rank].send_entry[i];
      if (se->ilo[0] == se->ihi[0] || se->ilo[1] == se->ihi[1] ||
          se->ilo[2] == se->ihi[2]) { // FIXME, we shouldn't even create these
        continue;
      }
      buf_size += se->len * (me - mb);
    }

    map_send.resize(buf_size);
    map_recv.resize(buf_size);

    uint off = 0;
    for (int i = 0; i < ri[sub->mpi_rank].n_send_entries; i++) {
      struct mrc_ddc_sendrecv_entry* se = &ri[sub->mpi_rank].send_entry[i];
      struct mrc_ddc_sendrecv_entry* re = &ri[sub->mpi_rank].recv_entry[i];
      if (se->ilo[0] == se->ihi[0] || se->ilo[1] == se->ihi[1] ||
          se->ilo[2] == se->ihi[2]) { // FIXME, we shouldn't even create these
        continue;
      }
      uint size = se->len * (me - mb);
      map_setup(map_send, off, mb, me, se->patch, se->ilo, se->ihi, mflds_gt,
                mflds_ib);
      map_setup(map_recv, off, mb, me, re->patch, re->ilo, re->ihi, mflds_gt,
                mflds_ib);
      off += size;
    }
  }

  static void map_setup(thrust::host_vector<uint>& map, uint off, int mb,
                        int me, int p, int ilo[3], int ihi[3],
                        storage_type& mflds_gt, const Int3& ib)
  {
    auto cur = &map[off];
    for (int m = mb; m < me; m++) {
      for (int iz = ilo[2]; iz < ihi[2]; iz++) {
        for (int iy = ilo[1]; iy < ihi[1]; iy++) {
          for (int ix = ilo[0]; ix < ihi[0]; ix++) {
            *cur++ = &mflds_gt(ix - ib[0], iy - ib[1], iz - ib[2], m, p) -
                     mflds_gt.data();
          }
        }
      }
    }
  }

  void clear()
  {
    maps_add_.clear();
    maps_fill_.clear();
  }

private:
  std::unordered_map<int, Maps> maps_add_;
  std::unordered_map<int, Maps> maps_fill_;
  int balance_generation_cnt_;
};
