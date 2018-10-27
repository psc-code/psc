
#ifndef CUDA_MPARTICLES_H
#define CUDA_MPARTICLES_H

#include "cuda_iface.h"
#include "cuda_mparticles_indexer.h"
#include "cuda_mparticles_sort.cuh"

#include "particles.hxx"
#include "psc_bits.h"
#include "cuda_bits.h"

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>

// ======================================================================
// cuda_mparticles_base

template<typename BS>
struct cuda_mparticles_base : cuda_mparticles_indexer<BS>
{
  cuda_mparticles_base(const Grid_t& grid);
  // copy constructor would work fine, but we don't want to copy everything
  // by accident
  cuda_mparticles_base(const cuda_mparticles<BS>&) = delete;

  void reserve_all(uint sinze);
  void resize_all(const uint *n_prts_by_patch);
  void get_size_all(uint *n_prts_by_patch);

  // protected:
  void to_device(float4 *xi4, float4 *pxi4, uint n_prts, uint off);
  void from_device(float4 *xi4, float4 *pxi4, uint n_prts, uint off);
  
  // per particle
  thrust::device_vector<float4> d_xi4;
  thrust::device_vector<float4> d_pxi4;

  // per block
  cuda_mparticles_sort2 by_block_;

  uint n_prts = 0;                       // total # of particles across all patches
  const Grid_t& grid_;
};

// ----------------------------------------------------------------------
// struct d_particle

struct d_particle {
  float xi[3];
  float kind_as_float;
  float pxi[3];
  float qni_wni;
};

#define LOAD_PARTICLE_POS(pp, d_xi4, n) do {				\
    float4 _xi4 = d_xi4[n];						\
    (pp).xi[0]         = _xi4.x;					\
    (pp).xi[1]         = _xi4.y;					\
    (pp).xi[2]         = _xi4.z;					\
    (pp).kind_as_float = _xi4.w;					\
} while (0)

#define LOAD_PARTICLE_MOM(pp, d_pxi4, n) do {				\
    float4 _pxi4 = d_pxi4[n];						\
    (pp).pxi[0]        = _pxi4.x;					\
    (pp).pxi[1]        = _pxi4.y;					\
    (pp).pxi[2]        = _pxi4.z;					\
    (pp).qni_wni       = _pxi4.w;					\
} while (0)

#define STORE_PARTICLE_POS(pp, d_xi4, n) do {				\
    float4 xi4 = { (pp).xi[0], (pp).xi[1], (pp).xi[2], (pp).kind_as_float }; \
    d_xi4[n] = xi4;							\
} while (0)

#define STORE_PARTICLE_MOM(pp, d_pxi4, n) do {				\
    float4 pxi4 = { (pp).pxi[0], (pp).pxi[1], (pp).pxi[2], (pp).qni_wni }; \
    d_pxi4[n] = pxi4;							\
} while (0)

// ----------------------------------------------------------------------
// cuda_mparticles

template<typename BS>
struct DMparticlesCuda;

template<typename _BS>
struct cuda_mparticles : cuda_mparticles_base<_BS>
{
  using BS = _BS;
  using particle_t = particle_cuda_t;
  using real_t = particle_t::real_t;
  using Real3 = Vec3<real_t>;
  using DMparticles = DMparticlesCuda<BS>;

  cuda_mparticles(const Grid_t& grid);

  void reserve_all(const uint *n_prts_by_patch);
  uint get_n_prts();
  void setup_internals();
  void inject_buf(const cuda_mparticles_prt *buf, const uint *buf_n_by_patch);
  void inject_buf(const particle_inject *buf, const uint *buf_n_by_patch);

  std::vector<cuda_mparticles_prt> get_particles(int beg, int end);
  std::vector<cuda_mparticles_prt> get_particles(int p);

  uint start(int p);
  
  template<typename F>
  void set_particles(uint p, F getter);

  void dump(const std::string& filename) const;
  void dump_by_patch(uint *n_prts_by_patch);

public:
  void find_block_indices_ids(thrust::device_vector<uint>& d_idx, thrust::device_vector<uint>& d_id);
  void find_cell_indices_ids(thrust::device_vector<uint>& d_idx, thrust::device_vector<uint>& d_id);
  void reorder();
  void reorder(const thrust::device_vector<uint>& d_id);
  void reorder_and_offsets(const thrust::device_vector<uint>& d_idx, const thrust::device_vector<uint>& d_id,
			   thrust::device_vector<uint>& d_off);
  void reorder_and_offsets_slow();
  void swap_alt();

  bool check_in_patch_unordered_slow();
  bool check_bidx_id_unordered_slow();
  bool check_ordered();
  bool check_bidx_after_push();

  void resize(uint n_prts);

public:
  // per particle
  thrust::device_vector<float4> d_alt_xi4;  // storage for out-of-place reordering of particle data
  thrust::device_vector<float4> d_alt_pxi4;

  std::vector<Real3> xb_by_patch; // lower left corner for each patch

  bool need_reorder = { false };            // particles haven't yet been put into their sorted order
};

template<typename BS_>
struct DMparticlesCuda : DParticleIndexer<BS_>
{
  using BS = BS_;
  using typename DParticleIndexer<BS>::real_t;
  
  static const int MAX_N_KINDS = 4;

  DMparticlesCuda(cuda_mparticles<BS>& cmprts)
    : DParticleIndexer<BS>{cmprts},
      dt_(cmprts.grid_.dt),
      fnqs_(cmprts.grid_.norm.fnqs),
      fnqxs_(cmprts.grid_.domain.dx[0] * fnqs_ / dt_),
      fnqys_(cmprts.grid_.domain.dx[1] * fnqs_ / dt_),
      fnqzs_(cmprts.grid_.domain.dx[2] * fnqs_ / dt_),
      dqs_(.5f * cmprts.grid_.norm.eta * dt_),
      xi4_(cmprts.d_xi4.data().get()), pxi4_(cmprts.d_pxi4.data().get()),
      alt_xi4_(cmprts.d_alt_xi4.data().get()), alt_pxi4_(cmprts.d_alt_pxi4.data().get()),
      off_(cmprts.by_block_.d_off.data().get()),
      bidx_(cmprts.by_block_.d_idx.data().get()),
      id_(cmprts.by_block_.d_id.data().get()),
      n_blocks_(cmprts.n_blocks)
  {
    auto& grid = cmprts.grid_;
    
    int n_kinds = grid.kinds.size();
    assert(n_kinds <= MAX_N_KINDS);
    for (int k = 0; k < n_kinds; k++) {
      dq_[k] = dqs_ * grid.kinds[k].q / grid.kinds[k].m;
      q_inv_[k] = 1.f / grid.kinds[k].q;
      q_[k] = grid.kinds[k].q;
      m_[k] = grid.kinds[k].m;
    }
  }

  __device__ real_t dt() const { return dt_; }
  __device__ real_t fnqs() const { return fnqs_; }
  __device__ real_t fnqxs() const { return fnqxs_; }
  __device__ real_t fnqys() const { return fnqys_; }
  __device__ real_t fnqzs() const { return fnqzs_; }
  __device__ real_t q_inv(int k) const { return q_inv_[k]; }
  __device__ real_t dq(int k) const { return dq_[k]; }
  __device__ real_t q(int k) const { return q_[k]; }
  __device__ real_t m(int k) const { return m_[k]; }

private:
  real_t dt_;
  real_t fnqs_;
  real_t fnqxs_, fnqys_, fnqzs_;
  real_t dqs_;
  real_t dq_[MAX_N_KINDS];
  real_t q_inv_[MAX_N_KINDS];
  real_t q_[MAX_N_KINDS];
  real_t m_[MAX_N_KINDS];
public:
  float4* xi4_;
  float4* pxi4_;
  float4* alt_xi4_;
  float4* alt_pxi4_;
  uint *off_;
  uint *bidx_;
  uint *id_;
  uint n_blocks_;
};

// ======================================================================
// cuda_mparticles implementation

// ----------------------------------------------------------------------
// set_particles

template<typename BS>
template<typename F>
void cuda_mparticles<BS>::set_particles(uint p, F getter)
{
  // FIXME, doing the copy here all the time would be nice to avoid
  // making sue we actually have a valid d_off would't hurt, either
  thrust::host_vector<uint> h_off(this->by_block_.d_off);

  uint off = h_off[p * this->n_blocks_per_patch];
  uint n_prts = h_off[(p+1) * this->n_blocks_per_patch] - off;
  
  thrust::host_vector<float4> xi4(n_prts);
  thrust::host_vector<float4> pxi4(n_prts);

  for (int n = 0; n < n_prts; n++) {
    struct cuda_mparticles_prt prt = getter(n);
    this->checkInPatchMod(prt.x);

    int kind = cuda_int_as_float(prt.kind);
    xi4[n].x  = prt.x[0];
    xi4[n].y  = prt.x[1];
    xi4[n].z  = prt.x[2];
    xi4[n].w  = kind;
    pxi4[n].x = prt.p[0];
    pxi4[n].y = prt.p[1];
    pxi4[n].z = prt.p[2];
    pxi4[n].w = prt.w * this->grid_.kinds[kind].q;
  }

  thrust::copy(xi4.begin(), xi4.end(), &this->d_xi4[off]);
  thrust::copy(pxi4.begin(), pxi4.end(), &this->d_pxi4[off]);
}



#endif
