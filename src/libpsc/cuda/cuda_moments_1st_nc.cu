
#include "cuda_mparticles.h"
#include "cuda_mfields.h"
#include "cuda_moments.cuh"

#define THREADS_PER_BLOCK (512)

// FIXME/TODO: we could do this w/o prior reordering, but currently the
// generic moment calculation code first reorders anyway (which it shouldn't)

// ======================================================================
// GCurr

class GCurr
{
public:
  DFields d_flds;

  __device__ GCurr(DFields _d_flds) :
    d_flds(_d_flds)
  {
  }

  __device__ void add(int m, int jx, int jy, int jz, float val)
  {
    float *addr = &d_flds(m, jx,jy,jz);
    atomicAdd(addr, val);
  }
};

// ======================================================================

// ----------------------------------------------------------------------
// rho_1st_nc_cuda_run

template<typename DMparticles, typename dim, bool REORDER>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
rho_1st_nc_cuda_run(DMparticles dmprts, DMFields dmflds)
{
  BlockSimple<typename DMparticles::BS, dim> current_block;
  if (!current_block.init(dmprts)) {
    return;
  }

  GCurr scurr(dmflds[current_block.p]);
  __syncthreads();

  int block_begin = dmprts.off_[current_block.bid];
  int block_end = dmprts.off_[current_block.bid + 1];
  for (int n : in_block_loop(block_begin, block_end)) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    if (REORDER) {
      uint id = dmprts.id_[n];
      LOAD_PARTICLE_POS(prt, dmprts.xi4_, id);
      LOAD_PARTICLE_MOM(prt, dmprts.pxi4_, id);
    } else {
      LOAD_PARTICLE_POS(prt, dmprts.xi4_, n);
      LOAD_PARTICLE_MOM(prt, dmprts.pxi4_, n);
    }

    float fnq = prt.qni_wni * dmprts.fnqs();
    
    int lf[3];
    float of[3];
    dmprts.find_idx_off_1st(prt.xi, lf, of, float(0.));

    if (dim::InvarX::value) { // FIXME, ugly...
      scurr.add(0, 0, lf[1]  , lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnq);
      scurr.add(0, 0, lf[1]+1, lf[2]  , (      of[1]) * (1.f - of[2]) * fnq);
      scurr.add(0, 0, lf[1]  , lf[2]+1, (1.f - of[1]) * (      of[2]) * fnq);
      scurr.add(0, 0, lf[1]+1, lf[2]+1, (      of[1]) * (      of[2]) * fnq);
    } else {
      scurr.add(0, lf[0]  , lf[1]  , lf[2]  , (1.f - of[0]) * (1.f - of[1]) * (1.f - of[2]) * fnq);
      scurr.add(0, lf[0]+1, lf[1]  , lf[2]  , (      of[0]) * (1.f - of[1]) * (1.f - of[2]) * fnq);
      scurr.add(0, lf[0]  , lf[1]+1, lf[2]  , (1.f - of[0]) * (      of[1]) * (1.f - of[2]) * fnq);
      scurr.add(0, lf[0]+1, lf[1]+1, lf[2]  , (      of[0]) * (      of[1]) * (1.f - of[2]) * fnq);
      scurr.add(0, lf[0]  , lf[1]  , lf[2]+1, (1.f - of[0]) * (1.f - of[1]) * (      of[2]) * fnq);
      scurr.add(0, lf[0]+1, lf[1]  , lf[2]+1, (      of[0]) * (1.f - of[1]) * (      of[2]) * fnq);
      scurr.add(0, lf[0]  , lf[1]+1, lf[2]+1, (1.f - of[0]) * (      of[1]) * (      of[2]) * fnq);
      scurr.add(0, lf[0]+1, lf[1]+1, lf[2]+1, (      of[0]) * (      of[1]) * (      of[2]) * fnq);
    }
  }
}

// ----------------------------------------------------------------------
// n_1st_cuda_run

template<typename BS, typename dim, bool REORDER>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
n_1st_cuda_run(DMparticlesCuda<BS> dmprts, DMFields dmflds)
{
  BlockSimple<BS, dim> current_block;
  if (!current_block.init(dmprts)) {
    return;
  }

  GCurr scurr(dmflds[current_block.p]);
  __syncthreads();

  int block_begin = dmprts.off_[current_block.bid];
  int block_end = dmprts.off_[current_block.bid + 1];
  for (int n : in_block_loop(block_begin, block_end)) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    if (REORDER) {
      uint id = dmprts.id_[n];
      LOAD_PARTICLE_POS(prt, dmprts.xi4_, id);
      LOAD_PARTICLE_MOM(prt, dmprts.pxi4_, id);
    } else {
      LOAD_PARTICLE_POS(prt, dmprts.xi4_, n);
      LOAD_PARTICLE_MOM(prt, dmprts.pxi4_, n);
    }

    int kind = __float_as_int(prt.kind_as_float);
    float wni = prt.qni_wni * dmprts.q_inv(kind);
    float fnq = wni * dmprts.fnqs();
    
    int lf[3];
    float of[3];
    dmprts.find_idx_off_1st(prt.xi, lf, of, float(-.5));

    if (dim::InvarX::value) { // FIXME, ugly...
      scurr.add(kind, 0, lf[1]  , lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnq);
      scurr.add(kind, 0, lf[1]+1, lf[2]  , (      of[1]) * (1.f - of[2]) * fnq);
      scurr.add(kind, 0, lf[1]  , lf[2]+1, (1.f - of[1]) * (      of[2]) * fnq);
      scurr.add(kind, 0, lf[1]+1, lf[2]+1, (      of[1]) * (      of[2]) * fnq);
    } else {
      scurr.add(kind, lf[0]  , lf[1]  , lf[2]  , (1.f - of[0]) * (1.f - of[1]) * (1.f - of[2]) * fnq);
      scurr.add(kind, lf[0]+1, lf[1]  , lf[2]  , (      of[0]) * (1.f - of[1]) * (1.f - of[2]) * fnq);
      scurr.add(kind, lf[0]  , lf[1]+1, lf[2]  , (1.f - of[0]) * (      of[1]) * (1.f - of[2]) * fnq);
      scurr.add(kind, lf[0]+1, lf[1]+1, lf[2]  , (      of[0]) * (      of[1]) * (1.f - of[2]) * fnq);
      scurr.add(kind, lf[0]  , lf[1]  , lf[2]+1, (1.f - of[0]) * (1.f - of[1]) * (      of[2]) * fnq);
      scurr.add(kind, lf[0]+1, lf[1]  , lf[2]+1, (      of[0]) * (1.f - of[1]) * (      of[2]) * fnq);
      scurr.add(kind, lf[0]  , lf[1]+1, lf[2]+1, (1.f - of[0]) * (      of[1]) * (      of[2]) * fnq);
      scurr.add(kind, lf[0]+1, lf[1]+1, lf[2]+1, (      of[0]) * (      of[1]) * (      of[2]) * fnq);
    }
  }
}

// ----------------------------------------------------------------------
// CudaMoments1stNcRho::operator()

template<typename CudaMparticles, typename dim>
void CudaMoments1stNcRho<CudaMparticles, dim>::operator()(CudaMparticles& cmprts, struct cuda_mfields *cmres)
{
  cmprts.reorder(); // FIXME/OPT?
  
  if (!cmprts.need_reorder) {
    invoke<false>(cmprts, cmres);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// CudaMoments1stNcRho::invoke

template<typename CudaMparticles, typename dim>
template<bool REORDER>
void CudaMoments1stNcRho<CudaMparticles, dim>::invoke(CudaMparticles& cmprts, struct cuda_mfields *cmres)
{
  dim3 dimGrid = BlockSimple<typename CudaMparticles::BS, dim>::dimGrid(cmprts);

  rho_1st_nc_cuda_run<typename CudaMparticles::DMparticles, dim, REORDER>
    <<<dimGrid, THREADS_PER_BLOCK>>>(cmprts, *cmres);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// CudaMoments1stNcN::operator()

template<typename CudaMparticles, typename dim>
void CudaMoments1stNcN<CudaMparticles, dim>::operator()(CudaMparticles& cmprts, struct cuda_mfields *cmres)
{
  cmprts.reorder(); // FIXME/OPT?

  if (!cmprts.need_reorder) {
    invoke<false>(cmprts, cmres);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// CudaMoments1stNcN::invoke

template<typename CudaMparticles, typename dim>
template<bool REORDER>
void CudaMoments1stNcN<CudaMparticles, dim>::invoke(CudaMparticles& cmprts, struct cuda_mfields *cmres)
{
  dim3 dimGrid = BlockSimple<typename CudaMparticles::BS, dim>::dimGrid(cmprts);

  n_1st_cuda_run<typename CudaMparticles::BS, dim, REORDER>
    <<<dimGrid, THREADS_PER_BLOCK>>>(cmprts, *cmres);
  cuda_sync_if_enabled();
}

template struct CudaMoments1stNcRho<cuda_mparticles<BS144>, dim_yz>;
template struct CudaMoments1stNcN<cuda_mparticles<BS144>, dim_yz>;

template struct CudaMoments1stNcN<cuda_mparticles<BS444>, dim_xyz>;
