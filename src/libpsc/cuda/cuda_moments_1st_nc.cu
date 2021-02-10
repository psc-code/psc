
#include "cuda_mparticles.cuh"
#include "cuda_mfields.h"
#include "cuda_moments.cuh"
#include "bs.hxx"
#include "pushp.hxx"

#define THREADS_PER_BLOCK (512)

// FIXME/TODO: we could do this w/o prior reordering, but currently the
// generic moment calculation code first reorders anyway (which it shouldn't)

// ======================================================================

template <typename DIM>
class Deposit
{
public:
  using R = float;

  GT_INLINE Deposit(DFields& dflds, R fnq) : dflds_(dflds), fnq_(fnq) {}

  __device__ void operator()(int m, int lf[3], R of[3], R val, dim_yz tag)
  {
    R what = fnq_ * val;

    atomicAdd(&dflds_(m, 0, lf[1], lf[2]),
              (1.f - of[1]) * (1.f - of[2]) * what);
    atomicAdd(&dflds_(m, 0, lf[1] + 1, lf[2]), (of[1]) * (1.f - of[2]) * what);
    atomicAdd(&dflds_(m, 0, lf[1], lf[2] + 1), (1.f - of[1]) * (of[2]) * what);
    atomicAdd(&dflds_(m, 0, lf[1] + 1, lf[2] + 1), (of[1]) * (of[2]) * what);
  }

  __device__ void operator()(int m, int lf[3], R of[3], R val, dim_xyz tag)
  {
    R what = fnq_ * val;

    atomicAdd(&dflds_(m, lf[0], lf[1], lf[2]),
              (1.f - of[0]) * (1.f - of[1]) * (1.f - of[2]) * what);
    atomicAdd(&dflds_(m, lf[0] + 1, lf[1], lf[2]),
              (of[0]) * (1.f - of[1]) * (1.f - of[2]) * what);
    atomicAdd(&dflds_(m, lf[0], lf[1] + 1, lf[2]),
              (1.f - of[0]) * (of[1]) * (1.f - of[2]) * what);
    atomicAdd(&dflds_(m, lf[0] + 1, lf[1] + 1, lf[2]),
              (of[0]) * (of[1]) * (1.f - of[2]) * what);
    atomicAdd(&dflds_(m, lf[0], lf[1], lf[2] + 1),
              (1.f - of[0]) * (1.f - of[1]) * (of[2]) * what);
    atomicAdd(&dflds_(m, lf[0] + 1, lf[1], lf[2] + 1),
              (of[0]) * (1.f - of[1]) * (of[2]) * what);
    atomicAdd(&dflds_(m, lf[0], lf[1] + 1, lf[2] + 1),
              (1.f - of[0]) * (of[1]) * (of[2]) * what);
    atomicAdd(&dflds_(m, lf[0] + 1, lf[1] + 1, lf[2] + 1),
              (of[0]) * (of[1]) * (of[2]) * what);
  }

  __device__ void operator()(int m, int lf[3], R of[3], R val)
  {
    (*this)(m, lf, of, val, DIM{});
  }

private:
  DFields& dflds_;
  R fnq_;
};

// ----------------------------------------------------------------------
// rho_1st_nc_cuda_run

template <typename DMparticles, typename dim, bool REORDER>
__global__ static void __launch_bounds__(THREADS_PER_BLOCK, 3)
  rho_1st_nc_cuda_run(DMparticles dmprts, DMFields dmflds)
{
  BlockSimple<typename DMparticles::BS, dim> current_block;
  if (!current_block.init(dmprts)) {
    return;
  }

  DFields dflds(dmflds[current_block.p]);

  __syncthreads();

  int block_begin = dmprts.off_[current_block.bid];
  int block_end = dmprts.off_[current_block.bid + 1];
  for (int n : in_block_loop(block_begin, block_end)) {
    if (n < block_begin) {
      continue;
    }
    const auto prt =
      REORDER ? dmprts.storage[dmprts.id_[n]] : dmprts.storage[n];

    float fnq = dmprts.prt_w(prt) * dmprts.fnqs();
    float q = dmprts.prt_q(prt);

    int lf[3];
    float of[3];
    dmprts.find_idx_off_1st(prt.x, lf, of, float(0.));

    Deposit<dim> deposit(dflds, fnq);

    deposit(0, lf, of, q);
  }
}

// ----------------------------------------------------------------------
// n_1st_cuda_run

template <typename BS, typename dim, bool REORDER>
__global__ static void __launch_bounds__(THREADS_PER_BLOCK, 3)
  n_1st_cuda_run(DMparticlesCuda<BS> dmprts, DMFields dmflds)
{
  BlockSimple<BS, dim> current_block;
  if (!current_block.init(dmprts)) {
    return;
  }

  DFields dflds(dmflds[current_block.p]);

  __syncthreads();

  int block_begin = dmprts.off_[current_block.bid];
  int block_end = dmprts.off_[current_block.bid + 1];
  for (int n : in_block_loop(block_begin, block_end)) {
    if (n < block_begin) {
      continue;
    }
    const auto prt =
      REORDER ? dmprts.storage[dmprts.id_[n]] : dmprts.storage[n];

    int kind = prt.kind;
    float fnq = dmprts.prt_w(prt) * dmprts.fnqs();
    float q = dmprts.prt_q(prt);

    int lf[3];
    float of[3];
    dmprts.find_idx_off_1st(prt.x, lf, of, float(-.5));

    Deposit<dim> deposit(dflds, fnq);

    deposit(kind, lf, of, 1.f);
  }
}

// ----------------------------------------------------------------------
// all_1st_cuda_run

template <typename BS, typename dim, bool REORDER>
__global__ static void __launch_bounds__(THREADS_PER_BLOCK, 3)
  all_1st_cuda_run(DMparticlesCuda<BS> dmprts, DMFields dmflds)
{
  BlockSimple<BS, dim> current_block;
  if (!current_block.init(dmprts)) {
    return;
  }

  DFields dflds(dmflds[current_block.p]);

  __syncthreads();

  int block_begin = dmprts.off_[current_block.bid];
  int block_end = dmprts.off_[current_block.bid + 1];
  for (int n : in_block_loop(block_begin, block_end)) {
    if (n < block_begin) {
      continue;
    }
    const auto prt =
      REORDER ? dmprts.storage[dmprts.id_[n]] : dmprts.storage[n];

    float fnq = dmprts.prt_w(prt) * dmprts.fnqs();
    float q = dmprts.prt_q(prt);
    float m = dmprts.prt_m(prt);

    int lf[3];
    float of[3];
    dmprts.find_idx_off_1st(prt.x, lf, of, float(-.5));

    AdvanceParticle<float, dim> advance{dmprts.dt()};
    auto v = advance.calc_v(prt.u);

    Deposit<dim> deposit(dflds, fnq);

    int n_moments = 13;
    int mm = prt.kind * n_moments;
    deposit(mm + 0, lf, of, q);
    deposit(mm + 1, lf, of, q * v[0]);
    deposit(mm + 2, lf, of, q * v[1]);
    deposit(mm + 3, lf, of, q * v[2]);
    deposit(mm + 4, lf, of, m * prt.u[0]);
    deposit(mm + 5, lf, of, m * prt.u[1]);
    deposit(mm + 6, lf, of, m * prt.u[2]);
    deposit(mm + 7, lf, of, m * prt.u[0] * v[0]);
    deposit(mm + 8, lf, of, m * prt.u[1] * v[1]);
    deposit(mm + 9, lf, of, m * prt.u[2] * v[2]);
    deposit(mm + 10, lf, of, m * prt.u[0] * v[1]);
    deposit(mm + 11, lf, of, m * prt.u[1] * v[2]);
    deposit(mm + 12, lf, of, m * prt.u[2] * v[0]);
  }
}

// ----------------------------------------------------------------------
// CudaMoments1stNcRho::operator()

template <typename CudaMparticles, typename dim>
void CudaMoments1stNcRho<CudaMparticles, dim>::operator()(
  CudaMparticles& cmprts, struct cuda_mfields* cmres)
{
  if (cmprts.n_prts == 0) {
    return;
  }
  cmprts.reorder(); // FIXME/OPT?

  if (!cmprts.need_reorder) {
    dim3 dimGrid =
      BlockSimple<typename CudaMparticles::BS, dim>::dimGrid(cmprts);

    rho_1st_nc_cuda_run<typename CudaMparticles::DMparticles, dim, false>
      <<<dimGrid, THREADS_PER_BLOCK>>>(cmprts, *cmres);
    cuda_sync_if_enabled();
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// CudaMoments1stN::operator()

template <typename CudaMparticles, typename dim>
void CudaMoments1stN<CudaMparticles, dim>::operator()(
  CudaMparticles& cmprts, struct cuda_mfields* cmres)
{
  static int pr, pr_1;
  if (!pr) {
    pr = prof_register("cuda_mom_n", 1, 0, 0);
    pr_1 = prof_register("cuda_mom_n_reorder", 1, 0, 0);
  }

  prof_start(pr);
  if (cmprts.n_prts == 0) {
    return;
  }

  prof_start(pr_1);
  cmprts.reorder(); // FIXME/OPT?
  prof_stop(pr_1);

  if (!cmprts.need_reorder) {
    dim3 dimGrid =
      BlockSimple<typename CudaMparticles::BS, dim>::dimGrid(cmprts);

    n_1st_cuda_run<typename CudaMparticles::BS, dim, false>
      <<<dimGrid, THREADS_PER_BLOCK>>>(cmprts, *cmres);
    cuda_sync_if_enabled();
  } else {
    assert(0);
  }
  prof_stop(pr);
}

// ----------------------------------------------------------------------
// CudaMoments1stAll::operator()

template <typename CudaMparticles, typename dim>
void CudaMoments1stAll<CudaMparticles, dim>::operator()(
  CudaMparticles& cmprts, struct cuda_mfields* cmres)
{
  static int pr, pr_1;
  if (!pr) {
    pr = prof_register("cuda_mom_all", 1, 0, 0);
    pr_1 = prof_register("cuda_mom_all_reorder", 1, 0, 0);
  }

  // prof_start(pr);
  if (cmprts.n_prts == 0) {
    return;
  }

  // prof_start(pr_1);
  cmprts.reorder(); // FIXME/OPT?
  // prof_stop(pr_1);

  if (!cmprts.need_reorder) {
    dim3 dimGrid =
      BlockSimple<typename CudaMparticles::BS, dim>::dimGrid(cmprts);

    all_1st_cuda_run<typename CudaMparticles::BS, dim, false>
      <<<dimGrid, THREADS_PER_BLOCK>>>(cmprts, *cmres);
    cuda_sync_if_enabled();
  } else {
    assert(0);
  }
  // prof_stop(pr);
}

template struct CudaMoments1stNcRho<cuda_mparticles<BS144>, dim_yz>;
template struct CudaMoments1stN<cuda_mparticles<BS144>, dim_yz>;
template struct CudaMoments1stAll<cuda_mparticles<BS144>, dim_yz>;

template struct CudaMoments1stNcRho<cuda_mparticles<BS444>, dim_xyz>;
template struct CudaMoments1stN<cuda_mparticles<BS444>, dim_xyz>;
template struct CudaMoments1stAll<cuda_mparticles<BS444>, dim_xyz>;
