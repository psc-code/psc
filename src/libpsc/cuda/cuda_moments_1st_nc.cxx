
#include "cuda_mparticles.hxx"
#include "cuda_moments.hxx"
#include "bs.hxx"
#include "pushp.hxx"
#include "fields.hxx"

#define THREADS_PER_BLOCK (512)

// FIXME/TODO: we could do this w/o prior reordering, but currently the
// generic moment calculation code first reorders anyway (which it shouldn't)

// ======================================================================

template <typename DIM>
class Deposit
{
public:
  using R = float;

  template <typename E>
  GT_INLINE void operator()(E& flds, int m, int lf[3], R of[3], R val,
                            dim_yz tag)
  {
    atomicAdd(&flds(m, 0, lf[1], lf[2]), (1.f - of[1]) * (1.f - of[2]) * val);
    atomicAdd(&flds(m, 0, lf[1] + 1, lf[2]), (of[1]) * (1.f - of[2]) * val);
    atomicAdd(&flds(m, 0, lf[1], lf[2] + 1), (1.f - of[1]) * (of[2]) * val);
    atomicAdd(&flds(m, 0, lf[1] + 1, lf[2] + 1), (of[1]) * (of[2]) * val);
  }

  template <typename E>
  GT_INLINE void operator()(E& flds, int m, int lf[3], R of[3], R val,
                            dim_xyz tag)
  {
    atomicAdd(&flds(m, lf[0], lf[1], lf[2]),
              (1.f - of[0]) * (1.f - of[1]) * (1.f - of[2]) * val);
    atomicAdd(&flds(m, lf[0] + 1, lf[1], lf[2]),
              (of[0]) * (1.f - of[1]) * (1.f - of[2]) * val);
    atomicAdd(&flds(m, lf[0], lf[1] + 1, lf[2]),
              (1.f - of[0]) * (of[1]) * (1.f - of[2]) * val);
    atomicAdd(&flds(m, lf[0] + 1, lf[1] + 1, lf[2]),
              (of[0]) * (of[1]) * (1.f - of[2]) * val);
    atomicAdd(&flds(m, lf[0], lf[1], lf[2] + 1),
              (1.f - of[0]) * (1.f - of[1]) * (of[2]) * val);
    atomicAdd(&flds(m, lf[0] + 1, lf[1], lf[2] + 1),
              (of[0]) * (1.f - of[1]) * (of[2]) * val);
    atomicAdd(&flds(m, lf[0], lf[1] + 1, lf[2] + 1),
              (1.f - of[0]) * (of[1]) * (of[2]) * val);
    atomicAdd(&flds(m, lf[0] + 1, lf[1] + 1, lf[2] + 1),
              (of[0]) * (of[1]) * (of[2]) * val);
  }

  template <typename E>
  GT_INLINE void operator()(E& flds, int m, int lf[3], R of[3], R val)
  {
    (*this)(flds, m, lf, of, val, DIM{});
  }
};

// ----------------------------------------------------------------------
// rho_1st_nc_cuda_run

template <typename DMparticles, typename dim, bool REORDER, typename E>
__global__ static void __launch_bounds__(THREADS_PER_BLOCK, 3)
  rho_1st_nc_cuda_run(DMparticles dmprts, E mflds_gt, Int3 ib)
{
  BlockSimple<typename DMparticles::BS, dim> current_block;
  if (!current_block.init(dmprts)) {
    return;
  }

  auto gt = view_patch(mflds_gt, current_block.p);
  auto flds = make_Fields3d<dim>(gt, ib);
  Deposit<dim> deposit;
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

    deposit(flds, 0, lf, of, q * fnq);
  }
}

// ----------------------------------------------------------------------
// n_1st_cuda_run

template <typename BS, typename dim, bool REORDER, typename E>
__global__ static void __launch_bounds__(THREADS_PER_BLOCK, 3)
  n_1st_cuda_run(DMparticlesCuda<BS> dmprts, E mflds_gt, Int3 ib)
{
  BlockSimple<BS, dim> current_block;
  if (!current_block.init(dmprts)) {
    return;
  }

  auto gt = view_patch(mflds_gt, current_block.p);
  auto flds = make_Fields3d<dim>(gt, ib);
  Deposit<dim> deposit;
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

    int lf[3];
    float of[3];
    dmprts.find_idx_off_1st(prt.x, lf, of, float(-.5));

    deposit(flds, kind, lf, of, fnq);
  }
}

// ----------------------------------------------------------------------
// all_1st_cuda_run

template <typename BS, typename dim, bool REORDER, typename E>
__global__ static void __launch_bounds__(THREADS_PER_BLOCK, 3)
  all_1st_cuda_run(DMparticlesCuda<BS> dmprts, E mflds_gt, Int3 ib)
{
  BlockSimple<BS, dim> current_block;
  if (!current_block.init(dmprts)) {
    return;
  }

  auto gt = view_patch(mflds_gt, current_block.p);
  auto flds = make_Fields3d<dim>(gt, ib);
  Deposit<dim> deposit;
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

    int n_moments = 13;
    int mm = prt.kind * n_moments;
    deposit(flds, mm + 0, lf, of, fnq * q);
    deposit(flds, mm + 1, lf, of, fnq * q * v[0]);
    deposit(flds, mm + 2, lf, of, fnq * q * v[1]);
    deposit(flds, mm + 3, lf, of, fnq * q * v[2]);
    deposit(flds, mm + 4, lf, of, fnq * m * prt.u[0]);
    deposit(flds, mm + 5, lf, of, fnq * m * prt.u[1]);
    deposit(flds, mm + 6, lf, of, fnq * m * prt.u[2]);
    deposit(flds, mm + 7, lf, of, fnq * m * prt.u[0] * v[0]);
    deposit(flds, mm + 8, lf, of, fnq * m * prt.u[1] * v[1]);
    deposit(flds, mm + 9, lf, of, fnq * m * prt.u[2] * v[2]);
    deposit(flds, mm + 10, lf, of, fnq * m * prt.u[0] * v[1]);
    deposit(flds, mm + 11, lf, of, fnq * m * prt.u[1] * v[2]);
    deposit(flds, mm + 12, lf, of, fnq * m * prt.u[2] * v[0]);
  }
}

// ----------------------------------------------------------------------
// CudaMoments1stNcRho::operator()

template <typename CudaMparticles, typename dim>
void CudaMoments1stNcRho<CudaMparticles, dim>::operator()(
  CudaMparticles& cmprts, MfieldsCuda::Storage& mres_gt, const Int3& mres_ib)
{
  if (cmprts.n_prts == 0) {
    return;
  }
  cmprts.reorder(); // FIXME/OPT?

  if (!cmprts.need_reorder) {
    dim3 dimGrid =
      BlockSimple<typename CudaMparticles::BS, dim>::dimGrid(cmprts);

    auto k_mres_gt = mres_gt.to_kernel();
    rho_1st_nc_cuda_run<typename CudaMparticles::DMparticles, dim, false>
      <<<dimGrid, THREADS_PER_BLOCK>>>(cmprts, k_mres_gt, mres_ib);
    cuda_sync_if_enabled();
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// CudaMoments1stN::operator()

template <typename CudaMparticles, typename dim>
void CudaMoments1stN<CudaMparticles, dim>::operator()(
  CudaMparticles& cmprts, MfieldsCuda::Storage& mres_gt, const Int3& mres_ib)
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

    auto k_mres_gt = mres_gt.to_kernel();
    n_1st_cuda_run<typename CudaMparticles::BS, dim, false>
      <<<dimGrid, THREADS_PER_BLOCK>>>(cmprts, k_mres_gt, mres_ib);
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
  CudaMparticles& cmprts, MfieldsCuda::Storage& mres_gt, const Int3& mres_ib)
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

    auto k_mres_gt = mres_gt.to_kernel();
    all_1st_cuda_run<typename CudaMparticles::BS, dim, false>
      <<<dimGrid, THREADS_PER_BLOCK>>>(cmprts, k_mres_gt, mres_ib);
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
