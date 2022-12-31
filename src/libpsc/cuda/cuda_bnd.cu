
#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "cuda_bits.h"

#include "psc.h"
#include "psc_fields_cuda.h"
#include "fields.hxx"

#include "cuda_bnd.cuh"

template <typename real_t>
__global__ static void k_scatter_add(const real_t* buf, const uint* map,
                                     real_t* flds, unsigned int size)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i < size) {
    atomicAdd(&flds[map[i]], buf[i]);
  }
}

template <typename real_t>
void ScatterAdd::operator()(const thrust::host_vector<uint>& map,
                            const thrust::host_vector<real_t>& buf,
                            thrust::host_vector<real_t>& h_flds)
{
  auto p = buf.begin();
  for (auto cur : map) {
    h_flds[cur] += *p++;
  }
}

template <typename real_t>
void ScatterAdd::operator()(const psc::device_vector<uint>& map,
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

// ======================================================================
// Scatter

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
void Scatter::operator()(const thrust::host_vector<uint>& map,
                         const thrust::host_vector<real_t>& buf,
                         thrust::host_vector<real_t>& h_flds)
{
  thrust::scatter(buf.begin(), buf.end(), map.begin(), h_flds.begin());
}

template <typename real_t>
void Scatter::operator()(const psc::device_vector<uint>& map,
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
  k_scatter<<<dimGrid, THREADS_PER_BLOCK>>>(buf.data().get(), map.data().get(),
                                            d_flds.get(), buf.size());
  cuda_sync_if_enabled();
#endif
}

// ======================================================================

template void ScatterAdd::operator()(const thrust::host_vector<uint>& map,
                                     const thrust::host_vector<float>& buf,
                                     thrust::host_vector<float>& h_flds);
template void ScatterAdd::operator()(const psc::device_vector<uint>& map,
                                     const psc::device_vector<float>& buf,
                                     thrust::device_ptr<float> d_flds);

template void Scatter::operator()(const thrust::host_vector<uint>& map,
                                  const thrust::host_vector<float>& buf,
                                  thrust::host_vector<float>& h_flds);
template void Scatter::operator()(const psc::device_vector<uint>& map,
                                  const psc::device_vector<float>& buf,
                                  thrust::device_ptr<float> d_flds);

// ======================================================================

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4

#define SW (2) // FIXME

template <bool lo, bool hi, typename E>
__global__ static void conducting_wall_H_y(E gt, Int3 ib)
{
  int iz = blockIdx.x * blockDim.x + threadIdx.x - SW;

  auto _d_flds = make_Fields3d<dim_xyz>(gt, ib);

  if (iz >= gt.shape(2) - SW)
    return;

  int my = gt.shape(1) - 2 * SW;

  if (lo) {
    _d_flds(HY, 0, -1, iz) = _d_flds(HY, 0, 1, iz);
    _d_flds(HX, 0, -1, iz) = -_d_flds(HX, 0, 0, iz);
    _d_flds(HZ, 0, -1, iz) = -_d_flds(HZ, 0, 0, iz);
  }

  if (hi) {
    _d_flds(HY, 0, my + 1, iz) = _d_flds(HY, 0, my - 1, iz);
    _d_flds(HX, 0, my, iz) = -_d_flds(HX, 0, my - 1, iz);
    _d_flds(HZ, 0, my, iz) = -_d_flds(HZ, 0, my - 1, iz);
  }
}

template <bool lo, bool hi, typename E>
__global__ static void conducting_wall_E_y(E gt, Int3 ib)
{
  int iz = blockIdx.x * blockDim.x + threadIdx.x - SW;

  auto _d_flds = make_Fields3d<dim_xyz>(gt, ib);

  if (iz >= gt.shape(2) - SW)
    return;

  int my = gt.shape(1) - 2 * SW;

  if (lo) {
    _d_flds(EX, 0, 0, iz) = 0.;
    _d_flds(EX, 0, -1, iz) = _d_flds(EX, 0, 1, iz);
    _d_flds(EY, 0, -1, iz) = -_d_flds(EY, 0, 0, iz);
    _d_flds(EZ, 0, 0, iz) = 0.;
    _d_flds(EZ, 0, -1, iz) = _d_flds(EZ, 0, 1, iz);
  }

  if (hi) {
    _d_flds(EX, 0, my, iz) = 0.;
    _d_flds(EX, 0, my + 1, iz) = _d_flds(EX, 0, my - 1, iz);
    _d_flds(EY, 0, my, iz) = -_d_flds(EY, 0, my - 1, iz);
    _d_flds(EZ, 0, my, iz) = 0.;
    _d_flds(EZ, 0, my + 1, iz) = _d_flds(EZ, 0, my - 1, iz);
  }
}

template <bool lo, bool hi, typename E>
__global__ static void conducting_wall_J_y(E gt, Int3 ib)
{
  int iz = blockIdx.x * blockDim.x + threadIdx.x - SW;

  auto _d_flds = make_Fields3d<dim_xyz>(gt, ib);

  if (iz >= gt.shape(2) - SW)
    return;

  int my = gt.shape(1) - 2 * SW;

  if (lo) {
    _d_flds(JYI, 0, 0, iz) -= _d_flds(JYI, 0, -1, iz);
    _d_flds(JYI, 0, -1, iz) = 0.;
    _d_flds(JXI, 0, 1, iz) += _d_flds(JXI, 0, -1, iz);
    _d_flds(JXI, 0, -1, iz) = 0.;
    _d_flds(JZI, 0, 1, iz) += _d_flds(JZI, 0, -1, iz);
    _d_flds(JZI, 0, -1, iz) = 0.;
  }

  if (hi) {
    _d_flds(JYI, 0, my - 1, iz) -= _d_flds(JYI, 0, my, iz);
    _d_flds(JYI, 0, my, iz) = 0.;
    _d_flds(JXI, 0, my - 1, iz) += _d_flds(JXI, 0, my + 1, iz);
    _d_flds(JXI, 0, my + 1, iz) = 0.;
    _d_flds(JZI, 0, my - 1, iz) += _d_flds(JZI, 0, my + 1, iz);
    _d_flds(JZI, 0, my + 1, iz) = 0.;
  }
}

template <bool lo, bool hi>
static void cuda_conducting_wall_H_y(MfieldsCuda& mflds, int p)
{
  int dimGrid = (mflds.gt().shape(2) + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  conducting_wall_H_y<lo, hi>
    <<<dimGrid, BLOCKSIZE_Z>>>(view_patch(mflds.gt(), p), -mflds.ibn());
  cuda_sync_if_enabled();
}

template <bool lo, bool hi>
static void cuda_conducting_wall_E_y(MfieldsCuda& mflds, int p)
{
  int dimGrid = (mflds.gt().shape(2) + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  conducting_wall_E_y<lo, hi>
    <<<dimGrid, BLOCKSIZE_Z>>>(view_patch(mflds.gt(), p), -mflds.ibn());
  cuda_sync_if_enabled();
}

template <bool lo, bool hi>
static void cuda_conducting_wall_J_y(MfieldsCuda& mflds, int p)
{
  int dimGrid = (mflds.gt().shape(2) + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z;
  conducting_wall_J_y<lo, hi>
    <<<dimGrid, BLOCKSIZE_Z>>>(view_patch(mflds.gt(), p), -mflds.ibn());
  cuda_sync_if_enabled();
}

void cuda_conducting_wall_H_lo_y(MfieldsCuda& mflds, int p)
{
  cuda_conducting_wall_H_y<true, false>(mflds, p);
}

void cuda_conducting_wall_H_hi_y(MfieldsCuda& mflds, int p)
{
  cuda_conducting_wall_H_y<false, true>(mflds, p);
}

void cuda_conducting_wall_H_lo_hi_y(MfieldsCuda& mflds, int p)
{
  cuda_conducting_wall_H_y<true, true>(mflds, p);
}

void cuda_conducting_wall_E_lo_y(MfieldsCuda& mflds, int p)
{
  cuda_conducting_wall_E_y<true, false>(mflds, p);
}

void cuda_conducting_wall_E_hi_y(MfieldsCuda& mflds, int p)
{
  cuda_conducting_wall_E_y<false, true>(mflds, p);
}

void cuda_conducting_wall_E_lo_hi_y(MfieldsCuda& mflds, int p)
{
  cuda_conducting_wall_E_y<true, true>(mflds, p);
}

void cuda_conducting_wall_J_lo_y(MfieldsCuda& mflds, int p)
{
  cuda_conducting_wall_J_y<true, false>(mflds, p);
}

void cuda_conducting_wall_J_hi_y(MfieldsCuda& mflds, int p)
{
  cuda_conducting_wall_J_y<false, true>(mflds, p);
}

void cuda_conducting_wall_J_lo_hi_y(MfieldsCuda& mflds, int p)
{
  cuda_conducting_wall_J_y<true, true>(mflds, p);
}
