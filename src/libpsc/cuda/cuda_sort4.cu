
#include "psc_cuda.h"
#include "particles_cuda.h"

#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/device_vector.h>

static void
do_sort(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  thrust::device_ptr<unsigned int> d_bidx(mprts_cuda->d_bidx);
  thrust::device_ptr<unsigned int> d_ids(mprts_cuda->d_ids);
  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_ids(d_ids, d_ids + mprts_cuda->nr_prts);

  thrust::host_vector<unsigned int> h_cnts(mprts_cuda->nr_total_blocks + 1);
  for (int n = 0; n < mprts_cuda->nr_prts; n++) {
    assert(h_bidx[n] <= mprts_cuda->nr_total_blocks);
    h_cnts[h_bidx[n]]++;
  }
  thrust::host_vector<unsigned int> h_sums(mprts_cuda->nr_total_blocks + 1);
  thrust::exclusive_scan(h_cnts.begin(), h_cnts.end(), h_sums.begin());
  thrust::host_vector<unsigned int> h_bidx2(mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_ids2(mprts_cuda->nr_prts);
  for (int n = 0; n < mprts_cuda->nr_prts; n++) {
    int nn = h_sums[h_bidx[n]]++;
    h_ids2[nn] = n;
    h_bidx2[nn] = h_bidx[n];
  }
  
  thrust::copy(h_ids2.begin(), h_ids2.end(), d_ids);
  thrust::copy(h_bidx2.begin(), h_bidx2.end(), d_bidx);
}

void
cuda_mprts_sort_pairs_device(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

#if 0
  thrust::device_ptr<unsigned int> d_ids(mprts_cuda->d_ids);
  thrust::device_ptr<unsigned int> d_bidx(mprts_cuda->d_bidx);
  thrust::sequence(d_ids, d_ids + mprts_cuda->nr_prts);
  thrust::stable_sort_by_key(d_bidx, d_bidx + mprts_cuda->nr_prts, d_ids);
#else
  do_sort(mprts);
#endif
  // d_ids contains the indices to read from, d_bidx is sorted (for now)
}
