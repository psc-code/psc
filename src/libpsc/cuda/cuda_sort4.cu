
#include "psc_cuda.h"
#include "particles_cuda.h"

#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/device_vector.h>

void
cuda_mprts_sort_pairs_device(struct psc_mparticles *mprts)
{
#if 0
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  thrust::device_ptr<unsigned int> d_ids(mprts_cuda->d_ids);
  thrust::device_ptr<unsigned int> d_bidx(mprts_cuda->d_bidx);
  thrust::sequence(d_ids, d_ids + mprts_cuda->nr_prts);
  thrust::stable_sort_by_key(d_bidx, d_bidx + mprts_cuda->nr_prts, d_ids);
#else
  cuda_mprts_do_sort(mprts);
#endif
  // d_ids contains the indices to read from, d_bidx is sorted (for now)
}
