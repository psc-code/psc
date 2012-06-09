
#include "psc_cuda.h"
#include "particles_cuda.h"

#include <thrust/sort.h>
#include <thrust/device_vector.h>

void
cuda_mprts_sort_pairs_device(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  thrust::device_ptr<unsigned int> td_bidx(mprts_cuda->d_bidx);
  thrust::device_ptr<unsigned int> td_ids(mprts_cuda->d_ids);
  thrust::sequence(td_ids, td_ids + mprts_cuda->nr_prts);
  thrust::stable_sort_by_key(td_bidx, td_bidx + mprts_cuda->nr_prts, td_ids);
}
