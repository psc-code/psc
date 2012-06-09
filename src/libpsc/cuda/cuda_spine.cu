
#include "psc_cuda.h"
#include "particles_cuda.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

void
cuda_mprts_spine_reduce(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  int nr_blocks;
  assert(mprts->nr_patches > 0);
  {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, 0);
    struct psc_particles_cuda *prts_cuda = psc_particles_cuda(prts);
    nr_blocks = prts_cuda->nr_blocks; // FIXME...
  }
  unsigned int nr_total_blocks = nr_blocks * mprts->nr_patches;

  thrust::device_ptr<unsigned int> d_spine_cnts(mprts_cuda->d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_bidx(mprts_cuda->d_bidx);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_bidx);

  thrust::fill(d_spine_cnts, d_spine_cnts + nr_total_blocks * (CUDA_BND_STRIDE + 1), 0);
  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + nr_total_blocks + 1);

  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    for (int n = 0; n < prts->n_part; n++) {
      assert((h_bidx[off + n] >= p * nr_blocks && h_bidx[off + n] < (p+1) * nr_blocks) ||
	     (h_bidx[off + n] == nr_total_blocks));
    }
    off += prts->n_part;
  }  
}

