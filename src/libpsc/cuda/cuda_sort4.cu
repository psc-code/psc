
#include "psc_cuda.h"
#include "particles_cuda.h"

#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/device_vector.h>

static void
do_sort(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  int nr_total_blocks = mprts_cuda->nr_total_blocks;

  thrust::device_ptr<unsigned int> d_bidx(mprts_cuda->d_bidx);
  thrust::device_ptr<unsigned int> d_ids(mprts_cuda->d_ids);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);
  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_ids(mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + nr_total_blocks + 1);

  thrust::host_vector<unsigned int> h_b_cnts(nr_total_blocks * (nr_total_blocks + 1));
  thrust::host_vector<unsigned int> h_cnts(nr_total_blocks + 1);

  for (int bid = 0; bid < nr_total_blocks; bid++) {
    for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
      if (h_bidx[n] < mprts_cuda->nr_total_blocks) {
	h_cnts[h_bidx[n]]++;
	h_b_cnts[h_bidx[n] * (nr_total_blocks + 1) + bid]++;
      }
    }
  }
  for (int n = mprts_cuda->nr_prts - mprts_cuda->nr_prts_recv; n < mprts_cuda->nr_prts; n++) {
    assert(h_bidx[n] < mprts_cuda->nr_total_blocks);
    h_cnts[h_bidx[n]]++;
    h_b_cnts[h_bidx[n] * (nr_total_blocks + 1) + nr_total_blocks]++;
  }

  thrust::host_vector<unsigned int> h_sums(nr_total_blocks + 1);
  thrust::host_vector<unsigned int> h_b_sums(nr_total_blocks * (nr_total_blocks + 1));

  thrust::exclusive_scan(h_cnts.begin(), h_cnts.end(), h_sums.begin());
  thrust::exclusive_scan(h_b_cnts.begin(), h_b_cnts.end(), h_b_sums.begin());

  thrust::host_vector<unsigned int> h_bidx2(mprts_cuda->nr_prts);

  mprintf("h_sums %d // %d %d\n", h_sums[1], h_b_sums[1], h_b_sums[nr_total_blocks + 1]);
  for (int bid = 0; bid < nr_total_blocks; bid++) {
    for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
      if (h_bidx[n] < nr_total_blocks) {
	int nn = h_sums[h_bidx[n]]++;
	int nn2 = h_b_sums[h_bidx[n] * (nr_total_blocks + 1) + bid]++;
	if (nn != nn2) {
	  mprintf("bid %d nn %d nn2 %d\n", bid, nn, nn2);
	}
	assert(nn == nn2);
	h_ids[nn] = n;
	h_bidx2[nn] = h_bidx[n];
      }
    }
  }
  for (int n = mprts_cuda->nr_prts - mprts_cuda->nr_prts_recv; n < mprts_cuda->nr_prts; n++) {
      int nn = h_sums[h_bidx[n]]++;
      int nn2 = h_b_sums[h_bidx[n] * (nr_total_blocks + 1) + nr_total_blocks]++;
      assert(nn == nn2);
      h_ids[nn] = n;
      h_bidx2[nn] = h_bidx[n];
  }

  thrust::copy(h_ids.begin(), h_ids.end(), d_ids);
  thrust::copy(h_bidx2.begin(), h_bidx2.end(), d_bidx);
}

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
  do_sort(mprts);
#endif
  // d_ids contains the indices to read from, d_bidx is sorted (for now)
}
