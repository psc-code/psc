
#include "psc_cuda.h"
#include "particles_cuda.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>

#define NBLOCKS_X 1
#define NBLOCKS_Y 8
#define NBLOCKS_Z 8

// spine
//     0   1   2   3   4   5   6   7   8   NEW
// b0 |   |   |   |   |   |   |   |   |   |   |
// b1 |   |   |   |   |   |   |   |   |   |   |
// b2 |   |   |   |   |   |   |   |   |   |   |
// ...
// bn |   |   |   |   |   |   |   |   |   |   |

//    |   |   |   |   |   |   |   |   |   |   |   |   | ... |   | # oob
//     b0  b1  b2  b3                                        bn

static unsigned int
bidx_to_key(unsigned int bidx, unsigned int bid, int nr_patches)
{
  if (bidx == NBLOCKS_X * NBLOCKS_Y * NBLOCKS_Z * nr_patches) {
    return CUDA_BND_S_OOB;
  } else {
    int b_diff = bid - bidx + NBLOCKS_Y + 1;
    int d1 = b_diff % NBLOCKS_Y;
    int d2 = b_diff / NBLOCKS_Y;
    return d2 * 3 + d1;
  }
}

void
cuda_mprts_spine_reduce(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  unsigned int nr_total_blocks = mprts_cuda->nr_total_blocks;
  unsigned int nr_blocks = mprts_cuda->nr_blocks;
  assert(nr_blocks == NBLOCKS_Y * NBLOCKS_Z);

  thrust::device_ptr<unsigned int> d_spine_cnts(mprts_cuda->d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);
  thrust::device_ptr<unsigned int> d_bidx(mprts_cuda->d_bidx);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);

  thrust::fill(d_spine_cnts, d_spine_cnts + 1 + nr_total_blocks * (CUDA_BND_STRIDE + 1), 0);

  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + nr_total_blocks + 1);
  thrust::host_vector<unsigned int> h_spine_cnts(d_spine_cnts, d_spine_cnts + 1 + nr_total_blocks * (CUDA_BND_STRIDE + 1));

  for (int p = 0; p < mprts->nr_patches; p++) {
    for (int b = 0; b < nr_blocks; b++) {
      unsigned int bid = b + p * nr_blocks;
      for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
	assert((h_bidx[n] >= p * nr_blocks && h_bidx[n] < (p+1) * nr_blocks) ||
	       (h_bidx[n] == nr_total_blocks));

	unsigned int key = bidx_to_key(h_bidx[n], bid, mprts->nr_patches);
	if (key < 9) {
	  int dy = key % 3;
	  int dz = key / 3;
	  int by = b % NBLOCKS_Y;
	  int bz = b / NBLOCKS_Y;
	  unsigned int bby = by + 1 - dy;
	  unsigned int bbz = bz + 1 - dz;
	  unsigned int bb = bbz * NBLOCKS_Y + bby;
	  if (bby < NBLOCKS_Y && bbz < NBLOCKS_Z) {
	    h_spine_cnts[(bb + p * nr_blocks) * 10 + key]++;
	  } else {
	    assert(0);
	  }
	} else if (key == CUDA_BND_S_OOB) {
	  h_spine_cnts[NBLOCKS_Y*NBLOCKS_Z*mprts->nr_patches * 10 + bid]++;
	}
      }
    }
  }  

  for (int p = 0; p < mprts->nr_patches; p++) {
    unsigned int n_oob = 0;
    for (int b = 0; b < nr_blocks; b++) {
      unsigned int bid = b + p * nr_blocks;
      n_oob += h_spine_cnts[NBLOCKS_Y*NBLOCKS_Z*mprts->nr_patches * 10 + bid];
    }
    mprintf("p%d: n_oob %d\n", p, n_oob);
  }

  thrust::copy(h_spine_cnts.begin(), h_spine_cnts.end(), d_spine_cnts);
  thrust::exclusive_scan(d_spine_cnts + nr_total_blocks * 10,
			 d_spine_cnts + nr_total_blocks * 10 + nr_total_blocks + 1,
			 d_spine_sums + nr_total_blocks * 10);
}

void
cuda_mprts_do_sort(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  int nr_blocks = mprts_cuda->nr_blocks;
  int nr_total_blocks = mprts_cuda->nr_total_blocks;

  thrust::device_ptr<unsigned int> d_bidx(mprts_cuda->d_bidx);
  thrust::device_ptr<unsigned int> d_ids(mprts_cuda->d_ids);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);
  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_ids(mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + nr_total_blocks + 1);

  thrust::host_vector<unsigned int> h_b_cnts(nr_total_blocks * (nr_total_blocks + 1));
  thrust::host_vector<unsigned int> h_b_sums(nr_total_blocks * (nr_total_blocks + 1));
  thrust::host_vector<unsigned int> h_bidx2(mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_spine_cnts(nr_total_blocks * (10 + 1));
  thrust::host_vector<unsigned int> h_spine_sums(nr_total_blocks * (10 + 1));

  for (int bid = 0; bid < nr_total_blocks; bid++) {
    int b = bid % nr_blocks;
    int p = bid / nr_blocks;
    for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
      unsigned int key = bidx_to_key(h_bidx[n], bid, mprts->nr_patches);
      if (key < 9) {
	h_b_cnts[h_bidx[n] * (nr_total_blocks + 1) + bid]++;
	
	int dy = key % 3;
	int dz = key / 3;
	int by = b % NBLOCKS_Y;
	int bz = b / NBLOCKS_Y;
	unsigned int bby = by + 1 - dy;
	unsigned int bbz = bz + 1 - dz;
	assert(bby < NBLOCKS_Y && bbz < NBLOCKS_Z);
	unsigned int bb = bbz * NBLOCKS_Y + bby;
	h_spine_cnts[(bb + p * nr_blocks) * 10 + key]++;
      } else { // OOB
	h_spine_cnts[nr_total_blocks * 10 + bid]++;
      }
    }
  }
  for (int n = mprts_cuda->nr_prts - mprts_cuda->nr_prts_recv; n < mprts_cuda->nr_prts; n++) {
    assert(h_bidx[n] < mprts_cuda->nr_total_blocks);
    h_b_cnts[h_bidx[n] * (nr_total_blocks + 1) + nr_total_blocks]++;
    h_spine_cnts[h_bidx[n] * 10 + CUDA_BND_S_NEW]++;
  }

#if 0
  for (int i = 0; i < 200; i++) {
    if (h_b_cnts[i] != 0) printf("b %d: %d\n", i, h_b_cnts[i]);
  }
  for (int i = 0; i < 200; i++) {
    if (h_spine_cnts[i] != 0) printf("s %d: %d\n", i, h_spine_cnts[i]);
  }
#endif

  thrust::exclusive_scan(h_b_cnts.begin(), h_b_cnts.end(), h_b_sums.begin());
  thrust::exclusive_scan(h_spine_cnts.begin(), h_spine_cnts.end(), h_spine_sums.begin());

  for (int bid = 0; bid < nr_total_blocks; bid++) {
    int b = bid % nr_blocks;
    int p = bid / nr_blocks;
    for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
      unsigned int key = bidx_to_key(h_bidx[n], bid, mprts->nr_patches);
      if (key < 9) {
	h_b_cnts[h_bidx[n] * (nr_total_blocks + 1) + bid]++;
	
	int dy = key % 3;
	int dz = key / 3;
	int by = b % NBLOCKS_Y;
	int bz = b / NBLOCKS_Y;
	unsigned int bby = by + 1 - dy;
	unsigned int bbz = bz + 1 - dz;
	assert(bby < NBLOCKS_Y && bbz < NBLOCKS_Z);
	unsigned int bb = bbz * NBLOCKS_Y + bby;
	int nn = h_b_sums[h_bidx[n] * (nr_total_blocks + 1) + bid]++;
	int nn2 = h_spine_sums[(bb + p * nr_blocks) * 10 + key]++;
	if (nn != nn2) {
	  mprintf("n %d nn %d nn2 %d\n", n, nn, nn2);
	}
	assert(nn == nn2);
	h_ids[nn] = n;
	h_bidx2[nn] = h_bidx[n];
      } else { // OOB
      }
    }
  }
  for (int n = mprts_cuda->nr_prts - mprts_cuda->nr_prts_recv; n < mprts_cuda->nr_prts; n++) {
      int nn = h_b_sums[h_bidx[n] * (nr_total_blocks + 1) + nr_total_blocks]++;
      int nn2 = h_spine_sums[h_bidx[n] * 10 + CUDA_BND_S_NEW]++;
      assert(nn == nn2);
      h_ids[nn] = n;
      h_bidx2[nn] = h_bidx[n];
  }

  thrust::copy(h_ids.begin(), h_ids.end(), d_ids);
  thrust::copy(h_bidx2.begin(), h_bidx2.end(), d_bidx);
}

