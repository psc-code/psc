
// ----------------------------------------------------------------------
// reorder_and_offsets_slow

template <typename BS>
void cuda_mparticles<BS>::reorder_and_offsets_slow()
{
  if (this->n_patches() == 0) {
    return;
  }

  HMparticlesCudaStorage h_storage(this->storage);
  HMparticlesCudaStorage h_alt_storage(this->alt_storage);
  thrust::host_vector<uint> h_off(this->by_block_.d_off);
  thrust::host_vector<uint> h_bidx(this->by_block_.d_idx);
  thrust::host_vector<uint> h_id(this->by_block_.d_id);

  for (int i = 0; i <= this->n_prts; i++) {
    //    uint bidx;
    uint block;
    if (i < this->n_prts) {
      h_alt_storage.store(h_storage[h_id[i]], i);
      // bidx = blockIndex(h_alt_storage[i], 0);
      block = h_bidx[i];
    } else {
      // bidx = n_blocks;
      block = this->n_blocks;
    }
    // if (i < 10) {
    //   printf("i %d bidx %d block %d x %g %g\n", bidx, block,
    //   h_alt_storage[i].x[1], h_alt_storage[i].x[2]);
    // }
    int prev_block = (i > 0) ? (int)h_bidx[i - 1] : -1;
    for (int b = prev_block + 1; b <= block; b++) {
      h_off[b] = i;
    }
  }

  h_alt_storage.resize(this->n_prts);
  thrust::copy(h_alt_storage.begin(), h_alt_storage.end(), alt_storage.begin());
  thrust::copy(h_off.begin(), h_off.end(), this->by_block_.d_off.begin());

  swap_alt();
  need_reorder = false;
}
