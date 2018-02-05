
// ----------------------------------------------------------------------
// reorder_and_offsets_slow

void cuda_mparticles::reorder_and_offsets_slow()
{
  if (n_patches == 0) {
    return;
  }

  thrust::host_vector<float4> h_xi4(d_xi4);
  thrust::host_vector<float4> h_pxi4(d_pxi4);
  thrust::host_vector<float4> h_alt_xi4(n_prts);
  thrust::host_vector<float4> h_alt_pxi4(n_prts);
  thrust::host_vector<uint> h_off(d_off);
  thrust::host_vector<uint> h_bidx(d_bidx);
  thrust::host_vector<uint> h_id(d_id);

  for (int i = 0; i <= n_prts; i++) {
    //    uint bidx;
    uint block;
    if (i < n_prts) {
      h_alt_xi4[i] = h_xi4[h_id[i]];
      h_alt_pxi4[i] = h_pxi4[h_id[i]];
      //bidx = blockIndex(h_alt_xi4[i], 0);
      block = h_bidx[i];
    } else {
      //bidx = n_blocks;
      block = n_blocks;
    }
    // if (i < 10) {
    //   printf("i %d bidx %d block %d xi4 %g %g\n", bidx, block, h_alt_xi4[i].y, h_alt_xi4[i].z);
    // }
    int prev_block = (i > 0) ? (int) h_bidx[i-1] : -1;
    for (int b = prev_block + 1; b <= block; b++) {
      h_off[b] = i;
    }
  }

  d_alt_xi4.resize(n_prts);
  d_alt_pxi4.resize(n_prts);
  thrust::copy(h_alt_xi4.begin(), h_alt_xi4.end(), d_alt_xi4.begin());
  thrust::copy(h_alt_pxi4.begin(), h_alt_pxi4.end(), d_alt_pxi4.begin());
  thrust::copy(h_off.begin(), h_off.end(), d_off.begin());
  
  swap_alt();
  need_reorder = false;
}

