
#pragma once

// ======================================================================
// cuda_mparticles_sort

struct cuda_mparticles_sort
{
  cuda_mparticles_sort(uint n_cells)
    : d_off(n_cells + 1)
  {}

  template<typename BS>
  void find_indices_ids(cuda_mparticles<BS>& cmprts)
  {
    d_cidx.resize(cmprts.n_prts);
    d_id.resize(cmprts.n_prts);
    cmprts.find_cell_indices_ids(d_cidx, d_id);
  }
  
  void stable_sort_cidx()
  {
    thrust::stable_sort_by_key(d_cidx.begin(), d_cidx.end(), d_id.begin());
  }

  void find_offsets()
  {
    int n_cells = d_off.size() - 1;
    thrust::counting_iterator<uint> search_begin(0);
    thrust::upper_bound(d_cidx.begin(), d_cidx.end(),
			search_begin, search_begin + n_cells,
			d_off.begin() + 1);
    // d_off[0] was set to zero during d_off initialization
  }
  
  template<typename BS>
  void reorder(cuda_mparticles<BS>& cmprts)
  {
    cmprts.reorder(d_id);
  }
  
  template<typename BS>
  void reorder_and_offsets(cuda_mparticles<BS>& cmprts)
  {
    cmprts.reorder_and_offsets_cidx(d_cidx, d_id, d_off);
  }
  
public:
  thrust::device_vector<uint> d_cidx;     // cell index (incl patch) per particle
  thrust::device_vector<uint> d_id;       // particle id used for reordering
  thrust::device_vector<uint> d_off;      // particles per cell
                                          // are at indices [offsets[cell] .. offsets[cell+1][
};


