
#ifndef CUDA_MFIELDS_H
#define CUDA_MFIELDS_H

#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "dim.hxx"

#include <kg/SArray.h>

#define MAX_BND_FIELDS (17)
#define MAX_BND_COMPONENTS (3)

#include <thrust/device_vector.h>

// ----------------------------------------------------------------------
// cuda_mfields_bnd_map
//
// the maps differs depending on the number of components in the field
// (incorporates specific strides)

struct cuda_mfields_bnd_map {
  int *h_map_out; // maps thread id to a particular offset for ghosts in the flds array 
  int *d_map_out;
  int nr_map_out; // number of entries in the map
  int *h_map_in; // maps thread id to a particular offset for ghosts in the flds array 
  int *d_map_in;
  int nr_map_in; // number of entries in the map
};

// ----------------------------------------------------------------------
// cuda_mfields_bnd

struct cuda_mfields_bnd {
  int n_patches;
  int im[3];
  int ib[3];
  struct cuda_mfields_bnd_patch *bnd_by_patch;
  thrust::device_vector<fields_cuda_real_t> d_buf;
  thrust::host_vector<fields_cuda_real_t> h_buf;
  int *h_nei_patch;
  int *d_nei_patch;
  struct cuda_mfields_bnd_map map[MAX_BND_FIELDS];
};

// ======================================================================
// cuda_mfields

struct DMFields;
struct DFields;

struct cuda_mfields
{
  using real_t = float;
  using fields_host_t = kg::SArray<real_t>;

  cuda_mfields(const Grid_t& grid, int n_fields, const Int3& ibn);
  cuda_mfields(const cuda_mfields&) = delete;

  void zero_comp(int m, dim_yz tag);
  void zero_comp(int m, dim_xyz tag);
  
  void axpy_comp_yz(int ym, float a, cuda_mfields *x, int xm);

  fields_host_t get_host_fields();
  void copy_to_device(int p, const fields_host_t& h_flds, int mb, int me);
  void copy_from_device(int p, fields_host_t& h_flds, int mb, int me);

  mrc_json_t to_json();
  void dump(const char *filename);

  real_t *data() { return d_flds_.data().get(); }
  operator DMFields();
  DFields operator[](int p);
  const Grid_t& grid() const { return grid_; }

  int index(int m, int i, int j, int k, int p) const
  {
    return (((((p)
	       * n_fields + m)
	      * im[2] + (k - ib[2]))
	     * im[1] + (j - ib[1]))
	    * im[0] + (i - ib[0]));
  }

  real_t get_value(int idx) const { return d_flds_[idx]; }
  void set_value(int idx, real_t val) { d_flds_[idx] = val; }
  
public:
  Int3 ib;
  Int3 im;
  int n_patches;
  int n_fields;
  int n_cells_per_patch;
  int n_cells;
private:
  thrust::device_vector<fields_cuda_real_t> d_flds_;
  const Grid_t& grid_;
};

// ======================================================================
// DFields

struct DFields
{
  using Storage = kg::StorageNoOwnership<float>;
  using real_t = typename Storage::value_type;
  
  __host__ __device__ DFields(real_t* d_flds, Int3 im, Int3 ib, int n_comps)
    : storage_{d_flds},
      im_{im},
      ib_{ib},
      n_comps_{n_comps}
  {}
  
  __device__ real_t  operator()(int m, int i, int j, int k) const { return storage_[index(m, {i,j,k})]; }
  __device__ real_t& operator()(int m, int i, int j, int k)       { return storage_[index(m, {i,j,k})]; }

  __host__ real_t *data() { return storage_.data(); }

  __device__ int im(int d) const { return im_[d]; }

private:
  __device__ int index(int m, Int3 idx) const
  {
#if 0
    if (idx[0] - ib_[0] < 0 || idx[0] - ib_[0] >= im_[0]) printf("!!! i %d\n", j);
    if (jdx[1] - ib_[1] < 0 || idx[1] - ib_[1] >= im_[1]) printf("!!! j %d\n", j);
    if (idx[2] - ib_[2] < 0 || idx[2] - ib_[2] >= im_[2]) printf("!!! k %d\n", k);
#endif
    return kg::layoutDataOffset<kg::LayoutSOA>(n_comps_, im_, m, idx);
  }

private:
  Storage storage_;
  Int3 im_;
  Int3 ib_;
  int n_comps_;
};

// ======================================================================
// DMFields

struct DMFields
{
  using real_t = float;
  
  __host__ DMFields(real_t* d_flds, uint stride, Int3 im, Int3 ib, int n_comps)
    : d_flds_(d_flds),
      stride_(stride),
      im_{im},
      ib_{ib},
      n_comps_{n_comps}
  {}

  __host__ __device__ DFields operator[](int p)
  {
    return DFields(d_flds_ + p * stride_, im_, ib_, n_comps_);
  }

  __device__ int im(int d) const { return im_[d]; }
  
private:
  real_t *d_flds_;
  uint stride_;
  Int3 im_;
  Int3 ib_;
  int n_comps_;
};

#endif
