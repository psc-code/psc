
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
// StoragNoOwnershipDevice
//
// FIXME, same as StorageNoOwnership


template <typename T>
class StorageNoOwnershipDevice
{
public:
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;

  __host__ __device__ StorageNoOwnershipDevice(pointer data) : data_{data} {}

  __device__ const_reference operator[](int offset) const { return data_[offset]; }
  __device__ reference operator[](int offset) { return data_[offset]; }

  // FIXME access to underlying storage might better be avoided?
  // use of this makes assumption that storage is contiguous
  const_pointer data() const { return data_; }
  pointer data() { return data_; }

private:
  pointer data_;
};

// ======================================================================
// DFields

struct DFields
{
  using Storage = StorageNoOwnershipDevice<float>;
  using real_t = typename Storage::value_type;
  
  __host__ __device__ DFields(real_t* d_flds, Int3 im, Int3 ib)
    : storage_{d_flds},
      im_{im},
      ib_{ib}
  {}
  
  __device__ real_t  operator()(int m, int i, int j, int k) const { return storage_[index(m, i,j,k)]; }
  __device__ real_t& operator()(int m, int i, int j, int k)       { return storage_[index(m, i,j,k)]; }

  __host__ real_t *data() { return storage_.data(); }

  __device__ int im(int d) const { return im_[d]; }

private:
  __device__ int index(int m, int i, int j, int k) const
  {
#if 0
    if (i - ib_[0] < 0 || i - ib_[0] >= im_[0]) printf("!!! i %d\n", j);
    if (j - ib_[1] < 0 || j - ib_[1] >= im_[1]) printf("!!! j %d\n", j);
    if (k - ib_[2] < 0 || k - ib_[2] >= im_[2]) printf("!!! k %d\n", k);
#endif
    return ((((m)
	      *im_[2] + (k - ib_[2]))
	     *im_[1] + (j - ib_[1]))
	    *im_[0] + (i - ib_[0]));
  }

private:
  Storage storage_;
  Int3 im_;
  Int3 ib_;
};

// ======================================================================
// DMFields

struct DMFields
{
  using real_t = float;
  
  __host__ DMFields(real_t* d_flds, uint stride, Int3 im, Int3 ib)
    : d_flds_(d_flds),
      stride_(stride),
      im_{im},
      ib_{ib}
  {}

  __host__ __device__ DFields operator[](int p)
  {
    return DFields(d_flds_ + p * stride_, im_, ib_);
  }

  __device__ int im(int d) const { return im_[d]; }
  
private:
  real_t *d_flds_;
  uint stride_;
  Int3 im_;
  Int3 ib_;
};

#endif
