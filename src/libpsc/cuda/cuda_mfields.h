
#ifndef CUDA_MFIELDS_H
#define CUDA_MFIELDS_H

#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "dim.hxx"

#include <kg/SArray.h>
#include <kg/SArrayView.h>

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
// MfieldsStorageDeviceVector

class MfieldsStorageDeviceVector
{
public:
  using value_type = float;
  
  MfieldsStorageDeviceVector(size_t size, uint stride)
    : d_flds_(size), stride_{stride}
  {}

  value_type* data() { return d_flds_.data().get(); }

  void set_value(int idx, const value_type& val) { d_flds_[idx] = val; }
  value_type get_value(int idx) const { return d_flds_[idx]; }

private:
  thrust::device_vector<value_type> d_flds_;
  uint stride_;
};

// ======================================================================
// cuda_mfields

struct DMFields;
using DFields = kg::SArrayView<float, kg::LayoutSOA>;

struct cuda_mfields;

template<>
struct MfieldsCRTPInnerTypes<cuda_mfields>
{
  using Storage = MfieldsStorageDeviceVector;
};

struct cuda_mfields : MfieldsCRTP<cuda_mfields>
{
  using Base = MfieldsCRTP<cuda_mfields>;
  using Storage = typename Base::Storage;
  using real_t = typename Storage::value_type;
  using fields_host_t = kg::SArray<real_t>;

  cuda_mfields(const Grid_t& grid, int n_fields, const Int3& ibn);
  cuda_mfields(const cuda_mfields&) = delete;

  void zero_comp(int m, dim_yz tag);
  void zero_comp(int m, dim_xyz tag);
  
  void axpy_comp_yz(int ym, float a, cuda_mfields *x, int xm);

  void copy_to_device(int p, const fields_host_t& h_flds, int mb, int me);
  void copy_from_device(int p, fields_host_t& h_flds, int mb, int me);

  mrc_json_t to_json();
  void dump(const char *filename);

  real_t *data() { return storage_.data(); }
  operator DMFields();
  DFields operator[](int p);
  const Grid_t& grid() const { return grid_; }

  int index(int m, int i, int j, int k, int p) const
  {
    return (((((p)
	       * n_comps() + m)
	      * im(2) + (k - ib(2)))
	     * im(1) + (j - ib(1)))
	    * im(0) + (i - ib(0)));
  }

  real_t get_value(int idx) const { return storage_.get_value(idx); }
  void set_value(int idx, real_t val) { storage_.set_value(idx, val); }
  
private:
  Storage storage_;
  const Grid_t& grid_;
};

cuda_mfields::fields_host_t get_host_fields(const cuda_mfields& cmflds);

// ======================================================================
// DMFieldsStorage

class DMFieldsStorage
{
public:
  using value_type = float;
  
  DMFieldsStorage(uint stride, value_type* data)
    : stride_{stride}, data_{data}
  {}
  
  KG_INLINE value_type* operator[](int p) { return data_ + p * stride_; }
  KG_INLINE const value_type* operator[](int p) const { return data_ + p * stride_; }

private:
  value_type *data_;
  uint stride_;
};

// ======================================================================
// DMFields

struct DMfields;

template <>
struct MfieldsCRTPInnerTypes<DMFields>
{
  using Storage = DMFieldsStorage;
};

struct DMFields : MfieldsCRTP<DMFields>
{
  using Base = MfieldsCRTP<DMFields>;

  using real_t = typename Base::Real;
  
  DMFields(real_t* d_flds, uint stride, Int3 im, Int3 ib, int n_comps, int n_patches)
    : Base{n_comps, ib, im, n_patches},
      storage_{stride, d_flds}
  {}

private:
  Storage storage_;
  
  KG_INLINE Storage& storageImpl() { return storage_; }
  KG_INLINE const Storage& storageImpl() const { return storage_; }

  friend class MfieldsCRTP<DMFields>;
};

#endif
