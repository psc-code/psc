
#ifndef CUDA_MFIELDS_H
#define CUDA_MFIELDS_H

#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "psc_fields_cuda.h"
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

using MfieldsStorageDeviceVector = MfieldsStorageVector<thrust::device_vector<float>>;

// ======================================================================
// MfieldsStorageDeviceRaw

class MfieldsStorageDeviceRaw
{
public:
  using value_type = float;
  using iterator = float*;
  using const_iterator = const float*;
  
  MfieldsStorageDeviceRaw(uint stride, value_type* data)
    : stride_{stride}, data_{data}
  {}
  
  KG_INLINE value_type* operator[](int p) { return data_ + p * stride_; }
  KG_INLINE const value_type* operator[](int p) const { return data_ + p * stride_; }

private:
  value_type *data_;
  uint stride_;
};

// ======================================================================
// DMfields

struct DMFields;

template <>
struct MfieldsCRTPInnerTypes<DMFields>
{
  using Storage = MfieldsStorageDeviceRaw;
};

struct DMFields : MfieldsCRTP<DMFields>
{
  using Base = MfieldsCRTP<DMFields>;
  using Storage = typename Base::Storage;
  using real_t = typename Base::Real;
  
  DMFields(const kg::Box3& box, int n_comps, int n_patches, real_t* d_flds)
    : Base{n_comps, box, n_patches},
      storage_{uint(n_comps * box.size()), d_flds}
  {}

private:
  Storage storage_;
  
  KG_INLINE Storage& storageImpl() { return storage_; }
  KG_INLINE const Storage& storageImpl() const { return storage_; }

  friend class MfieldsCRTP<DMFields>;
};

using DFields = DMFields::fields_view_t;

// ======================================================================
// cuda_mfields

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

  mrc_json_t to_json();
  void dump(const char *filename);

  real_t *data() { return storage_.data(); }
  operator DMFields();
  DFields operator[](int p) const; // FIXME, const correctness
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

  KG_INLINE Storage& storageImpl() { return storage_; }
  KG_INLINE const Storage& storageImpl() const { return storage_; }

  friend class MfieldsCRTP<cuda_mfields>;
};

HMFields hostMirror(const cuda_mfields& cmflds);
void copy(const cuda_mfields& cmflds, HMFields& hmflds);
void copy(const HMFields& hmflds, cuda_mfields& cmflds);

cuda_mfields::fields_host_t get_host_fields(const cuda_mfields& cmflds);
void copy_to_device(int p, const cuda_mfields::fields_host_t& h_flds, cuda_mfields& cmflds, int mb, int me);
void copy_from_device(int p, cuda_mfields::fields_host_t& h_flds, const cuda_mfields& cmflds, int mb, int me);

#endif
