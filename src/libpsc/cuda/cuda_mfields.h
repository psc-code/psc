
#ifndef CUDA_MFIELDS_H
#define CUDA_MFIELDS_H

#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "cuda_base.cuh"
#include "psc_fields_cuda.h"
#include "dim.hxx"

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
// MfieldsStorageDeviceRaw

class MfieldsStorageDeviceRaw
{
public:
  using value_type = float;
  using reference = float&;
  using const_reference = const float&;
  using iterator = float*;
  using const_iterator = const float*;
  
  MfieldsStorageDeviceRaw(uint stride, value_type* data)
    : stride_{stride}, data_{data}
  {}
  
  KG_INLINE reference operator[](size_t i) { return data_[i]; }
  KG_INLINE const_reference operator[](size_t i) const { return data_[i]; }

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

using MfieldsStorageDeviceVector = thrust::device_vector<float>;

struct cuda_mfields : CudaMfields<MfieldsStorageDeviceVector>
{
  using Base = CudaMfields<MfieldsStorageDeviceVector>;
  using Storage = typename Base::Storage;
  using real_t = typename Storage::value_type;
  using pointer = typename Storage::pointer;
  using const_pointer = typename Storage::const_pointer;

  cuda_mfields(const Grid_t& grid, int n_comps, const Int3& ibn)
    : Base{{-ibn, grid.ldims + 2*ibn}, n_comps, grid.n_patches()},
      grid_{grid}
  {
    cuda_base_init();
  }

  using Base::Base;
  cuda_mfields(const cuda_mfields&) = delete;

  void zero_comp(int m, dim_yz tag);
  void zero_comp(int m, dim_xyz tag);
  
  void copy_comp_yz(int m_to, cuda_mfields *from, int m_from);
  void copy_comp_xyz(int m_to, cuda_mfields *from, int m_from);
  void copy_comp(int m_to, cuda_mfields *from, int m_from);

  void axpy_comp_yz(int ym, float a, cuda_mfields *x, int xm);
  void axpy_comp_xyz(int ym, float a, cuda_mfields *x, int xm);
  void axpy_comp(int ym, float a, cuda_mfields *x, int xm);

  mrc_json_t to_json();
  void dump(const char *filename);

  pointer data() { return storage().data(); }
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

  real_t get_value(int idx) const { return storage()[idx]; }
  void set_value(int idx, real_t val) { storage()[idx] = val; }
  
private:
  const Grid_t& grid_;
};

HMFields hostMirror(const cuda_mfields& cmflds);
void copy(const cuda_mfields& cmflds, HMFields& hmflds);
void copy(const HMFields& hmflds, cuda_mfields& cmflds);

#endif
