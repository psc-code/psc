
#ifndef CUDA_MFIELDS_H
#define CUDA_MFIELDS_H

#include "cuda_bits.h"
#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "cuda_base.cuh"
#include "dim.hxx"

#define MAX_BND_FIELDS (17)
#define MAX_BND_COMPONENTS (3)

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

// ----------------------------------------------------------------------
// cuda_mfields_bnd_map
//
// the maps differs depending on the number of components in the field
// (incorporates specific strides)

struct cuda_mfields_bnd_map
{
  ~cuda_mfields_bnd_map() {}

  thrust::host_vector<int> h_map_out; // maps thread id to a particular offset
                                      // for ghosts in the flds array
  psc::device_vector<int> d_map_out;
  thrust::host_vector<int> h_map_in; // maps thread id to a particular offset
                                     // for ghosts in the flds array
  psc::device_vector<int> d_map_in;
};

// ----------------------------------------------------------------------
// cuda_mfields_bnd

struct cuda_mfields_bnd
{
  int n_patches;
  int im[3];
  int ib[3];
  struct cuda_mfields_bnd_patch* bnd_by_patch;
  psc::device_vector<fields_cuda_real_t> d_buf;
  thrust::host_vector<fields_cuda_real_t> h_buf;
  thrust::host_vector<int> h_nei_patch;
  psc::device_vector<int> d_nei_patch;
  struct cuda_mfields_bnd_map map[MAX_BND_FIELDS];
};

// ======================================================================
// DMfields

struct DMFields;

template <>
struct MfieldsCRTPInnerTypes<DMFields>
{
  using Storage = gt::gtensor_span_device<float, 5>;
};

struct DMFields : MfieldsCRTP<DMFields>
{
  using Base = MfieldsCRTP<DMFields>;
  using Storage = typename Base::Storage;
  using real_t = typename Base::Real;

  DMFields(const kg::Box3& box, int n_comps, int n_patches, real_t* d_flds)
    : Base{n_comps, box, n_patches},
      storage_{gt::adapt_device(
        d_flds, gt::shape(box.im(0), box.im(1), box.im(2), n_comps, n_patches))}
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

using MfieldsStorageDeviceVector = gt::gtensor<float, 5, gt::space::device>;

struct cuda_mfields;

template <>
struct MfieldsCRTPInnerTypes<cuda_mfields>
{
  using Storage = gt::gtensor<float, 5, gt::space::device>;
};

struct cuda_mfields : MfieldsCRTP<cuda_mfields>
{
  using Base = MfieldsCRTP<cuda_mfields>;
  using Storage = typename Base::Storage;
  using real_t = typename Storage::value_type;
  using pointer = typename Storage::pointer;
  using const_pointer = typename Storage::const_pointer;

  cuda_mfields(const Grid_t& grid, int n_comps, const Int3& ibn)
    : Base{n_comps, {-ibn, grid.ldims + 2 * ibn}, grid.n_patches()},
      storage_(
        {box().im(0), box().im(1), box().im(2), n_comps, grid.n_patches()}),
      grid_{grid}
  {
    cuda_base_init();
  }

  using Base::Base;
  cuda_mfields(const cuda_mfields&) = delete;

  mrc_json_t to_json();
  void dump(const char* filename);

  pointer data() { return storage().data(); }
  operator DMFields();
  DFields operator[](int p) const; // FIXME, const correctness
  const Grid_t& grid() const { return grid_; }

  int index(int m, int i, int j, int k, int p) const
  {
    return (
      ((((p)*n_comps() + m) * im(2) + (k - ib(2))) * im(1) + (j - ib(1))) *
        im(0) +
      (i - ib(0)));
  }

  real_t get_value(int idx) const { return storage().data()[idx]; }
  void set_value(int idx, real_t val) { storage().data()[idx] = val; }

private:
  const Grid_t& grid_;
  Storage storage_;

  KG_INLINE Storage& storageImpl() { return storage_; }
  KG_INLINE const Storage& storageImpl() const { return storage_; }

  friend class MfieldsCRTP<cuda_mfields>;
};

MfieldsSingle hostMirror(cuda_mfields& cmflds);
MfieldsSingle hostMirror(const cuda_mfields& cmflds);
void copy(const cuda_mfields& cmflds, MfieldsSingle& hmflds);
void copy(const MfieldsSingle& hmflds, cuda_mfields& cmflds);

#endif
