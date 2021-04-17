
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

// ======================================================================
// cuda_mfields

using MfieldsStorageDeviceVector = gt::gtensor<float, 5, gt::space::device>;

struct cuda_mfields;

template <>
struct MfieldsCRTPInnerTypes<cuda_mfields>
{
  using Storage = gt::gtensor<float, 5, gt::space::device>;
};

struct cuda_mfields
{
  using Base = MfieldsCRTP<cuda_mfields>;
  using Storage = typename Base::Storage;
  using real_t = typename Storage::value_type;
  using pointer = typename Storage::pointer;
  using const_pointer = typename Storage::const_pointer;

  cuda_mfields(const Grid_t& grid, int n_comps, const Int3& ibn)
    : storage_({grid.ldims[0] + 2 * ibn[0], grid.ldims[1] + 2 * ibn[1],
                grid.ldims[2] + 2 * ibn[2], n_comps, grid.n_patches()})
  {}

  cuda_mfields(const cuda_mfields&) = delete;

  Storage& storage() { return storage_; }
  const Storage& storage() const { return storage_; }

private:
  Storage storage_;
};

#endif
