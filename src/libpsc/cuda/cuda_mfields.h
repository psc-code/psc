
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

  cuda_mfields(const cuda_mfields&) = delete;

  pointer data() { return storage().data(); }
  const Grid_t& grid() const { return grid_; }

private:
  const Grid_t& grid_;
  Storage storage_;

  KG_INLINE Storage& storageImpl() { return storage_; }
  KG_INLINE const Storage& storageImpl() const { return storage_; }

  friend class MfieldsCRTP<cuda_mfields>;
};

#endif
