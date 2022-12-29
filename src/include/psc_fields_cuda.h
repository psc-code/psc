
#ifndef PSC_FIELDS_CUDA_H
#define PSC_FIELDS_CUDA_H

#include <mpi.h>
#include "cuda_bits.h"
#include "fields3d.hxx"
#include "fields_traits.hxx"

#include "psc_fields_single.h"

#include "mrc_json.h"

// ======================================================================
// MfieldsCuda

namespace psc
{
template <typename T, gt::size_type N>
using gtensor_device = gt::gtensor_container<psc::device_vector<T>, N>;
}

struct MfieldsCuda;

template <>
struct MfieldsCRTPInnerTypes<MfieldsCuda>
{
  using Storage = psc::gtensor_device<float, 5>;
};

struct MfieldsCuda
  : MfieldsBase
  , MfieldsCRTP<MfieldsCuda>
{
  using real_t = float;
  using Real = real_t;
  using Base = MfieldsCRTP<MfieldsCuda>;
  using Storage = typename Base::Storage;
  using space = gt::space::device;

  MfieldsCuda(const Grid_t& grid, int n_fields, Int3 ibn)
    : MfieldsBase{grid, n_fields, ibn},
      Base(n_fields, {-ibn, grid.ldims + 2 * ibn}, grid.n_patches()),
      storage_(gt::shape(Base::box().im(0), Base::box().im(1),
                         Base::box().im(2), n_fields, Base::n_patches()))
  {
    thrust::fill(storage_.data(), storage_.data() + storage_.size(), real_t{});
  }

  const Grid_t& grid() const { return _grid(); }

  void reset(const Grid_t& new_grid) override
  {
    *this = MfieldsCuda(new_grid, n_comps(), ibn());
  }

  auto gt()
  {
    return gt::adapt_device(storage_.data().get(), storage_.shape());
  }

  auto gt() const
  {
    return gt::adapt_device(storage_.data().get(), storage_.shape());
  }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

private:
  Storage storage_;

  Storage& storageImpl() { return storage_; }
  const Storage& storageImpl() const { return storage_; }

  friend class MfieldsCRTP<MfieldsCuda>;
};

inline MfieldsSingle hostMirror(MfieldsCuda& mflds)
{
  return MfieldsSingle{mflds.grid(), mflds.n_comps(), mflds.ibn()};
}

inline MfieldsSingle hostMirror(const MfieldsCuda& mflds)
{
  return MfieldsSingle{mflds.grid(), mflds.n_comps(), mflds.ibn()};
}

inline void copy(const MfieldsCuda& mflds, MfieldsSingle& hmflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("mflds to host", 1., 0, 0);
  }
  prof_start(pr);
  // gt::copy(mflds.gt(), hmflds.storage());
  thrust::copy(mflds.gt().data(), mflds.gt().data() + mflds.gt().size(),
               hmflds.storage().data());
  prof_stop(pr);
}

inline void copy(const MfieldsSingle& hmflds, MfieldsCuda& mflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("mflds from host", 1., 0, 0);
  }
  prof_start(pr);
  // gt::copy(hmflds.storage(), mflds.gt());
  thrust::copy(hmflds.storage().data(),
               hmflds.storage().data() + hmflds.storage().size(),
               mflds.gt().data());
  prof_stop(pr);
}

// ======================================================================
// MfieldsStateCuda

struct MfieldsStateCuda : MfieldsStateBase
{
  using real_t = MfieldsCuda::real_t;
  using Storage = MfieldsCuda::Storage;
  using space = gt::space::device;

  MfieldsStateCuda(const Grid_t& grid)
    : MfieldsStateBase{grid, NR_FIELDS, grid.ibn},
      mflds_{grid, NR_FIELDS, grid.ibn}
  {}

  int n_patches() const { return mflds_.n_patches(); }
  int n_comps() const { return mflds_.n_comps(); }
  const Int3& ib() const { return mflds_.ib(); }

  const Grid_t& grid() const { return *grid_; }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  MfieldsCuda& mflds() { return mflds_; }
  const MfieldsCuda& mflds() const { return mflds_; }

  auto gt() { return mflds_.gt(); }
  auto& storage() { return mflds_.storage(); }
  auto& storage() const { return mflds_.storage(); }

private:
  MfieldsCuda mflds_;
};

template <>
struct Mfields_traits<MfieldsCuda>
{
  static constexpr const char* name = "cuda";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

namespace detail
{
template <typename T>
struct Mfields_from_type_space<T, gt::space::device>
{
  static_assert(std::is_same<T, float>::value, "CUDA only supports float");
  using type = MfieldsCuda;
};
} // namespace detail

inline MfieldsSingle hostMirror(MfieldsStateCuda& mflds)
{
  return hostMirror(mflds.mflds());
}

inline MfieldsSingle hostMirror(const MfieldsStateCuda& mflds)
{
  return hostMirror(mflds.mflds());
}

inline void copy(const MfieldsStateCuda& mflds, MfieldsSingle& hmflds)
{
  copy(mflds.mflds(), hmflds);
}

inline void copy(const MfieldsSingle& hmflds, MfieldsStateCuda& mflds)
{
  copy(hmflds, mflds.mflds());
}

// ----------------------------------------------------------------------
// FIXME hacky workaround for lack of gt::view on device

template <typename E>
GT_INLINE auto view_patch(E&& gt, int p)
{
  return gt::adapt_device<4>(
    (&gt(0, 0, 0, 0, p)).get(),
    {gt.shape(0), gt.shape(1), gt.shape(2), gt.shape(3)});
}

#endif
