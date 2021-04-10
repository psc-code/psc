
#ifndef PSC_FIELDS_CUDA_H
#define PSC_FIELDS_CUDA_H

#include <mpi.h>
#include "fields3d.hxx"
#include "fields_traits.hxx"

#include "psc_fields_single.h"

#include "mrc_json.h"

struct cuda_mfields;

struct fields_cuda_t
{
  using real_t = float;
};

// ======================================================================
// CudaMfields

template <typename S>
struct CudaMfields;

template <typename S>
struct MfieldsCRTPInnerTypes<CudaMfields<S>>
{
  using Storage = S;
};

template <typename S>
struct CudaMfields : MfieldsCRTP<CudaMfields<S>>
{
  using Base = MfieldsCRTP<CudaMfields<S>>;
  using Storage = typename Base::Storage;
  using real_t = typename Base::Real;

  CudaMfields(const kg::Box3& box, int n_comps, int n_patches)
    : Base{n_comps, box, n_patches},
      storage_({box.im(0), box.im(1), box.im(2), n_comps, n_patches})
  {}

  auto gt() { return Base::storage().view(); }

private:
  Storage storage_;

  KG_INLINE Storage& storageImpl() { return storage_; }
  KG_INLINE const Storage& storageImpl() const { return storage_; }

  friend class MfieldsCRTP<CudaMfields>;
};

// ======================================================================
// MfieldsCuda

struct MfieldsCuda : MfieldsBase
{
  using real_t = fields_cuda_t::real_t;
  using Real = real_t;

  MfieldsCuda(const Grid_t& grid, int n_fields, Int3 ibn);
  MfieldsCuda(const MfieldsCuda&) = delete;
  MfieldsCuda(MfieldsCuda&&) = default;
  ~MfieldsCuda();

  cuda_mfields* cmflds() { return cmflds_; }
  const cuda_mfields* cmflds() const { return cmflds_; }

  int n_comps() const;
  int n_patches() const;
  const Grid_t& grid() const { return *grid_; }

  void reset(const Grid_t& new_grid) override;

  int index(int m, int i, int j, int k, int p) const;

  gt::gtensor_span_device<real_t, 5> gt();

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  cuda_mfields* cmflds_;
  const Grid_t* grid_;
};

MfieldsSingle hostMirror(MfieldsCuda& mflds);
MfieldsSingle hostMirror(const MfieldsCuda& mflds);
void copy(const MfieldsCuda& mflds, MfieldsSingle& hmflds);
void copy(const MfieldsSingle& hmflds, MfieldsCuda& mflds);

// ======================================================================
// MfieldsStateCuda

struct MfieldsStateCuda : MfieldsStateBase
{
  using real_t = MfieldsCuda::real_t;
  using space = gt::space::device;

  MfieldsStateCuda(const Grid_t& grid)
    : MfieldsStateBase{grid, NR_FIELDS, grid.ibn},
      mflds_{grid, NR_FIELDS, grid.ibn}
  {}

  /* void reset(const Grid_t& new_grid) override */
  /* { */
  /*   MfieldsStateBase::reset(new_grid); */
  /*   mflds_.reset(new_grid); */
  /* } */

  cuda_mfields* cmflds() { return mflds_.cmflds(); }

  int n_patches() const { return mflds_.n_patches(); };
  const Grid_t& grid() const { return *grid_; }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  MfieldsCuda& mflds() { return mflds_; }
  const MfieldsCuda& mflds() const { return mflds_; }

  auto gt() { return mflds_.gt(); }

private:
  MfieldsCuda mflds_;
};

template <>
struct Mfields_traits<MfieldsCuda>
{
  static constexpr const char* name = "cuda";
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

#endif
