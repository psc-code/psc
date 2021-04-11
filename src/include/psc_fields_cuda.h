
#ifndef PSC_FIELDS_CUDA_H
#define PSC_FIELDS_CUDA_H

#include <mpi.h>
#include "fields3d.hxx"
#include "fields_traits.hxx"
#include "cuda_mfields.h"

#include "psc_fields_single.h"

#include "mrc_json.h"

// ======================================================================
// MfieldsCuda

struct MfieldsCuda : MfieldsBase
{
  using real_t = float;
  using Real = real_t;

  MfieldsCuda(const Grid_t& grid, int n_fields, Int3 ibn)
    : MfieldsBase{grid, n_fields, ibn},
      grid_{&grid},
      cmflds_{new cuda_mfields(grid, n_fields, ibn)}
  {}
  MfieldsCuda(const MfieldsCuda&) = delete;
  MfieldsCuda(MfieldsCuda&&) = default;
  ~MfieldsCuda() { delete cmflds_; }

  cuda_mfields* cmflds() { return cmflds_; }
  const cuda_mfields* cmflds() const { return cmflds_; }

  int n_comps() const { return cmflds_->n_comps(); }
  int n_patches() const { return cmflds_->n_patches(); };
  const Grid_t& grid() const { return *grid_; }

  void reset(const Grid_t& new_grid) override
  {
    MfieldsBase::reset(new_grid);
    Int3 ibn = -cmflds()->ib();
    int n_comps = cmflds()->n_comps();
    delete cmflds_;
    cmflds_ = new cuda_mfields(new_grid, n_comps, ibn);
    grid_ = &new_grid;
  }

  gt::gtensor_span_device<real_t, 5> gt()
  {
    return gt::adapt_device(cmflds_->storage().data().get(),
                            cmflds_->storage().shape());
  }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  cuda_mfields* cmflds_;
  const Grid_t* grid_;
};

inline MfieldsSingle hostMirror(MfieldsCuda& mflds)
{
  return hostMirror(*mflds.cmflds());
}

inline MfieldsSingle hostMirror(const MfieldsCuda& mflds)
{
  return hostMirror(*mflds.cmflds());
}

inline void copy(const MfieldsCuda& mflds, MfieldsSingle& hmflds)
{
  copy(*mflds.cmflds(), hmflds);
}

inline void copy(const MfieldsSingle& hmflds, MfieldsCuda& mflds)
{
  copy(hmflds, *mflds.cmflds());
}

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
