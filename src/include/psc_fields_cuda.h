
#ifndef PSC_FIELDS_CUDA_H
#define PSC_FIELDS_CUDA_H

#include <mpi.h>
#include "fields3d.hxx"
#include "fields_traits.hxx"
#include <kg/SArray.h>

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
      storage_(n_patches * n_comps * box.size())
  {}

private:
  Storage storage_;
  
  KG_INLINE Storage& storageImpl() { return storage_; }
  KG_INLINE const Storage& storageImpl() const { return storage_; }

  friend class MfieldsCRTP<CudaMfields>;
};

using HMFields = CudaMfields<std::vector<float>>;

// ======================================================================
// MfieldsCuda

struct MfieldsCuda : MfieldsBase
{
  using real_t = fields_cuda_t::real_t;

  class Accessor
  {
  public:
    Accessor(MfieldsCuda& mflds, int idx);
    operator real_t() const;
    real_t operator=(real_t val);
    real_t operator+=(real_t val);

  private:
    MfieldsCuda& mflds_;
    int idx_;
  };
  
  class Patch
  {
  public:
    Patch(MfieldsCuda& mflds, int p);
    
    Accessor operator()(int m, int i, int j, int k);
    
  private:
    MfieldsCuda& mflds_;
    int p_;
  };

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
  void zero_comp(int m);

  void zero();
  void copy_comp(int m_to, MfieldsCuda& from, int m_from);
  void axpy_comp(int ym, float a, MfieldsCuda& x, int xm);

  int index(int m, int i, int j, int k, int p) const;
  Patch operator[](int p) { return { *this, p }; }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }
  
  cuda_mfields* cmflds_;
  const Grid_t* grid_;
};

HMFields hostMirror(const MfieldsCuda& mflds);
void copy(const MfieldsCuda& mflds, HMFields& hmflds);
void copy(const HMFields& hmflds, MfieldsCuda& mflds);

// ======================================================================
// MfieldsStateCuda

struct MfieldsStateCuda : MfieldsStateBase
{
  using real_t = MfieldsCuda::real_t;
  
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

  MfieldsCuda::Patch operator[](int p) { return mflds_[p]; }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  MfieldsCuda& mflds() { return mflds_; }
  const MfieldsCuda& mflds() const { return mflds_; }
  
private:
  MfieldsCuda mflds_;
};

template<>
struct Mfields_traits<MfieldsCuda>
{
  static constexpr const char* name = "cuda";
};

inline HMFields hostMirror(const MfieldsStateCuda& mflds)
{
  return hostMirror(mflds.mflds());
}

inline void copy(const MfieldsStateCuda& mflds, HMFields& hmflds)
{
  copy(mflds.mflds(), hmflds);
}

inline void copy(const HMFields& hmflds, MfieldsStateCuda& mflds)
{
  copy(hmflds, mflds.mflds());
}

// ----------------------------------------------------------------------

#endif
