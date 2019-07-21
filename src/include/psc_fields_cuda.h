
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

struct MfieldsCuda : MfieldsBase
{
  using real_t = fields_cuda_t::real_t;
  using fields_host_t = kg::SArray<real_t>;

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
  void write_as_mrc_fld(mrc_io *io, const std::string& name, const std::vector<std::string>& comp_names) override;

  void zero();
  void axpy_comp_yz(int ym, float a, MfieldsCuda& x, int xm);

  int index(int m, int i, int j, int k, int p) const;
  Patch operator[](int p) { return { *this, p }; }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }
  
  cuda_mfields* cmflds_;
  const Grid_t* grid_;
};

MfieldsCuda::fields_host_t get_host_fields(const MfieldsCuda& mflds);
void copy_to_device(int p, const MfieldsCuda::fields_host_t& h_flds, MfieldsCuda& mflds, int mb, int me);
void copy_from_device(int p, MfieldsCuda::fields_host_t& h_flds, const MfieldsCuda& mflds, int mb, int me);

// ======================================================================
// MfieldsStateCuda

struct MfieldsStateCuda : MfieldsStateBase
{
  using real_t = MfieldsCuda::real_t;
  using fields_host_t = MfieldsCuda::fields_host_t;
  
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

inline MfieldsStateCuda::fields_host_t get_host_fields(const MfieldsStateCuda& mflds)
{
  return ::get_host_fields(mflds.mflds());
}

inline void copy_to_device(int p, const MfieldsStateCuda::fields_host_t& h_flds, MfieldsStateCuda& mflds, int mb, int me)
{
  copy_to_device(p, h_flds, mflds.mflds(), mb, me);
}

inline void copy_from_device(int p, MfieldsStateCuda::fields_host_t& h_flds, const MfieldsStateCuda& mflds, int mb, int me)
{
  copy_from_device(p, h_flds, mflds.mflds(), mb, me);
}

// ----------------------------------------------------------------------

#endif
