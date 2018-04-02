
#ifndef PSC_FIELDS_CUDA_H
#define PSC_FIELDS_CUDA_H

#include <mpi.h>
#include "fields3d.hxx"
#include "fields_traits.hxx"

#include "psc_fields_single.h"

#include "mrc_json.h"

struct fields_cuda_t
{
  using real_t = float;
};

struct MfieldsCuda : MfieldsBase
{
  using fields_t = fields_cuda_t;
  using real_t = fields_cuda_t::real_t;

  class Accessor
  {
  public:
    Accessor(MfieldsCuda& mflds, int idx);
    operator real_t() const;
    real_t operator=(real_t val);

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

  void zero_comp(int m) override;
  void set_comp(int m, double val) override { assert(0); }
  void scale_comp(int m, double val) override { assert(0); }
  void axpy_comp(int m_y, double alpha, MfieldsBase& x, int m_x) override { assert(0); }
  double max_comp(int m) override { assert(0); return 0.; }

  void zero();
  void axpy_comp_yz(int ym, float a, MfieldsCuda& x, int xm);

  fields_single_t get_host_fields();
  void copy_to_device(int p, fields_single_t h_flds, int mb, int me);
  void copy_from_device(int p, fields_single_t h_flds, int mb, int me);

  int index(int m, int i, int j, int k, int p) const;
  Patch operator[](int p) { return { *this, p }; }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }
  
  struct cuda_mfields *cmflds;
};

template<>
struct Mfields_traits<MfieldsCuda>
{
  static constexpr const char* name = "cuda";
};

// ----------------------------------------------------------------------

#endif
