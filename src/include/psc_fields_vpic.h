
#ifndef PSC_FIELDS_VPIC_H
#define PSC_FIELDS_VPIC_H

#include "fields3d.hxx"
#include "fields_traits.hxx"

#include "../libpsc/vpic/vpic_iface.h" // FIXME path

struct fields_vpic_t : fields3d<float, LayoutAOS>
{
  using Base = fields3d<float, LayoutAOS>;

  using Base::Base;
};

struct MfieldsVpic : MfieldsBase
{
  using fields_t = fields_vpic_t;
  using Base = MfieldsBase;

  using Base::Base;

  MfieldsVpic(const Grid_t& grid, int n_fields, Int3 ibn);
  ~MfieldsVpic();

  inline fields_vpic_t operator[](int p)
  {
    // FIXME hacky...
    if (n_comps() == VPIC_MFIELDS_N_COMP) {
      int ib[3], im[3];
      float* data = Simulation_mflds_getData(sim, vmflds_fields, ib, im);
      return fields_vpic_t(ib, im, VPIC_MFIELDS_N_COMP, data);
    } else if (n_comps() == VPIC_HYDRO_N_COMP) {
      int ib[3], im[3];
      float* data = Simulation_hydro_getData(sim, vmflds_hydro, ib, im);
      return fields_vpic_t(ib, im, VPIC_HYDRO_N_COMP, data);
    } else {
      assert(0);
    }
  }

  double synchronize_tang_e_norm_b()
  {
    double err;
    TIC err = vmflds_fields->synchronize_tang_e_norm_b(); TOC(synchronize_tang_e_norm_b, 1);
    return err;
  }
  
  void compute_div_b_err()
  {
    TIC vmflds_fields->compute_div_b_err(); TOC(compute_div_b_err, 1);
  }
  
  double compute_rms_div_b_err()
  {
    double err;
    TIC err = vmflds_fields->compute_rms_div_b_err(); TOC(compute_rms_div_b_err, 1);
    return err;
  }

  void clean_div_b()
  {
    TIC vmflds_fields->clean_div_b(); TOC(clean_div_b, 1);
  }

  void compute_div_e_err()
  {
    TIC vmflds_fields->compute_div_e_err(); TOC(compute_div_e_err, 1);
  }

  double compute_rms_div_e_err()
  {
    double err;
    TIC err = vmflds_fields->compute_rms_div_e_err(); TOC(compute_rms_div_e_err, 1);
    return err;
  }

  void clean_div_e()
  {
    TIC vmflds_fields->clean_div_e(); TOC(clean_div_e, 1);
  }

  void clear_rhof()
  {
    TIC vmflds_fields->clear_rhof(); TOC(clear_rhof, 1);
  }

  void synchronize_rho()
  {
    TIC vmflds_fields->synchronize_rho(); TOC(synchronize_rho, 1);
  }
  
  void compute_rhob()
  {
    TIC vmflds_fields->compute_rhob(); TOC(compute_rhob, 1);
  }
  
  void compute_curl_b()
  {
    TIC vmflds_fields->compute_curl_b(); TOC(compute_curl_b, 1);
  }

  void accumulate_rho_p(Particles* vmprts);

  void zero_comp(int m) override { assert(0); }
  void set_comp(int m, double val) override { assert(0); }
  void scale_comp(int m, double val) override { assert(0); }
  void axpy_comp(int m_y, double alpha, MfieldsBase& x_base, int m_x) override { assert(0); }
  double max_comp(int m) override { assert(0); }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  Simulation *sim;
  FieldArray *vmflds_fields;
  HydroArray *vmflds_hydro;
};

template<>
struct Mfields_traits<MfieldsVpic>
{
  static constexpr const char* name = "vpic";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#include "../libpsc/vpic/vpic_iface.h"

#endif
