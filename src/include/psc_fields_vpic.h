
#ifndef PSC_FIELDS_VPIC_H
#define PSC_FIELDS_VPIC_H

#include "fields3d.hxx"

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
    return Simulation_mflds_synchronize_tang_e_norm_b(sim, vmflds_fields);
  }
  
  void compute_div_b_err()
  {
    Simulation_mflds_compute_div_b_err(sim, vmflds_fields);
  }
  
  double compute_rms_div_b_err()
  {
    return Simulation_mflds_compute_rms_div_b_err(sim, vmflds_fields);
  }

  void clean_div_b()
  {
    Simulation_mflds_clean_div_b(sim, vmflds_fields);
  }

  void compute_div_e_err()
  {
    Simulation_mflds_compute_div_e_err(sim, vmflds_fields);
  }

  double compute_rms_div_e_err()
  {
    return Simulation_mflds_compute_rms_div_e_err(sim, vmflds_fields);
  }

  void clean_div_e()
  {
    Simulation_mflds_clean_div_e(sim, vmflds_fields);
  }

  void clear_rhof()
  {
    Simulation_mflds_clear_rhof(sim, vmflds_fields);
  }

  void synchronize_rho()
  {
    Simulation_mflds_synchronize_rho(sim, vmflds_fields);
  }
  
  void compute_rhob()
  {
    Simulation_mflds_compute_rhob(sim, vmflds_fields);
  }
  
  void compute_curl_b()
  {
    Simulation_mflds_compute_curl_b(sim, vmflds_fields);
  }

  void accumulate_rho_p(Particles* vmprts);

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  Simulation *sim;
  FieldArray *vmflds_fields;
  HydroArray *vmflds_hydro;
};

using PscMfieldsVpic = PscMfields<MfieldsVpic>;

template<>
struct fields_traits<fields_vpic_t>
{
  static constexpr const char* name = "vpic";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#include "vpic_iface.h"

#endif
