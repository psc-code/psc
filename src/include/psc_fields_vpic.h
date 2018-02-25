
#ifndef PSC_FIELDS_VPIC_H
#define PSC_FIELDS_VPIC_H

#include "fields3d.hxx"

#include "../libpsc/vpic/vpic_iface.h" // FIXME path

struct fields_vpic_t : fields3d<float, LayoutAOS>
{
  using mfields_t = PscMfields<fields_vpic_t>;
  using Base = fields3d<float, LayoutAOS>;

  using Base::Base;
};

struct psc_mfields_vpic : psc_mfields_base
{
  using fields_t = fields_vpic_t;
  using Base = psc_mfields_base;

  using Base::Base;

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

  Simulation *sim;
  FieldArray *vmflds_fields;
  HydroArray *vmflds_hydro;
};

using mfields_vpic_t = PscMfields<psc_mfields_vpic>;

template<>
inline fields_vpic_t mfields_vpic_t::operator[](int p)
{
  fields_vpic_t psc_mfields_vpic_get_field_t(struct psc_mfields *mflds, int p);
  return psc_mfields_vpic_get_field_t(mflds_, p);
}

template<>
struct fields_traits<fields_vpic_t>
{
  static constexpr const char* name = "vpic";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#include "vpic_iface.h"

BEGIN_C_DECLS

// ----------------------------------------------------------------------

#define psc_mfields_vpic(mflds) ({					\
      mrc_to_subobj(mflds, struct psc_mfields_vpic);			\
    })


END_C_DECLS

#endif
