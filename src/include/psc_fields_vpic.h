
#ifndef PSC_FIELDS_VPIC_H
#define PSC_FIELDS_VPIC_H

#include "fields3d.hxx"

#include "../libpsc/vpic/vpic_iface.h" // FIXME path

struct fields_vpic_t : fields3d<float, LayoutAOS>
{
  using mfields_t = mfields_base<fields_vpic_t>;
  using Base = fields3d<float, LayoutAOS>;

  using Base::Base;
};

struct psc_mfields_vpic : psc_mfields_base
{
  using fields_t = fields_vpic_t;
  using Base = psc_mfields_base;

  using Base::Base;
  
  Simulation *sim;
  FieldArray *vmflds_fields;
  HydroArray *vmflds_hydro;
};

using mfields_vpic_t = mfields_base<psc_mfields_vpic>;

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
