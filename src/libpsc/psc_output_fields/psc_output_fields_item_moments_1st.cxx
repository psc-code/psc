
#include "psc.h"

// ======================================================================
// !!! These moments are shifted to (n+.5) * dt, rather than n * dt,
// since the "single/double" particles are shifted that way.

#include "fields_item_moments_1st.hxx"

#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "psc_fields_c.h"

template<typename Moment_t>
using Ops = FieldsItemMomentOps<Moment_t>;
  
#define MAKE_POFI_OPS(MP, MF, TYPE)					\
  Ops<Moment_n_1st<MP, MF>> psc_output_fields_item_n_1st_##TYPE##_ops;	\
  Ops<Moment_v_1st<MP, MF>> psc_output_fields_item_v_1st_##TYPE##_ops;	\
  Ops<Moment_p_1st<MP, MF>> psc_output_fields_item_p_1st_##TYPE##_ops;	\
  Ops<Moment_vv_1st<MP, MF>> psc_output_fields_item_vv_1st_##TYPE##_ops; \
  Ops<Moment_T_1st<MP, MF>> psc_output_fields_item_T_1st_##TYPE##_ops;	\
  Ops<Moment_Tvv_1st<MP, MF>> psc_output_fields_item_Tvv_1st_##TYPE##_ops; \

MAKE_POFI_OPS(MparticlesSingle, MfieldsC, single);
MAKE_POFI_OPS(MparticlesDouble, MfieldsC, double);
