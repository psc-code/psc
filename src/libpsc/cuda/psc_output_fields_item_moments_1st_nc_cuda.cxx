
#include "psc_output_fields_item_private.h"
#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

#include "fields_item_moments_1st_cuda.hxx"

FieldsItemOps<FieldsItemMoment<Moment_rho_1st_nc_cuda<MparticlesCuda<BS144>, dim_yz>>> psc_output_fields_item_rho_1st_nc_cuda_ops;
FieldsItemOps<FieldsItemMoment<Moment_n_1st_cuda<MparticlesCuda<BS144>, dim_yz>>> psc_output_fields_item_n_1st_cuda_ops;

