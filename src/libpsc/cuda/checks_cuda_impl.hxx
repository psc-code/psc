
#pragma once

#include "fields_item_moments_1st_cuda.hxx"

template <typename MP, typename D>
using ChecksCuda =
  ChecksCommon<MP, MfieldsCuda::Storage, Moment_rho_1st_nc_cuda<D>>;
