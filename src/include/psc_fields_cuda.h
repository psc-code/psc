
#ifndef PSC_FIELDS_CUDA_H
#define PSC_FIELDS_CUDA_H

#include <mpi.h>
#include "fields3d.hxx"
#include "fields_traits.hxx"

struct fields_cuda_t
{
  using real_t = float;
  using mfields_t = mfields_base<fields_cuda_t>;
};

struct psc_mfields_cuda
{
  using fields_t = fields_cuda_t;
  
  struct cuda_mfields *cmflds;
};

using mfields_cuda_t = mfields_base<psc_mfields_cuda>;

template<>
struct fields_traits<fields_cuda_t>
{
  static constexpr const char* name = "cuda";
};

// ----------------------------------------------------------------------

#define psc_mfields_cuda(pf) mrc_to_subobj(pf, struct psc_mfields_cuda)

#endif
