
#include "cuda_mfields.h"
#include "cuda_bits.h"
#include "cuda_base.cuh"

#include "fields.hxx"

#include <cstdio>
#include <cassert>

// ======================================================================
// cuda_mfields

// ----------------------------------------------------------------------
// cast to DMFields

cuda_mfields::operator DMFields()
{
  return DMFields{box(), n_comps(), n_patches(), storage().data().get()};
}

// ----------------------------------------------------------------------
// operator[]

DFields cuda_mfields::operator[](int p) const
{
  return static_cast<DMFields>(const_cast<cuda_mfields&>(*this))[p];
}

MfieldsSingle hostMirror(cuda_mfields& cmflds)
{
  return MfieldsSingle{cmflds.grid(), cmflds.n_comps(), -cmflds.ib()};
}

MfieldsSingle hostMirror(const cuda_mfields& cmflds)
{
  return MfieldsSingle{cmflds.grid(), cmflds.n_comps(), -cmflds.ib()};
}

void copy(const cuda_mfields& cmflds, MfieldsSingle& hmflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("cmflds to host", 1., 0, 0);
  }
  prof_start(pr);
  thrust::copy(cmflds.storage().data(),
               cmflds.storage().data() + cmflds.storage().size(),
               hmflds.storage().data());
  prof_stop(pr);
}

void copy(const MfieldsSingle& hmflds, cuda_mfields& cmflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("cmflds from host", 1., 0, 0);
  }
  prof_start(pr);
  thrust::copy(hmflds.storage().data(),
               hmflds.storage().data() + hmflds.storage().size(),
               cmflds.storage().data());
  prof_stop(pr);
}
