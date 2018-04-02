
#include "psc_bnd_private.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "bnd.hxx"

#include <mrc_ddc_private.h>
#include <mrc_profile.h>

#include "bnd_cuda_impl.hxx"

// ======================================================================
// psc_bnd: subclass "cuda"

struct psc_bnd_ops_cuda : psc_bnd_ops {
  using PscBnd_t = PscBndWrapper<BndCuda>;
  psc_bnd_ops_cuda() {
    name                    = "cuda";
    size                    = PscBnd_t::size;
    setup                   = PscBnd_t::setup;
    destroy                 = PscBnd_t::destroy;
  }
} psc_bnd_cuda_ops;

