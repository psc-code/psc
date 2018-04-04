
#include "psc_marder_private.h"
#include "psc_bnd.h"
#include "psc_output_fields_item.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

#include "marder.hxx"
#include "fields_item.hxx"
#include "marder_cuda_impl.hxx"

#include <stdlib.h>

// ======================================================================
// psc_marder: subclass "cuda"

struct psc_marder_ops_cuda : psc_marder_ops {
  using Wrapper = MarderWrapper<MarderCuda>;
  psc_marder_ops_cuda() {
    name                  = "cuda";
    size                  = Wrapper::size;
    setup                 = Wrapper::setup;
    destroy               = Wrapper::destroy;
  }
} psc_marder_cuda_ops;

