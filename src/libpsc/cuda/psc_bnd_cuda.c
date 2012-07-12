
#include "psc_bnd_cuda.h"
#include "psc_cuda.h"
#include "psc_bnd_cuda_fields.h"
#include "particles_cuda.h"

// ======================================================================
// psc_bnd: subclass "cuda"

struct psc_bnd_ops psc_bnd_cuda_ops = {
  .name                    = "cuda",
  .create_ddc              = psc_bnd_fld_cuda_create,
  .add_ghosts              = psc_bnd_fld_cuda_add_ghosts,
  .fill_ghosts             = psc_bnd_fld_cuda_fill_ghosts,
};

