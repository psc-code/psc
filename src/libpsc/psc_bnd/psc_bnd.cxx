
#include "psc_bnd_private.h"

#include <mrc_io.h>
#include <mrc_ddc.h>
#include <mrc_profile.h>

// ======================================================================
// psc_bnd

// ----------------------------------------------------------------------
// psc_set_psc

void
psc_bnd_set_psc(struct psc_bnd *bnd, struct psc *psc)
{
  bnd->psc = psc;
}

// ----------------------------------------------------------------------
// psc_bnd_write

static void
_psc_bnd_write(struct psc_bnd *bnd, struct mrc_io *io)
{
  mrc_io_write_ref(io, bnd, "psc", bnd->psc);
}

// ----------------------------------------------------------------------
// psc_bnd_read

static void
_psc_bnd_read(struct psc_bnd *bnd, struct mrc_io *io)
{
  bnd->psc = mrc_io_read_ref(io, bnd, "psc", psc);

  psc_bnd_setup(bnd);
}

// ======================================================================
// psc_bnd_init

extern struct psc_bnd_ops psc_bnd_c_ops;
extern struct psc_bnd_ops psc_bnd_single_ops;
extern struct psc_bnd_ops psc_bnd_cuda_ops;
extern struct psc_bnd_ops psc_bnd_vpic_ops;

static void
psc_bnd_init()
{
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_single_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_cuda_ops);
#endif
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_vpic_ops);
#endif
}

// ======================================================================
// psc_bnd class

struct mrc_class_psc_bnd_ : mrc_class_psc_bnd
{
  mrc_class_psc_bnd_()
  {
    name             = "psc_bnd";
    size             = sizeof(struct psc_bnd);
    init             = psc_bnd_init;
    write            = _psc_bnd_write;
    read             = _psc_bnd_read;
  }
} mrc_class_psc_bnd;

