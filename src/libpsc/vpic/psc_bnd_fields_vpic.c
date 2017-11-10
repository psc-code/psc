
#include "psc_bnd_fields_private.h"

#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_bnd_fields_vpic_create_ddc

static void
psc_bnd_fields_vpic_create_ddc(struct psc_bnd *bnd)
{
}

// ----------------------------------------------------------------------
// psc_bnd_fields: subclass "vpic"

struct psc_bnd_fields_ops psc_bnd_fields_vpic_ops = {
  .name                  = "vpic",
  //  .create_ddc            = psc_bnd_vpic_create_ddc,
};

