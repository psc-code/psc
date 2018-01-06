
#include "psc_bnd_private.h"

#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_bnd_vpic_create_ddc

static void
psc_bnd_vpic_create_ddc(struct psc_bnd *bnd)
{
}

// ----------------------------------------------------------------------
// psc_bnd_vpic_fill_ghosts

static void
psc_bnd_vpic_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds,
			 int mb, int me)
{  
}

// ----------------------------------------------------------------------
// psc_bnd_vpic_add_ghosts

static void
psc_bnd_vpic_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds,
			int mb, int me)
{  
}

// ----------------------------------------------------------------------
// psc_bnd: subclass "vpic"

struct psc_bnd_ops psc_bnd_vpic_ops = {
  .name                  = "vpic",
  .create_ddc            = psc_bnd_vpic_create_ddc,
  .fill_ghosts           = psc_bnd_vpic_fill_ghosts,
  .add_ghosts            = psc_bnd_vpic_add_ghosts,
};

