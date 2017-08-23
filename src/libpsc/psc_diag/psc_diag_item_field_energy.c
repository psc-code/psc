
#include "psc_diag_item_private.h"
#include "psc_fields_c.h"

// ----------------------------------------------------------------------
// psc_diag_item_field_energy_run

static void
psc_diag_item_field_energy_run(struct psc_diag_item *item,
			       struct psc *psc, double *EH2)
{
  struct psc_mfields *mflds = psc_mfields_get_as(psc->flds, "c", EX, HX + 3);
  psc_foreach_patch(psc, p) {
    double fac = psc->patch[p].dx[0] * psc->patch[p].dx[1] * psc->patch[p].dx[2];
    struct psc_fields *pf = psc_mfields_get_patch(mflds, p);
    // FIXME, this doesn't handle non-periodic b.c. right
    psc_foreach_3d(psc, p, ix, iy, iz, 0, 0) {
      EH2[0] +=	sqr(F3_C(pf, EX, ix,iy,iz)) * fac;
      EH2[1] +=	sqr(F3_C(pf, EY, ix,iy,iz)) * fac;
      EH2[2] +=	sqr(F3_C(pf, EZ, ix,iy,iz)) * fac;
      EH2[3] +=	sqr(F3_C(pf, HX, ix,iy,iz)) * fac;
      EH2[4] +=	sqr(F3_C(pf, HY, ix,iy,iz)) * fac;
      EH2[5] +=	sqr(F3_C(pf, HZ, ix,iy,iz)) * fac;
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, psc->flds, 0, 0);
}

// ======================================================================
// psc_diag_item_field_energy

struct psc_diag_item_ops psc_diag_item_field_energy_ops = {
  .name      = "field_energy",
  .run       = psc_diag_item_field_energy_run,
  .nr_values = 6,
  .title     = { "EX2", "EY2", "EZ2", "HX2", "HY2", "HZ2" },
};

