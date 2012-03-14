
#include "psc_diag_item_private.h"

// ----------------------------------------------------------------------
// psc_diag_item_em_energy

static void
psc_diag_item_em_energy_run(struct psc *psc, double *EH2)
{
  mfields_c_t *flds = psc_mfields_get_c(psc->flds, EX, HX + 3);

  double fac = psc->dx[0] * psc->dx[1] * psc->dx[2];
  psc_foreach_patch(psc, p) {
    fields_c_t *pf = psc_mfields_get_patch_c(flds, p);
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

  psc_mfields_put_c(flds, psc->flds, 0, 0);
}

struct psc_diag_item psc_diag_item_em_energy = {
  .run = psc_diag_item_em_energy_run,
  .n_values = 6,
  .names = { "EX2", "EY2", "EZ2", "HX2", "HY2", "HZ2" },
};

