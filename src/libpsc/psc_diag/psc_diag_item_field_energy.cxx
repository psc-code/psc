
#include "psc_diag_item_private.h"
#include "psc_fields_as_c.h"
#include "fields.hxx"

using Fields = Fields3d<mfields_t::fields_t, dim_xyz>;

// ----------------------------------------------------------------------
// psc_diag_item_field_energy_run

static void
psc_diag_item_field_energy_run(struct psc_diag_item *item,
			       struct psc *psc, double *EH2)
{
  const Grid_t& grid = psc->grid();
  mfields_t mf = psc->flds->get_as<mfields_t>(EX, HX + 3);
  psc_foreach_patch(psc, p) {
    double fac = grid.dx[0] * grid.dx[1] * grid.dx[2];
    Fields F(mf[p]);
    // FIXME, this doesn't handle non-periodic b.c. right
    psc_foreach_3d(psc, p, ix, iy, iz, 0, 0) {
      EH2[0] +=	sqr(F(EX, ix,iy,iz)) * fac;
      EH2[1] +=	sqr(F(EY, ix,iy,iz)) * fac;
      EH2[2] +=	sqr(F(EZ, ix,iy,iz)) * fac;
      EH2[3] +=	sqr(F(HX, ix,iy,iz)) * fac;
      EH2[4] +=	sqr(F(HY, ix,iy,iz)) * fac;
      EH2[5] +=	sqr(F(HZ, ix,iy,iz)) * fac;
    } foreach_3d_end;
  }
  mf.put_as(psc->flds, 0, 0);
}

// ======================================================================
// psc_diag_item_field_energy

struct psc_diag_item_ops_fe : psc_diag_item_ops {
  psc_diag_item_ops_fe() {
    name      = "field_energy";
    run       = psc_diag_item_field_energy_run;
    nr_values = 6;
    title[0]  = "EX2";
    title[1]  = "EXY";
    title[2]  = "EXZ";
    title[3]  = "HX2";
    title[4]  = "HXY";
    title[5]  = "HXZ";
  }
} psc_diag_item_field_energy_ops;

