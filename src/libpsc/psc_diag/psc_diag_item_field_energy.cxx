
#include "psc_diag_item_private.h"
#include "psc_fields_as_c.h"
#include "fields.hxx"

using Fields = Fields3d<MfieldsC::fields_t, dim_xyz>;

// ----------------------------------------------------------------------
// psc_diag_item_field_energy_run

static void
psc_diag_item_field_energy_run(struct psc_diag_item *item,
			       MparticlesBase& mprts,
			       MfieldsStateBase& mflds_base, double *EH2)
{
  const Grid_t& grid = mprts.grid();
  auto& mf = mflds_base.get_as<MfieldsStateDouble>(EX, HX + 3);
  for (int p = 0; p < grid.n_patches(); p++) {
    double fac = grid.domain.dx[0] * grid.domain.dx[1] * grid.domain.dx[2];
    Fields F(mf[p]);
    // FIXME, this doesn't handle non-periodic b.c. right
    grid.Foreach_3d(0, 0, [&](int ix, int iy, int iz) {
	EH2[0] +=	sqr(F(EX, ix,iy,iz)) * fac;
	EH2[1] +=	sqr(F(EY, ix,iy,iz)) * fac;
	EH2[2] +=	sqr(F(EZ, ix,iy,iz)) * fac;
	EH2[3] +=	sqr(F(HX, ix,iy,iz)) * fac;
	EH2[4] +=	sqr(F(HY, ix,iy,iz)) * fac;
	EH2[5] +=	sqr(F(HZ, ix,iy,iz)) * fac;
      });
  }
  mflds_base.put_as(mf, 0, 0);
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

