
#include "psc_output_fields_item_private.h"

// ======================================================================

#define define_dxdydz(dx, dy, dz)					\
  int dx _mrc_unused = (ppsc->domain.gdims[0] == 1) ? 0 : 1;		\
  int dy _mrc_unused = (ppsc->domain.gdims[1] == 1) ? 0 : 1;		\
  int dz _mrc_unused = (ppsc->domain.gdims[2] == 1) ? 0 : 1

// ======================================================================

static void
calc_dive_nc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	     struct psc_mparticles *mprts, struct psc_mfields *mres_base)
{
  define_dxdydz(dx, dy, dz);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, EX, EX + 3);
  struct psc_mfields *mres = psc_mfields_get_as(mres_base, FIELDS_TYPE, 0, 0);
  for (int p = 0; p < mres_base->nr_patches; p++) {
    fields_t flds = fields_t_mflds(mflds, p);
    fields_t res = fields_t_mflds(mres, p);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      _F3(res, 0, ix,iy,iz) = 
	((_F3(flds, EX, ix,iy,iz) - _F3(flds, EX, ix-dx,iy,iz)) / ppsc->patch[p].dx[0] +
	 (_F3(flds, EY, ix,iy,iz) - _F3(flds, EY, ix,iy-dy,iz)) / ppsc->patch[p].dx[1] +
	 (_F3(flds, EZ, ix,iy,iz) - _F3(flds, EZ, ix,iy,iz-dz)) / ppsc->patch[p].dx[2]);
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
  psc_mfields_put_as(mres, mres_base, 0, 1);
}

struct psc_output_fields_item_ops psc_output_fields_item_dive_ops = {
  .name      = "dive_" FIELDS_TYPE,
  .nr_comp   = 1,
  .fld_names = { "dive" },
  .run_all   = calc_dive_nc,
};

