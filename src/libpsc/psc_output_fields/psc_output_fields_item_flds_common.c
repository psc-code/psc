
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
  for (int p = 0; p < mres_base->nr_patches; p++) {
    struct psc_fields *flds_base = psc_mfields_get_patch(mflds_base, p);
    struct psc_fields *res_base = psc_mfields_get_patch(mres_base, p);
    struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, EX, EX + 3);
    struct psc_fields *res = psc_fields_get_as(res_base, FIELDS_TYPE, 0, 0);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      F3(res, 0, ix,iy,iz) = 
	((F3(flds, EX, ix,iy,iz) - F3(flds, EX, ix-dx,iy,iz)) / ppsc->patch[p].dx[0] +
	 (F3(flds, EY, ix,iy,iz) - F3(flds, EY, ix,iy-dy,iz)) / ppsc->patch[p].dx[1] +
	 (F3(flds, EZ, ix,iy,iz) - F3(flds, EZ, ix,iy,iz-dz)) / ppsc->patch[p].dx[2]);
    } foreach_3d_end;
    psc_fields_put_as(flds, flds_base, 0, 0);
    psc_fields_put_as(res, res_base, 0, 1);
  }
}

struct psc_output_fields_item_ops psc_output_fields_item_dive_ops = {
  .name      = "dive_" FIELDS_TYPE,
  .nr_comp   = 1,
  .fld_names = { "dive" },
  .run_all   = calc_dive_nc,
};

