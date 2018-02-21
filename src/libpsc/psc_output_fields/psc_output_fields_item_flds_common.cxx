
#include "psc_output_fields_item_private.h"

#include "fields.hxx"

// ======================================================================

#define define_dxdydz(dx, dy, dz)					\
  int dx _mrc_unused = (ppsc->domain.gdims[0] == 1) ? 0 : 1;		\
  int dy _mrc_unused = (ppsc->domain.gdims[1] == 1) ? 0 : 1;		\
  int dz _mrc_unused = (ppsc->domain.gdims[2] == 1) ? 0 : 1

// ======================================================================

template<typename MF>
struct output_flds_ops
{
  using mfields_t = MF;
  using fields_t = typename mfields_t::fields_t;
  using real_t = typename mfields_t::real_t;
  using Fields = Fields3d<fields_t>;

  static void
  calc_dive_nc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	       struct psc_mparticles *mprts, struct psc_mfields *mres_base)
  {
    define_dxdydz(dx, dy, dz);
    mfields_t mf = mflds_base->get_as<mfields_t>(EX, EX + 3);
    mfields_t mf_res = mres_base->get_as<mfields_t>(0, 0);
    for (int p = 0; p < mf_res->n_patches(); p++) {
      Fields F(mf[p]), R(mf_res[p]);
      psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
	R(0, ix,iy,iz) = 
	  ((F(EX, ix,iy,iz) - F(EX, ix-dx,iy,iz)) / ppsc->grid().dx[0] +
	   (F(EY, ix,iy,iz) - F(EY, ix,iy-dy,iz)) / ppsc->grid().dx[1] +
	   (F(EZ, ix,iy,iz) - F(EZ, ix,iy,iz-dz)) / ppsc->grid().dx[2]);
      } foreach_3d_end;
    }
    mf.put_as(mflds_base, 0, 0);
    mf_res.put_as(mres_base, 0, 1);
  }
};

struct psc_output_fields_item_dive_ops : psc_output_fields_item_ops {
  psc_output_fields_item_dive_ops() {
    name      = "dive_" FIELDS_TYPE;
    nr_comp   = 1;
    fld_names[0] = "dive";
    run_all   = output_flds_ops<mfields_t>::calc_dive_nc;
  }
} psc_output_fields_item_dive_ops;

