#include "psc_output_fields_item_private.h"

#include "fields.hxx"
#include "psc_fields_as_c.h"

using Fields = Fields3d<mfields_t::fields_t>;

// FIXME, we're assuming that the result fields are "c" type

// ======================================================================

#define define_dxdydz(dx, dy, dz)					\
  int dx _mrc_unused = (ppsc->domain.gdims[0] == 1) ? 0 : 1;		\
  int dy _mrc_unused = (ppsc->domain.gdims[1] == 1) ? 0 : 1;		\
  int dz _mrc_unused = (ppsc->domain.gdims[2] == 1) ? 0 : 1

// ======================================================================

#define JX_NC(ix,iy,iz) (.5f * (F(JXI,ix,iy,iz) + F(JXI,ix-dx,iy,iz)))
#define JY_NC(ix,iy,iz) (.5f * (F(JYI,ix,iy,iz) + F(JYI,ix,iy-dy,iz)))
#define JZ_NC(ix,iy,iz) (.5f * (F(JZI,ix,iy,iz) + F(JZI,ix,iy,iz-dz)))

static void
calc_j_nc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(JXI, JXI + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = JX_NC(ix,iy,iz);
      R(1, ix,iy,iz) = JY_NC(ix,iy,iz);
      R(2, ix,iy,iz) = JZ_NC(ix,iy,iz);
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_j_nc : psc_output_fields_item_ops {
  psc_output_fields_item_ops_j_nc() {
    name      = "j_nc";
    nr_comp   = 3;
    fld_names = { "jx_nc", "jy_nc", "jz_nc" };
    run_all   = calc_j_nc;
  }
} psc_output_fields_item_j_nc_ops;

// ======================================================================

#define JX_CC(ix,iy,iz) (.25f * (F(JXI, ix,iy,iz   ) + F(JXI, ix,iy+dy,iz   ) + \
				 F(JXI, ix,iy,iz+dz) + F(JXI, ix,iy+dy,iz+dz)))
#define JY_CC(ix,iy,iz) (.25f * (F(JYI, ix,iy,iz   ) + F(JYI, ix+dx,iy,iz   ) + \
				 F(JYI, ix,iy,iz+dz) + F(JYI, ix+dx,iy,iz+dz)))
#define JZ_CC(ix,iy,iz) (.25f * (F(JZI, ix,iy,iz   ) + F(JZI, ix+dx,iy   ,iz) + \
				 F(JZI, ix,iy+dy,iz) + F(JZI, ix+dx,iy+dy,iz)))

static void
calc_j(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
       struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(JXI, JXI + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = JX_CC(ix,iy,iz);
      R(1, ix,iy,iz) = JY_CC(ix,iy,iz);
      R(2, ix,iy,iz) = JZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_j : psc_output_fields_item_ops {
  psc_output_fields_item_ops_j() {
    name      = "j";
    nr_comp   = 3;
    fld_names[0] = "jx";
    fld_names[1] = "jy";
    fld_names[2] = "jz";
    run_all   = calc_j;
  }
} psc_output_fields_item_j_ops;

// ======================================================================

static void
calc_j_ec(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(JXI, JXI + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = F(JXI, ix,iy,iz);
      R(1, ix,iy,iz) = F(JYI, ix,iy,iz);
      R(2, ix,iy,iz) = F(JZI, ix,iy,iz);
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_j_ec : psc_output_fields_item_ops {
  psc_output_fields_item_ops_j_ec() {
    name      = "j_ec";
    nr_comp   = 3;
    fld_names[0] = "jx_ec";
    fld_names[1] = "jy_ec";
    fld_names[2] = "jz_ec";
    run_all   = calc_j_ec;
  }
} psc_output_fields_item_j_ec_ops;

// ======================================================================

#define EX_NC(ix,iy,iz) (.5f * (F( EX,ix,iy,iz) + F( EX,ix-dx,iy,iz)))
#define EY_NC(ix,iy,iz) (.5f * (F( EY,ix,iy,iz) + F( EY,ix,iy-dy,iz)))
#define EZ_NC(ix,iy,iz) (.5f * (F( EZ,ix,iy,iz) + F( EZ,ix,iy,iz-dz)))

static void
calc_E_nc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(EX, EX + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = EX_NC(ix,iy,iz);
      R(1, ix,iy,iz) = EY_NC(ix,iy,iz);
      R(2, ix,iy,iz) = EZ_NC(ix,iy,iz);
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_e_nc : psc_output_fields_item_ops {
  psc_output_fields_item_ops_e_nc() {
    name      = "e_nc";
    nr_comp   = 3;
    fld_names[0] = "ex_nc";
    fld_names[1] = "ey_nc";
    fld_names[2] = "ez_nc";
    run_all   = calc_E_nc;
  }
} psc_output_fields_item_e_nc_ops;

// ======================================================================

#define EX_CC(ix,iy,iz) (.25f * (F( EX,ix,iy,iz   ) + F( EX,ix,iy+dy,iz   ) + \
				 F( EX,ix,iy,iz+dz) + F( EX,ix,iy+dy,iz+dz)))
#define EY_CC(ix,iy,iz) (.25f * (F( EY,ix,iy,iz   ) + F( EY,ix+dx,iy,iz   ) + \
				 F( EY,ix,iy,iz+dz) + F( EY,ix+dx,iy,iz+dz)))
#define EZ_CC(ix,iy,iz) (.25f * (F( EZ,ix,iy   ,iz) + F( EZ,ix+dx,iy   ,iz) + \
				 F( EZ,ix,iy+dy,iz) + F( EZ,ix+dx,iy+dy,iz)))

static void
calc_E_cc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(EX, EX + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = EX_CC(ix,iy,iz);
      R(1, ix,iy,iz) = EY_CC(ix,iy,iz);
      R(2, ix,iy,iz) = EZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_e : psc_output_fields_item_ops {
  psc_output_fields_item_ops_e() {
    name      = "e";
    nr_comp   = 3;
    fld_names[0] = "ex";
    fld_names[1] = "ey";
    fld_names[2] = "ez";
    run_all   = calc_E_cc;
  }
} psc_output_fields_item_e_ops;

// ======================================================================

static void
calc_E_ec(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(EX, EX + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = F(EX, ix,iy,iz);
      R(1, ix,iy,iz) = F(EY, ix,iy,iz);
      R(2, ix,iy,iz) = F(EZ, ix,iy,iz);
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_e_ec : psc_output_fields_item_ops {
  psc_output_fields_item_ops_e_ec() {
    name      = "e_ec";
    nr_comp   = 3;
    fld_names[0] = "ex_ec";
    fld_names[1] = "ey_ec";
    fld_names[2] = "ez_ec";
    run_all   = calc_E_ec;
  }
} psc_output_fields_item_e_ec_ops;

// ======================================================================

#define HX_NC(ix,iy,iz) (.25f*(F(HX,ix,iy,iz   ) + F(HX,ix,iy-dy,iz   ) + \
			       F(HX,ix,iy,iz-dz) + F(HX,ix,iy-dy,iz-dz)))
#define HY_NC(ix,iy,iz) (.25f*(F(HY,ix,iy,iz   ) + F(HY,ix-dx,iy,iz   ) + \
			       F(HY,ix,iy,iz-dz) + F(HY,ix-dx,iy,iz-dz)))
#define HZ_NC(ix,iy,iz) (.25f*(F(HZ,ix,iy   ,iz) + F(HZ,ix-dx,iy   ,iz) + \
			       F(HZ,ix,iy-dy,iz) + F(HZ,ix-dx,iy-dy,iz)))

static void
calc_H_nc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(HX, HX + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = HX_NC(ix,iy,iz);
      R(1, ix,iy,iz) = HY_NC(ix,iy,iz);
      R(2, ix,iy,iz) = HZ_NC(ix,iy,iz);
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_h_nc : psc_output_fields_item_ops {
  psc_output_fields_item_ops_h_nc() {
    name      = "h_nc";
    nr_comp   = 3;
    fld_names[0] = "hx_nc";
    fld_names[1] = "hy_nc";
    fld_names[2] = "hz_nc";
    run_all   = calc_H_nc;
  }
} psc_output_fields_item_h_nc_ops;

// ======================================================================

#define HX_CC(ix,iy,iz) (.5f*(F(HX,ix,iy,iz) + F(HX,ix+dx,iy,iz)))
#define HY_CC(ix,iy,iz) (.5f*(F(HY,ix,iy,iz) + F(HY,ix,iy+dy,iz)))
#define HZ_CC(ix,iy,iz) (.5f*(F(HZ,ix,iy,iz) + F(HZ,ix,iy,iz+dz)))

static void
calc_H_cc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(HX, HX + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = HX_CC(ix,iy,iz);
      R(1, ix,iy,iz) = HY_CC(ix,iy,iz);
      R(2, ix,iy,iz) = HZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_h : psc_output_fields_item_ops {
  psc_output_fields_item_ops_h() {
    name      = "h";
    nr_comp   = 3;
    fld_names[0] = "hx";
    fld_names[1] = "hy";
    fld_names[2] = "hz";
    run_all   = calc_H_cc;
  }
} psc_output_fields_item_h_ops;

// ======================================================================

static void
calc_H_fc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(HX, HX + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = F(HX, ix,iy,iz);
      R(1, ix,iy,iz) = F(HY, ix,iy,iz);
      R(2, ix,iy,iz) = F(HZ, ix,iy,iz);
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_h_fc : psc_output_fields_item_ops {
  psc_output_fields_item_ops_h_fc() {
    name      = "h_fc";
    nr_comp   = 3;
    fld_names[0] = "hx_fc";
    fld_names[1] = "hy_fc";
    fld_names[2] = "hz_fc";
    run_all   = calc_H_fc;
  }
} psc_output_fields_item_h_fc_ops;

// ======================================================================

static void
calc_jdote_cc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	      struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(JXI, EX + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = JX_CC(ix,iy,iz) * EX_CC(ix,iy,iz);
      R(1, ix,iy,iz) = JY_CC(ix,iy,iz) * EY_CC(ix,iy,iz);
      R(2, ix,iy,iz) = JZ_CC(ix,iy,iz) * EZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_jdote : psc_output_fields_item_ops {
  psc_output_fields_item_ops_jdote() {
    name      = "jdote";
    nr_comp   = 3;
    fld_names[0] = "jxex";
    fld_names[1] = "jyey";
    fld_names[2] = "jzez";
    run_all   = calc_jdote_cc;
  }
} psc_output_fields_item_jdote_ops;

// ======================================================================

static void
calc_poyn_cc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	     struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(EX, HX + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = (EY_CC(ix,iy,iz) * HZ_CC(ix,iy,iz) - 
			       EZ_CC(ix,iy,iz) * HY_CC(ix,iy,iz));
      R(1, ix,iy,iz) = (EZ_CC(ix,iy,iz) * HX_CC(ix,iy,iz) -
			       EX_CC(ix,iy,iz) * HZ_CC(ix,iy,iz));
      R(2, ix,iy,iz) = (EX_CC(ix,iy,iz) * HY_CC(ix,iy,iz) -
			       EY_CC(ix,iy,iz) * HX_CC(ix,iy,iz));
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_poyn : psc_output_fields_item_ops {
  psc_output_fields_item_ops_poyn() {
    name      = "poyn";
    nr_comp   = 3;
    fld_names[0] = "poynx";
    fld_names[1] = "poyny";
    fld_names[2] = "poynz";
    run_all   = calc_poyn_cc;
  }
} psc_output_fields_item_poyn_ops;

// ======================================================================

static void
calc_E2_cc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	   struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(EX, EX + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = sqr(EX_CC(ix,iy,iz));
      R(1, ix,iy,iz) = sqr(EY_CC(ix,iy,iz));
      R(2, ix,iy,iz) = sqr(EZ_CC(ix,iy,iz));
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_e2 : psc_output_fields_item_ops {
  psc_output_fields_item_ops_e2() {
    name      = "e2";
    nr_comp   = 3;
    fld_names[0] = "ex2";
    fld_names[1] = "ey2";
    fld_names[2] = "ez2";
    run_all   = calc_E2_cc;
  }
} psc_output_fields_item_e2_ops;

// ======================================================================

static void
calc_H2_cc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	   struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(HX, HX + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = sqr(HX_CC(ix,iy,iz));
      R(1, ix,iy,iz) = sqr(HY_CC(ix,iy,iz));
      R(2, ix,iy,iz) = sqr(HZ_CC(ix,iy,iz));
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_h2 : psc_output_fields_item_ops {
  psc_output_fields_item_ops_h2() {
    name      = "h2";
    nr_comp   = 3;
    fld_names[0] = "hx2";
    fld_names[1] = "hy2";
    fld_names[2] = "hz2";
    run_all   = calc_H2_cc;
  }
} psc_output_fields_item_h2_ops;

// ======================================================================

static void
calc_divb(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(HX, HX + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = 
	((F(HX, ix+dx,iy,iz) - F(HX, ix,iy,iz)) / ppsc->grid().dx[0] +
	 (F(HY, ix,iy+dy,iz) - F(HY, ix,iy,iz)) / ppsc->grid().dx[1] +
	 (F(HZ, ix,iy,iz+dz) - F(HZ, ix,iy,iz)) / ppsc->grid().dx[2]);
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_divb : psc_output_fields_item_ops {
  psc_output_fields_item_ops_divb() {
    name      = "divb";
    nr_comp   = 1;
    fld_names[0] = "divb";
    run_all   = calc_divb;
  }
} psc_output_fields_item_divb_ops;

// ======================================================================

static void
calc_divj_nc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	     struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  mfields_t mf = mflds_base->get_as<mfields_t>(JXI, JXI + 3);
  mfields_t mf_res(mres);
  for (int p = 0; p < mf_res->n_patches(); p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = 
	((F(JXI, ix,iy,iz) - F(JXI, ix-dx,iy,iz)) / ppsc->grid().dx[0] +
	 (F(JYI, ix,iy,iz) - F(JYI, ix,iy-dy,iz)) / ppsc->grid().dx[1] +
	 (F(JZI, ix,iy,iz) - F(JZI, ix,iy,iz-dz)) / ppsc->grid().dx[2]);
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}

struct psc_output_fields_item_ops_divj : psc_output_fields_item_ops {
  psc_output_fields_item_ops_divj() {
    name      = "divj";
    nr_comp   = 1;
    fld_names[0] = "divj";
    run_all   = calc_divj_nc;
  }
} psc_output_fields_item_divj_ops;

