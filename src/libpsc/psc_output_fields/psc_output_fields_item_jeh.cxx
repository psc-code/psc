#include "psc_output_fields_item_private.h"

#include "fields.hxx"
#include "psc_fields_as_c.h"

using Fields = Fields3d<fields_t>;

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
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, JXI, JXI + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = JX_NC(ix,iy,iz);
      R(1, ix,iy,iz) = JY_NC(ix,iy,iz);
      R(2, ix,iy,iz) = JZ_NC(ix,iy,iz);
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_j_nc_ops = {
  .name      = "j_nc",
  .nr_comp   = 3,
  .fld_names = { "jx_nc", "jy_nc", "jz_nc" },
  .run_all   = calc_j_nc,
};

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
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, JXI, JXI + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = JX_CC(ix,iy,iz);
      R(1, ix,iy,iz) = JY_CC(ix,iy,iz);
      R(2, ix,iy,iz) = JZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_j_ops = {
  .name      = "j",
  .nr_comp   = 3,
  .fld_names = { "jx", "jy", "jz" },
  .run_all   = calc_j,
};

// ======================================================================

static void
calc_j_ec(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, JXI, JXI + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = F(JXI, ix,iy,iz);
      R(1, ix,iy,iz) = F(JYI, ix,iy,iz);
      R(2, ix,iy,iz) = F(JZI, ix,iy,iz);
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_j_ec_ops = {
  .name      = "j_ec",
  .nr_comp   = 3,
  .fld_names = { "jx_ec", "jy_ec", "jz_ec" },
  .run_all   = calc_j_ec,
};

// ======================================================================

#define EX_NC(ix,iy,iz) (.5f * (F( EX,ix,iy,iz) + F( EX,ix-dx,iy,iz)))
#define EY_NC(ix,iy,iz) (.5f * (F( EY,ix,iy,iz) + F( EY,ix,iy-dy,iz)))
#define EZ_NC(ix,iy,iz) (.5f * (F( EZ,ix,iy,iz) + F( EZ,ix,iy,iz-dz)))

static void
calc_E_nc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, EX, EX + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = EX_NC(ix,iy,iz);
      R(1, ix,iy,iz) = EY_NC(ix,iy,iz);
      R(2, ix,iy,iz) = EZ_NC(ix,iy,iz);
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_e_nc_ops = {
  .name      = "e_nc",
  .nr_comp   = 3,
  .fld_names = { "ex_nc", "ey_nc", "ez_nc" },
  .run_all   = calc_E_nc,
};

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
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, EX, EX + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = EX_CC(ix,iy,iz);
      R(1, ix,iy,iz) = EY_CC(ix,iy,iz);
      R(2, ix,iy,iz) = EZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_e_ops = {
  .name      = "e",
  .nr_comp   = 3,
  .fld_names = { "ex", "ey", "ez" },
  .run_all   = calc_E_cc,
};

// ======================================================================

static void
calc_E_ec(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, EX, EX + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = F(EX, ix,iy,iz);
      R(1, ix,iy,iz) = F(EY, ix,iy,iz);
      R(2, ix,iy,iz) = F(EZ, ix,iy,iz);
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_e_ec_ops = {
  .name      = "e_ec",
  .nr_comp   = 3,
  .fld_names = { "ex_ec", "ey_ec", "ez_ec" },
  .run_all   = calc_E_ec,
};

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
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, HX, HX + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = HX_NC(ix,iy,iz);
      R(1, ix,iy,iz) = HY_NC(ix,iy,iz);
      R(2, ix,iy,iz) = HZ_NC(ix,iy,iz);
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_h_nc_ops = {
  .name      = "h_nc",
  .nr_comp   = 3,
  .fld_names = { "hx_nc", "hy_nc", "hz_nc" },
  .run_all   = calc_H_nc,
};

// ======================================================================

#define HX_CC(ix,iy,iz) (.5f*(F(HX,ix,iy,iz) + F(HX,ix+dx,iy,iz)))
#define HY_CC(ix,iy,iz) (.5f*(F(HY,ix,iy,iz) + F(HY,ix,iy+dy,iz)))
#define HZ_CC(ix,iy,iz) (.5f*(F(HZ,ix,iy,iz) + F(HZ,ix,iy,iz+dz)))

static void
calc_H_cc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, HX, HX + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = HX_CC(ix,iy,iz);
      R(1, ix,iy,iz) = HY_CC(ix,iy,iz);
      R(2, ix,iy,iz) = HZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_h_ops = {
  .name      = "h",
  .nr_comp   = 3,
  .fld_names = { "hx", "hy", "hz" },
  .run_all   = calc_H_cc,
};

// ======================================================================

static void
calc_H_fc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, HX, HX + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = F(HX, ix,iy,iz);
      R(1, ix,iy,iz) = F(HY, ix,iy,iz);
      R(2, ix,iy,iz) = F(HZ, ix,iy,iz);
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_h_fc_ops = {
  .name      = "h_fc",
  .nr_comp   = 3,
  .fld_names = { "hx_fc", "hy_fc", "hz_fc" },
  .run_all   = calc_H_fc,
};

// ======================================================================

static void
calc_jdote_cc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	      struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, JXI, EX + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = JX_CC(ix,iy,iz) * EX_CC(ix,iy,iz);
      R(1, ix,iy,iz) = JY_CC(ix,iy,iz) * EY_CC(ix,iy,iz);
      R(2, ix,iy,iz) = JZ_CC(ix,iy,iz) * EZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_jdote_ops = {
  .name      = "jdote",
  .nr_comp   = 3,
  .fld_names = { "jxex", "jyey", "jzez" },
  .run_all   = calc_jdote_cc,
};

// ======================================================================

static void
calc_poyn_cc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	     struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, EX, HX + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
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
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_poyn_ops = {
  .name      = "poyn",
  .nr_comp   = 3,
  .fld_names = { "poynx", "poyny", "poynz" },
  .run_all   = calc_poyn_cc,
};

// ======================================================================

static void
calc_E2_cc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	   struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, EX, EX + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = sqr(EX_CC(ix,iy,iz));
      R(1, ix,iy,iz) = sqr(EY_CC(ix,iy,iz));
      R(2, ix,iy,iz) = sqr(EZ_CC(ix,iy,iz));
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_e2_ops = {
  .name      = "e2",
  .nr_comp   = 3,
  .fld_names = { "ex2", "ey2", "ez2" },
  .run_all   = calc_E2_cc,
};

// ======================================================================

static void
calc_H2_cc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	   struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, HX, HX + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = sqr(HX_CC(ix,iy,iz));
      R(1, ix,iy,iz) = sqr(HY_CC(ix,iy,iz));
      R(2, ix,iy,iz) = sqr(HZ_CC(ix,iy,iz));
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_h2_ops = {
  .name      = "h2",
  .nr_comp   = 3,
  .fld_names = { "hx2", "hy2", "hz2" },
  .run_all   = calc_H2_cc,
};

// ======================================================================

static void
calc_divb(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	  struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, HX, HX + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = 
	((F(HX, ix+dx,iy,iz) - F(HX, ix,iy,iz)) / ppsc->patch[p].dx[0] +
	 (F(HY, ix,iy+dy,iz) - F(HY, ix,iy,iz)) / ppsc->patch[p].dx[1] +
	 (F(HZ, ix,iy,iz+dz) - F(HZ, ix,iy,iz)) / ppsc->patch[p].dx[2]);
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_divb_ops = {
  .name      = "divb",
  .nr_comp   = 1,
  .fld_names = { "divb" },
  .run_all   = calc_divb,
};

// ======================================================================

static void
calc_divj_nc(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
	     struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  define_dxdydz(dx, dy, dz);
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, JXI, JXI + 3);
  mfields_t mf(mflds), mf_res(mres);
  for (int p = 0; p < mres->nr_patches; p++) {
    Fields F(mf[p]), R(mf_res[p]);
    psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
      R(0, ix,iy,iz) = 
	((F(JXI, ix,iy,iz) - F(JXI, ix-dx,iy,iz)) / ppsc->patch[p].dx[0] +
	 (F(JYI, ix,iy,iz) - F(JYI, ix,iy-dy,iz)) / ppsc->patch[p].dx[1] +
	 (F(JZI, ix,iy,iz) - F(JZI, ix,iy,iz-dz)) / ppsc->patch[p].dx[2]);
    } foreach_3d_end;
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

struct psc_output_fields_item_ops psc_output_fields_item_divj_ops = {
  .name      = "divj",
  .nr_comp   = 1,
  .fld_names = { "divj" },
  .run_all   = calc_divj_nc,
};

