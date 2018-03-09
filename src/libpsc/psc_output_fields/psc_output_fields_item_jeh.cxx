#include "psc_output_fields_item_private.h"

#include "fields.hxx"
#include "psc_fields_as_c.h"

using Fields = Fields3d<mfields_t::fields_t>;

template<typename Item>
struct ItemFields
{
  constexpr static char const* name = Item::name;
  constexpr static int n_comps = Item::n_comps; 
  static fld_names_t fld_names() { return Item::fld_names(); }
 
  static void run(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
		  struct psc_mparticles *mprts, struct psc_mfields *mres)
  {
    mfields_t mf = mflds_base->get_as<mfields_t>(JXI, JXI + 3);
    mfields_t mf_res(mres);
    for (int p = 0; p < mf_res->n_patches(); p++) {
      Fields F(mf[p]), R(mf_res[p]);
      psc_foreach_3d(ppsc, p, ix, iy, iz, 0, 0) {
	Item::set(R, F, ix,iy,iz);
      } foreach_3d_end;
    }
    mf.put_as(mflds_base, 0, 0);
  }
};

template<typename Item_t>
struct FieldsItemOps : psc_output_fields_item_ops {
  FieldsItemOps() {
    name      = Item_t::name;
    nr_comp   = Item_t::n_comps;
    fld_names = Item_t::fld_names();
    run_all   = Item_t::run;
  }
};

template<typename Item_t>
using FieldsItemFieldsOps = FieldsItemOps<ItemFields<Item_t>>;

// ======================================================================

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

struct Item_j_nc
{
  constexpr static char const* name = "j_nc";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "jx_nc", "jy_nc", "jz_nc" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = JX_NC(i,j,k);
    R(1, i,j,k) = JY_NC(i,j,k);
    R(2, i,j,k) = JZ_NC(i,j,k);
  }
};

FieldsItemFieldsOps<Item_j_nc> psc_output_fields_item_j_nc_ops;

// ======================================================================

#define JX_CC(ix,iy,iz) (.25f * (F(JXI, ix,iy,iz   ) + F(JXI, ix,iy+dy,iz   ) + \
				 F(JXI, ix,iy,iz+dz) + F(JXI, ix,iy+dy,iz+dz)))
#define JY_CC(ix,iy,iz) (.25f * (F(JYI, ix,iy,iz   ) + F(JYI, ix+dx,iy,iz   ) + \
				 F(JYI, ix,iy,iz+dz) + F(JYI, ix+dx,iy,iz+dz)))
#define JZ_CC(ix,iy,iz) (.25f * (F(JZI, ix,iy,iz   ) + F(JZI, ix+dx,iy   ,iz) + \
				 F(JZI, ix,iy+dy,iz) + F(JZI, ix+dx,iy+dy,iz)))

struct Item_j_cc
{
  constexpr static char const* name = "j";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "jx", "jy", "jz" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = JX_CC(i,j,k);
    R(1, i,j,k) = JY_CC(i,j,k);
    R(2, i,j,k) = JZ_CC(i,j,k);
  }
};

FieldsItemFieldsOps<Item_j_cc> psc_output_fields_item_j_ops;

// ======================================================================

struct Item_j_ec
{
  constexpr static char const* name = "j_ec";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "jx_ec", "jy_ec", "jz_ec" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = F(JXI, i,j,k);
    R(1, i,j,k) = F(JYI, i,j,k);
    R(2, i,j,k) = F(JZI, i,j,k);
  }
};

FieldsItemFieldsOps<Item_j_ec> psc_output_fields_item_j_ec_ops;

// ======================================================================

#define EX_NC(ix,iy,iz) (.5f * (F( EX,ix,iy,iz) + F( EX,ix-dx,iy,iz)))
#define EY_NC(ix,iy,iz) (.5f * (F( EY,ix,iy,iz) + F( EY,ix,iy-dy,iz)))
#define EZ_NC(ix,iy,iz) (.5f * (F( EZ,ix,iy,iz) + F( EZ,ix,iy,iz-dz)))

struct Item_e_nc
{
  constexpr static char const* name = "e_nc";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "ex_nc", "ey_nc", "ez_nc" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = EX_NC(i,j,k);
    R(1, i,j,k) = EY_NC(i,j,k);
    R(2, i,j,k) = EZ_NC(i,j,k);
  }
};

FieldsItemFieldsOps<Item_e_nc> psc_output_fields_item_e_nc_ops;

// ======================================================================

#define EX_CC(ix,iy,iz) (.25f * (F( EX,ix,iy,iz   ) + F( EX,ix,iy+dy,iz   ) + \
				 F( EX,ix,iy,iz+dz) + F( EX,ix,iy+dy,iz+dz)))
#define EY_CC(ix,iy,iz) (.25f * (F( EY,ix,iy,iz   ) + F( EY,ix+dx,iy,iz   ) + \
				 F( EY,ix,iy,iz+dz) + F( EY,ix+dx,iy,iz+dz)))
#define EZ_CC(ix,iy,iz) (.25f * (F( EZ,ix,iy   ,iz) + F( EZ,ix+dx,iy   ,iz) + \
				 F( EZ,ix,iy+dy,iz) + F( EZ,ix+dx,iy+dy,iz)))

struct Item_e_cc
{
  constexpr static char const* name = "e";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "ex", "ey", "ez" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = EX_CC(i,j,k);
    R(1, i,j,k) = EY_CC(i,j,k);
    R(2, i,j,k) = EZ_CC(i,j,k);
  }
};

FieldsItemFieldsOps<Item_e_cc> psc_output_fields_item_e_ops;

// ======================================================================

struct Item_e_ec
{
  constexpr static char const* name = "e_ec";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "ex_ec", "ey_ec", "ez_ec" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = F(EX, i,j,k);
    R(1, i,j,k) = F(EY, i,j,k);
    R(2, i,j,k) = F(EZ, i,j,k);
  }
};

FieldsItemFieldsOps<Item_e_ec> psc_output_fields_item_e_ec_ops;

// ======================================================================

#define HX_NC(ix,iy,iz) (.25f*(F(HX,ix,iy,iz   ) + F(HX,ix,iy-dy,iz   ) + \
			       F(HX,ix,iy,iz-dz) + F(HX,ix,iy-dy,iz-dz)))
#define HY_NC(ix,iy,iz) (.25f*(F(HY,ix,iy,iz   ) + F(HY,ix-dx,iy,iz   ) + \
			       F(HY,ix,iy,iz-dz) + F(HY,ix-dx,iy,iz-dz)))
#define HZ_NC(ix,iy,iz) (.25f*(F(HZ,ix,iy   ,iz) + F(HZ,ix-dx,iy   ,iz) + \
			       F(HZ,ix,iy-dy,iz) + F(HZ,ix-dx,iy-dy,iz)))

struct Item_h_nc
{
  constexpr static char const* name = "h_nc";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "hx_nc", "hy_nc", "hz_nc" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = HX_NC(i,j,k);
    R(1, i,j,k) = HY_NC(i,j,k);
    R(2, i,j,k) = HZ_NC(i,j,k);
  }
};

FieldsItemFieldsOps<Item_h_nc> psc_output_fields_item_h_nc_ops;

// ======================================================================

#define HX_CC(ix,iy,iz) (.5f*(F(HX,ix,iy,iz) + F(HX,ix+dx,iy,iz)))
#define HY_CC(ix,iy,iz) (.5f*(F(HY,ix,iy,iz) + F(HY,ix,iy+dy,iz)))
#define HZ_CC(ix,iy,iz) (.5f*(F(HZ,ix,iy,iz) + F(HZ,ix,iy,iz+dz)))

struct Item_h_cc
{
  constexpr static char const* name = "h";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "hx", "hy", "hz" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = HX_CC(i,j,k);
    R(1, i,j,k) = HY_CC(i,j,k);
    R(2, i,j,k) = HZ_CC(i,j,k);
  }
};

FieldsItemFieldsOps<Item_h_cc> psc_output_fields_item_h_ops;

// ======================================================================

struct Item_h_fc
{
  constexpr static char const* name = "h_fc";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "hx_fc", "hy_fc", "hz_fc" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = F(HX, i,j,k);
    R(1, i,j,k) = F(HY, i,j,k);
    R(2, i,j,k) = F(HZ, i,j,k);
  }
};

FieldsItemFieldsOps<Item_h_fc> psc_output_fields_item_h_fc_ops;

// ======================================================================

struct Item_jdote
{
  constexpr static char const* name = "jdote";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "jxex", "jyey", "jzez" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = JX_CC(i,j,k) * EX_CC(i,j,k);
    R(1, i,j,k) = JY_CC(i,j,k) * EY_CC(i,j,k);
    R(2, i,j,k) = JZ_CC(i,j,k) * EZ_CC(i,j,k);
  }
};

FieldsItemFieldsOps<Item_jdote> psc_output_fields_item_jdote_ops;

// ======================================================================

struct Item_poyn
{
  constexpr static char const* name = "poyn";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "poynx", "poyny", "poynz" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = (EY_CC(i,j,k) * HZ_CC(i,j,k) - 
		   EZ_CC(i,j,k) * HY_CC(i,j,k));
    R(1, i,j,k) = (EZ_CC(i,j,k) * HX_CC(i,j,k) -
		   EX_CC(i,j,k) * HZ_CC(i,j,k));
    R(2, i,j,k) = (EX_CC(i,j,k) * HY_CC(i,j,k) -
		   EY_CC(i,j,k) * HX_CC(i,j,k));
  }
};

FieldsItemFieldsOps<Item_poyn> psc_output_fields_item_poyn_ops;

// ======================================================================

struct Item_e2
{
  constexpr static char const* name = "e2";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "ex2", "ey2", "ez2" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = sqr(EX_CC(i,j,k));
    R(1, i,j,k) = sqr(EY_CC(i,j,k));
    R(2, i,j,k) = sqr(EZ_CC(i,j,k));
  }
};

FieldsItemFieldsOps<Item_e2> psc_output_fields_item_e2_ops;

// ======================================================================

struct Item_h2
{
  constexpr static char const* name = "h2";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "hx2", "hy2", "hz2" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = sqr(HX_CC(i,j,k));
    R(1, i,j,k) = sqr(HY_CC(i,j,k));
    R(2, i,j,k) = sqr(HZ_CC(i,j,k));
  }
};

FieldsItemFieldsOps<Item_h2> psc_output_fields_item_h2_ops;

// ======================================================================

struct Item_divb
{
  constexpr static char const* name = "divb";
  constexpr static int n_comps = 1;
  static fld_names_t fld_names() { return { "divb" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = 
      ((F(HX, i+dx,j,k) - F(HX, i,j,k)) / ppsc->grid().dx[0] +
       (F(HY, i,j+dy,k) - F(HY, i,j,k)) / ppsc->grid().dx[1] +
       (F(HZ, i,j,k+dz) - F(HZ, i,j,k)) / ppsc->grid().dx[2]);
  }
};

FieldsItemFieldsOps<Item_divb> psc_output_fields_item_divb_ops;

// ======================================================================

struct Item_divj
{
  constexpr static char const* name = "divj";
  constexpr static int n_comps = 1;
  static fld_names_t fld_names() { return { "divj" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = 
      ((F(JXI, i,j,k) - F(JXI, i-dx,j,k)) / ppsc->grid().dx[0] +
       (F(JYI, i,j,k) - F(JYI, i,j-dy,k)) / ppsc->grid().dx[1] +
       (F(JZI, i,j,k) - F(JZI, i,j,k-dz)) / ppsc->grid().dx[2]);
  }
};

FieldsItemFieldsOps<Item_divj> psc_output_fields_item_divj_ops;

