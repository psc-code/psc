#include "psc_output_fields_item_private.h"

#include "fields.hxx"
#include "fields_item.hxx"
#include "psc_fields_as_c.h"
#include "psc_fields_cuda.h"

using MfieldsState_t = MfieldsStateDouble;
using Mfields_t = MfieldsC;

// ======================================================================

using Fields = Fields3d<Mfields_t::fields_t>;
using FieldsState = Fields3d<MfieldsState_t::fields_t>;

// FIXME, we're assuming that the result fields are "c" type

// ======================================================================

#define define_dxdydz(dx, dy, dz)		       \
  int dx _mrc_unused = (grid.isInvar(0)) ? 0 : 1;      \
  int dy _mrc_unused = (grid.isInvar(1)) ? 0 : 1;      \
  int dz _mrc_unused = (grid.isInvar(2)) ? 0 : 1

// ======================================================================

#define JX_NC(ix,iy,iz) (.5f * (F(JXI,ix,iy,iz) + F(JXI,ix-dx,iy,iz)))
#define JY_NC(ix,iy,iz) (.5f * (F(JYI,ix,iy,iz) + F(JYI,ix,iy-dy,iz)))
#define JZ_NC(ix,iy,iz) (.5f * (F(JZI,ix,iy,iz) + F(JZI,ix,iy,iz-dz)))

struct Item_j_nc
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "j_nc";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "jx_nc", "jy_nc", "jz_nc" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
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
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "j";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "jx", "jy", "jz" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
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
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "j_ec";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "jx_ec", "jy_ec", "jz_ec" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
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
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "e_nc";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "ex_nc", "ey_nc", "ez_nc" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
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
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "e";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "ex", "ey", "ez" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
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
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "e_ec";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "ex_ec", "ey_ec", "ez_ec" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
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
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "h_nc";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "hx_nc", "hy_nc", "hz_nc" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
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
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "h";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "hx", "hy", "hz" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
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
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "h_fc";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "hx_fc", "hy_fc", "hz_fc" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
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
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "jdote";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "jxex", "jyey", "jzez" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
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
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "poyn";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "poynx", "poyny", "poynz" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
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
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "e2";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "ex2", "ey2", "ez2" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
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
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "h2";
  constexpr static int n_comps = 3;
  static fld_names_t fld_names() { return { "hx2", "hy2", "hz2" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
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
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "divb";
  constexpr static int n_comps = 1;
  static fld_names_t fld_names() { return { "divb" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = 
      ((F(HX, i+dx,j,k) - F(HX, i,j,k)) / grid.domain.dx[0] +
       (F(HY, i,j+dy,k) - F(HY, i,j,k)) / grid.domain.dx[1] +
       (F(HZ, i,j,k+dz) - F(HZ, i,j,k)) / grid.domain.dx[2]);
  }
};

FieldsItemFieldsOps<Item_divb> psc_output_fields_item_divb_ops;

// ======================================================================

struct Item_divj
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;
  
  constexpr static char const* name = "divj";
  constexpr static int n_comps = 1;
  static fld_names_t fld_names() { return { "divj" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = 
      ((F(JXI, i,j,k) - F(JXI, i-dx,j,k)) / grid.domain.dx[0] +
       (F(JYI, i,j,k) - F(JYI, i,j-dy,k)) / grid.domain.dx[1] +
       (F(JZI, i,j,k) - F(JZI, i,j,k-dz)) / grid.domain.dx[2]);
  }
};

FieldsItemFieldsOps<Item_divj> psc_output_fields_item_divj_ops;

