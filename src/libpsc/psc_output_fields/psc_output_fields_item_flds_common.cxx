
#include "psc_output_fields_item_private.h"

#include "fields.hxx"
#include "fields_item.hxx"

#include <string>

// ======================================================================

#define define_dxdydz(dx, dy, dz)					\
  int dx _mrc_unused = (ppsc->grid().gdims[0] == 1) ? 0 : 1;		\
  int dy _mrc_unused = (ppsc->grid().gdims[1] == 1) ? 0 : 1;		\
  int dz _mrc_unused = (ppsc->grid().gdims[2] == 1) ? 0 : 1

// ======================================================================

template<class MF>
struct Item_dive
{
  using mfields_t = MF;
  using fields_t = typename mfields_t::fields_t;
  using Fields = Fields3d<fields_t>;
  
  constexpr static char const* name = "dive";
  constexpr static int n_comps = 1;
  static fld_names_t fld_names() { return { "dive" }; }
  
  static void set(Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = 
      ((F(EX, i,j,k) - F(EX, i-dx,j,k)) / ppsc->grid().dx[0] +
       (F(EY, i,j,k) - F(EY, i,j-dy,k)) / ppsc->grid().dx[1] +
       (F(EZ, i,j,k) - F(EZ, i,j,k-dz)) / ppsc->grid().dx[2]);
  }
};

FieldsItemFieldsOps<Item_dive<mfields_t>> psc_output_fields_item_dive_ops;


