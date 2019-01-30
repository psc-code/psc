
#pragma once

#include "fields.hxx"
#include "fields_item.hxx"

#define define_dxdydz(dx, dy, dz)			       \
  int dx _mrc_unused = (grid.isInvar(0)) ? 0 : 1;	       \
  int dy _mrc_unused = (grid.isInvar(1)) ? 0 : 1;	       \
  int dz _mrc_unused = (grid.isInvar(2)) ? 0 : 1

// ======================================================================
// Item_dive

template<typename _MfieldsState, typename _Mfields>
struct Item_dive
{
  using MfieldsState = _MfieldsState;
  using Mfields = _Mfields;
  using fields_t = typename Mfields::fields_t;
  using Fields = Fields3d<typename Mfields::fields_t>;
  using FieldsState = Fields3d<typename MfieldsState::fields_t>;
  
  constexpr static char const* name = "dive";
  constexpr static int n_comps = 1;
  static fld_names_t fld_names() { return { "dive" }; }

  static void set(const Grid_t& grid, Fields& R, FieldsState&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = ((F(EX, i,j,k) - F(EX, i-dx,j,k)) / grid.domain.dx[0] +
		   (F(EY, i,j,k) - F(EY, i,j-dy,k)) / grid.domain.dx[1] +
		   (F(EZ, i,j,k) - F(EZ, i,j,k-dz)) / grid.domain.dx[2]);
  }
};

// ======================================================================
// Item_divj

// FIXME, almost same as dive

template<typename _MfieldsState, typename _Mfields>
struct Item_divj
{
  using MfieldsState = _MfieldsState;
  using Mfields = _Mfields;
  using fields_t = typename Mfields::fields_t;
  using Fields = Fields3d<fields_t>;
  
  constexpr static char const* name = "divj";
  constexpr static int n_comps = 1;
  static fld_names_t fld_names() { return { "divj" }; }
  
  static void set(const Grid_t& grid, Fields& R, Fields&F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i,j,k) = ((F(JXI, i,j,k) - F(JXI, i-dx,j,k)) / grid.domain.dx[0] +
		   (F(JYI, i,j,k) - F(JYI, i,j-dy,k)) / grid.domain.dx[1] +
		   (F(JZI, i,j,k) - F(JZI, i,j,k-dz)) / grid.domain.dx[2]);
  }
};

#undef define_dxdydz
