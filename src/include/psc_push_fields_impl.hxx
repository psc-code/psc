
#ifndef PSC_PUSH_FIELS_IMPL_HXX
#define PSC_PUSH_FIELS_IMPL_HXX

#include "fields.hxx"
#include "fields_traits.hxx"

#include "push_fields.hxx"
#include "psc.h" // FIXME, for foreach_3d macro

// ----------------------------------------------------------------------
// Foreach_3d

template<class F>
static void Foreach_3d(const Grid_t& grid, F f, int l, int r)
{
  grid.Foreach_3d(l, r, [&](int i, int j, int k) {
      f.x(i,j,k);
      f.y(i,j,k);
      f.z(i,j,k);
    });
}

// ----------------------------------------------------------------------

template<typename Fields>
class PushBase
{
public:
  using real_t = typename Fields::real_t;
  using fields_t = typename Fields::fields_t;
  using dim = typename Fields::dim;
  
  PushBase(const Grid_t& grid, double dt_fac)
  {
    dth = dt_fac * grid.dt;

    // FIXME, it'd be even better to not even calculate derivates
    // that will be multiplied by 0 
    cnx = dim::InvarX::value ? 0 : dth / grid.domain.dx[0];
    cny = dim::InvarY::value ? 0 : dth / grid.domain.dx[1];
    cnz = dim::InvarZ::value ? 0 : dth / grid.domain.dx[2];
  }

protected:
  real_t dth;
  real_t cnx, cny, cnz;
};
  
template<typename Fields>
class PushE : PushBase<Fields>
{
public:
  using Base = PushBase<Fields>;
  using typename Base::real_t;
  using typename Base::fields_t;
  
  PushE(const fields_t& flds, double dt_fac)
    : Base(flds.grid(), dt_fac),
      F(flds)
  {
  }
  
  void x(int i, int j, int k)
  {
    F(EX, i,j,k) += (cny * (F(HZ, i,j,k) - F(HZ, i,j-1,k)) - cnz * (F(HY, i,j,k) - F(HY, i,j,k-1)) -
    		     dth * F(JXI, i,j,k));
  }

  void y(int i, int j, int k)
  {
    F(EY, i,j,k) += (cnz * (F(HX, i,j,k) - F(HX, i,j,k-1)) - cnx * (F(HZ, i,j,k) - F(HZ, i-1,j,k)) -
		     dth * F(JYI, i,j,k));
  }

  void z(int i, int j, int k)
  {
    F(EZ, i,j,k) += (cnx * (F(HY, i,j,k) - F(HY, i-1,j,k)) - cny * (F(HX, i,j,k) - F(HX, i,j-1,k)) -
		     dth * F(JZI, i,j,k));
  }

protected:
  Fields F;
  using Base::dth;
  using Base::cnx;
  using Base::cny;
  using Base::cnz;
};

template<typename Fields>
class PushH : PushBase<Fields>
{
public:
  using Base = PushBase<Fields>;
  using typename Base::real_t;
  using typename Base::fields_t;
  
  PushH(const fields_t& flds, double dt_fac)
    : Base(flds.grid(), dt_fac),
      F(flds)
  {}
  
  void x(int i, int j, int k)
  {
    F(HX, i,j,k) -= (cny * (F(EZ, i,j+1,k) - F(EZ, i,j,k)) - cnz * (F(EY, i,j,k+1) - F(EY, i,j,k)));
  }

  void y(int i, int j, int k)
  {
    F(HY, i,j,k) -= (cnz * (F(EX, i,j,k+1) - F(EX, i,j,k)) - cnx * (F(EZ, i+1,j,k) - F(EZ, i,j,k)));
  }

  void z(int i, int j, int k)
  {
    F(HZ, i,j,k) -= (cnx * (F(EY, i+1,j,k) - F(EY, i,j,k)) - cny * (F(EX, i,j+1,k) - F(EX, i,j,k)));
  }

protected:
  Fields F;
  using Base::dth;
  using Base::cnx;
  using Base::cny;
  using Base::cnz;
};

// ======================================================================
// class PushFields

template<typename _MfieldsState>
class PushFields : public PushFieldsBase
{
  using MfieldsState = _MfieldsState;

public:
  // ----------------------------------------------------------------------
  // push_E

  template<typename dim>
  void push_E(MfieldsState& mflds, double dt_fac, dim tag)
  {
    using Fields = Fields3d<typename MfieldsState::fields_t, dim>;
    
    for (int p = 0; p < mflds.n_patches(); p++) {
      PushE<Fields> push_E(mflds[p], dt_fac);
      Foreach_3d(mflds.grid(), push_E, 1, 2);
    }
  }
  
  // ----------------------------------------------------------------------
  // push_H

  template<typename dim>
  void push_H(MfieldsState& mflds, double dt_fac, dim tag)
  {
    using Fields = Fields3d<typename MfieldsState::fields_t, dim>;

    for (int p = 0; p < mflds.n_patches(); p++) {
      PushH<Fields> push_H(mflds[p], dt_fac);
      Foreach_3d(mflds.grid(), push_H, 2, 1);
    }
  }

  // ----------------------------------------------------------------------
  // push_E
  //
  // E-field propagation E^(n)    , H^(n), j^(n) 
  //                  -> E^(n+0.5), H^(n), j^(n)
  // Ex^{n}[-.5:+.5][-1:1][-1:1] -> Ex^{n+.5}[-.5:+.5][-1:1][-1:1]
  // using Hx^{n}[-1:1][-1.5:1.5][-1.5:1.5]
  //       jx^{n+1}[-.5:.5][-1:1][-1:1]
  
  void push_E(MfieldsStateBase& mflds_base, double dt_fac) override
  {
    auto& mflds = mflds_base.get_as<MfieldsState>(JXI, HX + 3);
    
    const auto& grid = mflds.grid();
    using Bool3 = Vec3<bool>;
    Bool3 invar{grid.isInvar(0), grid.isInvar(1), grid.isInvar(2)};

    if (invar == Bool3{false, false, false}) {
      push_E(mflds, dt_fac, dim_xyz{});
    } else if (invar == Bool3{true, false, false}) {
      push_E(mflds, dt_fac, dim_yz{});
    } else if (invar == Bool3{false, true, false}) {
      push_E(mflds, dt_fac, dim_xz{});
    } else {
      assert(0);
    }

    mflds_base.put_as(mflds, EX, EX + 3);
  }

  // ----------------------------------------------------------------------
  // push_H
  //
  // B-field propagation E^(n+0.5), H^(n    ), j^(n), m^(n+0.5)
  //                  -> E^(n+0.5), H^(n+0.5), j^(n), m^(n+0.5)
  // Hx^{n}[:][-.5:.5][-.5:.5] -> Hx^{n+.5}[:][-.5:.5][-.5:.5]
  // using Ex^{n+.5}[-.5:+.5][-1:1][-1:1]

  void push_H(MfieldsStateBase& mflds_base, double dt_fac) override
  {
    auto& mflds = mflds_base.get_as<MfieldsState>(JXI, HX + 3);
    
    const auto& grid = mflds.grid();
    using Bool3 = Vec3<bool>;
    Bool3 invar{grid.isInvar(0), grid.isInvar(1), grid.isInvar(2)};

    if (invar == Bool3{false, false, false}) {
      push_H(mflds, dt_fac, dim_xyz{});
    } else if (invar == Bool3{true, false, false}) {
      push_H(mflds, dt_fac, dim_yz{});
    } else if (invar == Bool3{false, true, false}) {
      push_H(mflds, dt_fac, dim_xz{});
    } else {
      assert(0);
    }

    mflds_base.put_as(mflds, HX, HX + 3);
  }    
};

#endif

