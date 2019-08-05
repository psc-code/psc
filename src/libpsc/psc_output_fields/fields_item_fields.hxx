
#pragma once

#include "fields_item.hxx"
#include "psc_fields_c.h"

template <typename E>
struct isSpaceCuda : std::false_type
{};

// ======================================================================
// evalMfields

template <typename E,
          typename std::enable_if<!isSpaceCuda<E>::value, int>::type = 0>
MfieldsC evalMfields(const MFexpression<E>& xp)
{
  const auto& exp = xp.derived();
  MfieldsC mflds{exp.grid(), exp.n_comps(), exp.ibn()};

  for (int p = 0; p < mflds.n_patches(); p++) {
    auto flds = mflds[p];
    for (int m = 0; m < exp.n_comps(); m++) {
      mflds.Foreach_3d(0, 0, [&](int i, int j, int k) {
        flds(m, i, j, k) = exp(m, {i, j, k}, p);
      });
    }
  }
  return mflds;
}

#ifdef USE_CUDA

template <typename E,
          typename std::enable_if<isSpaceCuda<E>::value, int>::type = 0>
MfieldsC evalMfields(const MFexpression<E>& xp)
{
  const auto& exp = xp.derived().result();
  MfieldsC mflds{exp.grid(), exp.n_comps(), exp.ibn()};

  auto h_exp = hostMirror(exp);
  copy(exp, h_exp);

  for (int p = 0; p < mflds.n_patches(); p++) {
    auto flds = mflds[p];
    for (int m = 0; m < exp.n_comps(); m++) {
      mflds.Foreach_3d(0, 0, [&](int i, int j, int k) {
        flds(m, i, j, k) = h_exp[p](m, i, j, k);
      });
    }
  }
  return mflds;
}

#endif

// ======================================================================
// AdaptMfields

template <typename EXP>
class AdaptMfields
{
  class Patch;

public:
  explicit AdaptMfields(EXP& exp) : exp_{exp} {}

  Patch operator[](int p) const { return {*this, p}; }

  int n_patches() const { return exp_.grid().n_patches(); }

private:
  EXP& exp_;
};

template <typename EXP>
class AdaptMfields<EXP>::Patch
{
public:
  using Real = typename EXP::Real;

  Patch(const AdaptMfields& parent, int p) : parent_{parent}, p_{p} {}

  Real operator()(int m, int i, int j, int k) const
  {
    return parent_.exp_(m, {i, j, k}, p_);
  }

private:
  const AdaptMfields<EXP>& parent_;
  int p_;
};

template <typename E>
AdaptMfields<E> adaptMfields(MFexpression<E>& xp)
{
  return AdaptMfields<E>{xp.derived()};
}

// ======================================================================

using MfieldsState_t = MfieldsStateDouble;
using Mfields_t = MfieldsC;

// ======================================================================

#define define_dxdydz(dx, dy, dz)                                              \
  int dx _mrc_unused = (grid.isInvar(0)) ? 0 : 1;                              \
  int dy _mrc_unused = (grid.isInvar(1)) ? 0 : 1;                              \
  int dz _mrc_unused = (grid.isInvar(2)) ? 0 : 1

// ======================================================================

#define JX_NC(ix, iy, iz) (.5f * (F(JXI, ix, iy, iz) + F(JXI, ix - dx, iy, iz)))
#define JY_NC(ix, iy, iz) (.5f * (F(JYI, ix, iy, iz) + F(JYI, ix, iy - dy, iz)))
#define JZ_NC(ix, iy, iz) (.5f * (F(JZI, ix, iy, iz) + F(JZI, ix, iy, iz - dz)))

struct Item_j_nc
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "j_nc";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names()
  {
    return {"jx_nc", "jy_nc", "jz_nc"};
  }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) = JX_NC(i, j, k);
    R(1, i, j, k) = JY_NC(i, j, k);
    R(2, i, j, k) = JZ_NC(i, j, k);
  }
};

// ======================================================================

#define JX_CC(ix, iy, iz)                                                      \
  (.25f * (F(JXI, ix, iy, iz) + F(JXI, ix, iy + dy, iz) +                      \
           F(JXI, ix, iy, iz + dz) + F(JXI, ix, iy + dy, iz + dz)))
#define JY_CC(ix, iy, iz)                                                      \
  (.25f * (F(JYI, ix, iy, iz) + F(JYI, ix + dx, iy, iz) +                      \
           F(JYI, ix, iy, iz + dz) + F(JYI, ix + dx, iy, iz + dz)))
#define JZ_CC(ix, iy, iz)                                                      \
  (.25f * (F(JZI, ix, iy, iz) + F(JZI, ix + dx, iy, iz) +                      \
           F(JZI, ix, iy + dy, iz) + F(JZI, ix + dx, iy + dy, iz)))

struct Item_j_cc
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "j";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names() { return {"jx", "jy", "jz"}; }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) = JX_CC(i, j, k);
    R(1, i, j, k) = JY_CC(i, j, k);
    R(2, i, j, k) = JZ_CC(i, j, k);
  }
};

// ======================================================================

struct Item_j_ec
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "j_ec";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names()
  {
    return {"jx_ec", "jy_ec", "jz_ec"};
  }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) = F(JXI, i, j, k);
    R(1, i, j, k) = F(JYI, i, j, k);
    R(2, i, j, k) = F(JZI, i, j, k);
  }
};

// ======================================================================

#define EX_NC(ix, iy, iz) (.5f * (F(EX, ix, iy, iz) + F(EX, ix - dx, iy, iz)))
#define EY_NC(ix, iy, iz) (.5f * (F(EY, ix, iy, iz) + F(EY, ix, iy - dy, iz)))
#define EZ_NC(ix, iy, iz) (.5f * (F(EZ, ix, iy, iz) + F(EZ, ix, iy, iz - dz)))

struct Item_e_nc
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "e_nc";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names()
  {
    return {"ex_nc", "ey_nc", "ez_nc"};
  }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) = EX_NC(i, j, k);
    R(1, i, j, k) = EY_NC(i, j, k);
    R(2, i, j, k) = EZ_NC(i, j, k);
  }
};

// ======================================================================

#define EX_CC(ix, iy, iz)                                                      \
  (.25f * (F(EX, ix, iy, iz) + F(EX, ix, iy + dy, iz) +                        \
           F(EX, ix, iy, iz + dz) + F(EX, ix, iy + dy, iz + dz)))
#define EY_CC(ix, iy, iz)                                                      \
  (.25f * (F(EY, ix, iy, iz) + F(EY, ix + dx, iy, iz) +                        \
           F(EY, ix, iy, iz + dz) + F(EY, ix + dx, iy, iz + dz)))
#define EZ_CC(ix, iy, iz)                                                      \
  (.25f * (F(EZ, ix, iy, iz) + F(EZ, ix + dx, iy, iz) +                        \
           F(EZ, ix, iy + dy, iz) + F(EZ, ix + dx, iy + dy, iz)))

struct Item_e_cc
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "e";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names() { return {"ex", "ey", "ez"}; }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) = EX_CC(i, j, k);
    R(1, i, j, k) = EY_CC(i, j, k);
    R(2, i, j, k) = EZ_CC(i, j, k);
  }
};

// ======================================================================

struct Item_e_ec
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "e_ec";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names()
  {
    return {"ex_ec", "ey_ec", "ez_ec"};
  }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) = F(EX, i, j, k);
    R(1, i, j, k) = F(EY, i, j, k);
    R(2, i, j, k) = F(EZ, i, j, k);
  }
};

// ======================================================================

#define HX_NC(ix, iy, iz)                                                      \
  (.25f * (F(HX, ix, iy, iz) + F(HX, ix, iy - dy, iz) +                        \
           F(HX, ix, iy, iz - dz) + F(HX, ix, iy - dy, iz - dz)))
#define HY_NC(ix, iy, iz)                                                      \
  (.25f * (F(HY, ix, iy, iz) + F(HY, ix - dx, iy, iz) +                        \
           F(HY, ix, iy, iz - dz) + F(HY, ix - dx, iy, iz - dz)))
#define HZ_NC(ix, iy, iz)                                                      \
  (.25f * (F(HZ, ix, iy, iz) + F(HZ, ix - dx, iy, iz) +                        \
           F(HZ, ix, iy - dy, iz) + F(HZ, ix - dx, iy - dy, iz)))

struct Item_h_nc
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "h_nc";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names()
  {
    return {"hx_nc", "hy_nc", "hz_nc"};
  }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) = HX_NC(i, j, k);
    R(1, i, j, k) = HY_NC(i, j, k);
    R(2, i, j, k) = HZ_NC(i, j, k);
  }
};

// ======================================================================

#define HX_CC(ix, iy, iz) (.5f * (F(HX, ix, iy, iz) + F(HX, ix + dx, iy, iz)))
#define HY_CC(ix, iy, iz) (.5f * (F(HY, ix, iy, iz) + F(HY, ix, iy + dy, iz)))
#define HZ_CC(ix, iy, iz) (.5f * (F(HZ, ix, iy, iz) + F(HZ, ix, iy, iz + dz)))

struct Item_h_cc
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "h";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names() { return {"hx", "hy", "hz"}; }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) = HX_CC(i, j, k);
    R(1, i, j, k) = HY_CC(i, j, k);
    R(2, i, j, k) = HZ_CC(i, j, k);
  }
};

// ======================================================================

struct Item_h_fc
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "h_fc";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names()
  {
    return {"hx_fc", "hy_fc", "hz_fc"};
  }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) = F(HX, i, j, k);
    R(1, i, j, k) = F(HY, i, j, k);
    R(2, i, j, k) = F(HZ, i, j, k);
  }
};

// ======================================================================

struct Item_jdote
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "jdote";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names()
  {
    return {"jxex", "jyey", "jzez"};
  }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) = JX_CC(i, j, k) * EX_CC(i, j, k);
    R(1, i, j, k) = JY_CC(i, j, k) * EY_CC(i, j, k);
    R(2, i, j, k) = JZ_CC(i, j, k) * EZ_CC(i, j, k);
  }
};

// ======================================================================

struct Item_poyn
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "poyn";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names()
  {
    return {"poynx", "poyny", "poynz"};
  }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) =
      (EY_CC(i, j, k) * HZ_CC(i, j, k) - EZ_CC(i, j, k) * HY_CC(i, j, k));
    R(1, i, j, k) =
      (EZ_CC(i, j, k) * HX_CC(i, j, k) - EX_CC(i, j, k) * HZ_CC(i, j, k));
    R(2, i, j, k) =
      (EX_CC(i, j, k) * HY_CC(i, j, k) - EY_CC(i, j, k) * HX_CC(i, j, k));
  }
};

// ======================================================================

struct Item_e2
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "e2";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names() { return {"ex2", "ey2", "ez2"}; }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) = sqr(EX_CC(i, j, k));
    R(1, i, j, k) = sqr(EY_CC(i, j, k));
    R(2, i, j, k) = sqr(EZ_CC(i, j, k));
  }
};

// ======================================================================

struct Item_h2
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "h2";
  constexpr static int n_comps = 3;
  static std::vector<std::string> fld_names() { return {"hx2", "hy2", "hz2"}; }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) = sqr(HX_CC(i, j, k));
    R(1, i, j, k) = sqr(HY_CC(i, j, k));
    R(2, i, j, k) = sqr(HZ_CC(i, j, k));
  }
};

// ======================================================================

struct Item_divb
{
  using MfieldsState = MfieldsState_t;
  using Mfields = Mfields_t;

  constexpr static char const* name = "divb";
  constexpr static int n_comps = 1;
  static std::vector<std::string> fld_names() { return {"divb"}; }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) =
      ((F(HX, i + dx, j, k) - F(HX, i, j, k)) / grid.domain.dx[0] +
       (F(HY, i, j + dy, k) - F(HY, i, j, k)) / grid.domain.dx[1] +
       (F(HZ, i, j, k + dz) - F(HZ, i, j, k)) / grid.domain.dx[2]);
  }
};

// ======================================================================
// Item_jeh
//
// Main fields in their natural staggering

template <typename MfieldsState>
class Item_jeh : public MFexpression<Item_jeh<MfieldsState>>
{
public:
  using Real = typename MfieldsState::real_t;

  static char const* name() { return "jeh"; }
  static int n_comps() { return 9; }
  static std::vector<std::string> comp_names()
  {
    return {"jx_ec", "jy_ec", "jz_ec", "ex_ec", "ey_ec",
            "ez_ec", "hx_fc", "hy_fc", "hz_fc"};
  }

  Item_jeh(MfieldsState& mflds) : mflds_{mflds} {}

  Real operator()(int m, Int3 ijk, int p) const
  {
    switch (m) {
      case 0: return mflds_[p](JXI, ijk[0], ijk[1], ijk[2]);
      case 1: return mflds_[p](JYI, ijk[0], ijk[1], ijk[2]);
      case 2: return mflds_[p](JZI, ijk[0], ijk[1], ijk[2]);
      case 3: return mflds_[p](EX, ijk[0], ijk[1], ijk[2]);
      case 4: return mflds_[p](EY, ijk[0], ijk[1], ijk[2]);
      case 5: return mflds_[p](EZ, ijk[0], ijk[1], ijk[2]);
      case 6: return mflds_[p](HX, ijk[0], ijk[1], ijk[2]);
      case 7: return mflds_[p](HY, ijk[0], ijk[1], ijk[2]);
      case 8: return mflds_[p](HZ, ijk[0], ijk[1], ijk[2]);
      default: std::abort();
    }
  }

  const Grid_t& grid() const { return mflds_.grid(); }
  Int3 ibn() const { return {}; }
  int n_patches() const { return grid().n_patches(); }

private:
  MfieldsState& mflds_;
};

// ======================================================================
// Item_dive

template <typename MfieldsState>
class Item_dive : public MFexpression<Item_dive<MfieldsState>>
{
public:
  using Real = typename MfieldsState::real_t;

  static char const* name() { return "dive"; }
  static int n_comps() { return 1; }
  static std::vector<std::string> comp_names() { return {"dive"}; }

  Item_dive(MfieldsState& mflds) : mflds_{mflds} {}

  Real operator()(int m, Int3 ijk, int p) const
  {
    const auto& grid = mflds_.grid();
    define_dxdydz(dx, dy, dz);

    return ((mflds_[p](EX, ijk[0], ijk[1], ijk[2]) -
             mflds_[p](EX, ijk[0] - dx, ijk[1], ijk[2])) /
              grid.domain.dx[0] +
            (mflds_[p](EY, ijk[0], ijk[1], ijk[2]) -
             mflds_[p](EY, ijk[0], ijk[1] - dy, ijk[2])) /
              grid.domain.dx[1] +
            (mflds_[p](EZ, ijk[0], ijk[1], ijk[2]) -
             mflds_[p](EZ, ijk[0], ijk[1], ijk[2] - dz)) /
              grid.domain.dx[2]);
  }

  const Grid_t& grid() const { return mflds_.grid(); }
  Int3 ibn() const { return {}; }
  int n_patches() const { return grid().n_patches(); }

private:
  MfieldsState& mflds_;
};

// ======================================================================
// Item_divj

// FIXME, almost same as dive

template <typename _MfieldsState, typename _Mfields>
struct Item_divj
{
  using MfieldsState = _MfieldsState;
  using Mfields = _Mfields;

  constexpr static char const* name = "divj";
  constexpr static int n_comps = 1;
  static std::vector<std::string> fld_names() { return {"divj"}; }

  template <typename FE>
  static void set(const Grid_t& grid, FE& R, FE& F, int i, int j, int k)
  {
    define_dxdydz(dx, dy, dz);
    R(0, i, j, k) =
      ((F(JXI, i, j, k) - F(JXI, i - dx, j, k)) / grid.domain.dx[0] +
       (F(JYI, i, j, k) - F(JYI, i, j - dy, k)) / grid.domain.dx[1] +
       (F(JZI, i, j, k) - F(JZI, i, j, k - dz)) / grid.domain.dx[2]);
  }
};

#undef define_dxdydz
