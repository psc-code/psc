
#ifndef FIELDS3D_HXX
#define FIELDS3D_HXX

#include "psc.h"

#include "grid.hxx"
#include <kg/SArrayView.h>

#include <mrc_profile.h>
#include <psc/gtensor.h>

#include <type_traits>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <typeindex>
#include <list>
#include <string>

// ======================================================================
// SArrayView

namespace detail
{

GT_INLINE gt::shape_type<4> strides(const Int3& im, int n_comps)
{
  return gt::shape(1, im[0], im[0] * im[1], im[0] * im[1] * im[2]);
}

} // namespace detail

template <typename T, typename S>
struct SArrayView
{
  using Storage = gt::gtensor_span<T, 4, S>;
  using value_type = typename Storage::value_type;
  using reference = typename Storage::reference;
  using const_reference = typename Storage::const_reference;
  using pointer = typename Storage::pointer;
  using const_pointer = typename Storage::const_pointer;

  KG_INLINE SArrayView(const Int3& ib, const Storage& storage)
    : ib_(ib), storage_(storage)
  {}

  KG_INLINE const Int3& ib() const { return ib_; }
  KG_INLINE Storage& storage() { return storage_; }
  KG_INLINE const Storage& storage() const { return storage_; }

private:
  Storage storage_;
  Int3 ib_; //> lower bounds per direction
};

template <typename Derived>
class MFexpression
{
public:
  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }
};

// ======================================================================
// MfieldsBase

struct MfieldsBase
{
  using convert_func_t = void (*)(MfieldsBase&, MfieldsBase&, int, int);
  using Convert = std::unordered_map<std::type_index, convert_func_t>;

  MfieldsBase(const Grid_t& grid, int n_fields, Int3 ibn)
    : grid_(&grid), n_fields_(n_fields), ibn_(ibn)
  {
    instances.push_back(this);
  }

  MfieldsBase(const MfieldsBase&) = delete;
  MfieldsBase& operator=(const MfieldsBase&) = delete;

  MfieldsBase(MfieldsBase&& o) : MfieldsBase(*o.grid_, o.n_fields_, o.ibn_) {}

  MfieldsBase& operator=(MfieldsBase&& o) = default;

  virtual ~MfieldsBase() { instances.remove(this); }

  virtual void reset(const Grid_t& grid) { grid_ = &grid; }

  int _n_comps() const { return n_fields_; }
  Int3 ibn() const { return ibn_; }

  const Grid_t& _grid() const { return *grid_; }

  template <typename MF>
  MF& get_as(int mb, int me)
  {
    // If we're already the subtype, nothing to be done
    if (typeid(*this) == typeid(MF)) {
      return *dynamic_cast<MF*>(this);
    }

    static int pr;
    if (!pr) {
      pr = prof_register("Mfields_get_as", 1., 0, 0);
    }
    prof_start(pr);

    // mprintf("get_as %s (%s) %d %d\n", type, psc_mfields_type(mflds_base), mb,
    // me);

    auto& mflds = *new MF{_grid(), n_fields_, ibn()};

    MfieldsBase::convert(*this, mflds, mb, me);

    prof_stop(pr);
    return mflds;
  }

  template <typename MF>
  void put_as(MF& mflds, int mb, int me)
  {
    // If we're already the subtype, nothing to be done
    if (typeid(*this) == typeid(mflds)) {
      return;
    }

    static int pr;
    if (!pr) {
      pr = prof_register("Mfields_put_as", 1., 0, 0);
    }
    prof_start(pr);

    MfieldsBase::convert(mflds, *this, mb, me);
    delete &mflds;

    prof_stop(pr);
  }

  virtual const Convert& convert_to()
  {
    static const Convert convert_to_;
    return convert_to_;
  }
  virtual const Convert& convert_from()
  {
    static const Convert convert_from_;
    return convert_from_;
  }
  static void convert(MfieldsBase& mf_from, MfieldsBase& mf_to, int mb, int me);

  static std::list<MfieldsBase*> instances;

private:
  int n_fields_;

protected:
  const Grid_t* grid_;
  Int3 ibn_;
};

// ======================================================================
// MfieldsStateBase

struct MfieldsStateBase
{
  using convert_func_t = void (*)(MfieldsStateBase&, MfieldsStateBase&, int,
                                  int);
  using Convert = std::unordered_map<std::type_index, convert_func_t>;

  MfieldsStateBase(const Grid_t& grid, int n_fields, Int3 ibn)
    : grid_(&grid), n_fields_(n_fields), ibn_(ibn)
  {
    instances.push_back(this);
  }

  virtual ~MfieldsStateBase() { instances.remove(this); }

  virtual void reset(const Grid_t& grid) { grid_ = &grid; }

  int _n_patches() const { return grid_->n_patches(); }
  int _n_comps() const { return n_fields_; }
  Int3 ibn() const { return ibn_; }

  const Grid_t& _grid() { return *grid_; }

  virtual const Convert& convert_to()
  {
    static const Convert convert_to_;
    return convert_to_;
  }
  virtual const Convert& convert_from()
  {
    static const Convert convert_from_;
    return convert_from_;
  }
  static void convert(MfieldsStateBase& mf_from, MfieldsStateBase& mf_to,
                      int mb, int me);

  static std::list<MfieldsStateBase*> instances;

  template <typename MF>
  MF& get_as(int mb, int me)
  {
    // If we're already the subtype, nothing to be done
    if (typeid(*this) == typeid(MF)) {
      return *dynamic_cast<MF*>(this);
    }

    static int pr;
    if (!pr) {
      pr = prof_register("Mfields_get_as", 1., 0, 0);
    }
    prof_start(pr);

    // mprintf("get_as %s (%s) %d %d\n", type, psc_mfields_type(mflds_base), mb,
    // me);

    auto& mflds = *new MF{_grid()};

    MfieldsStateBase::convert(*this, mflds, mb, me);

    prof_stop(pr);
    return mflds;
  }

  template <typename MF>
  void put_as(MF& mflds, int mb, int me)
  {
    // If we're already the subtype, nothing to be done
    if (typeid(*this) == typeid(mflds)) {
      return;
    }

    static int pr;
    if (!pr) {
      pr = prof_register("Mfields_put_as", 1., 0, 0);
    }
    prof_start(pr);

    MfieldsStateBase::convert(mflds, *this, mb, me);
    delete &mflds;

    prof_stop(pr);
  }

protected:
  int n_fields_;
  const Grid_t* grid_;
  Int3 ibn_;
};

// ======================================================================
// MfieldsCRTP

template <typename C>
struct MfieldsCRTPInnerTypes;

template <typename D>
class MfieldsCRTP
{
public:
  using Derived = D;

  using InnerTypes = MfieldsCRTPInnerTypes<D>;
  using Storage = typename InnerTypes::Storage;
  using Real = typename Storage::value_type;
  using fields_view_t = SArrayView<Real, typename Storage::space_type>;

  KG_INLINE const kg::Box3& box() const { return box_; }
  KG_INLINE const Int3& ib() const { return box_.ib(); }
  KG_INLINE const Int3& im() const { return box_.im(); }
  KG_INLINE int ib(int d) const { return box_.ib(d); }
  KG_INLINE int im(int d) const { return box_.im(d); }
  KG_INLINE int n_comps() const { return n_fields_; }
  KG_INLINE int n_patches() const { return n_patches_; }

  MfieldsCRTP(int n_fields, const kg::Box3& box, int n_patches)
    : n_fields_(n_fields), box_{box}, n_patches_{n_patches}
  {}

  void reset(int n_patches)
  {
    n_patches_ = n_patches;
    storage().resize(
      {box().im(0), box().im(1), box().im(2), n_comps(), n_patches});
    storage().view(_all, _all, _all, _all, _all) = Real{};
  }

  KG_INLINE fields_view_t operator[](int p)
  {
    return fields_view_t(
      box().ib(), gt::gtensor_span<Real, 4, typename Storage::space_type>(
                    &storage()(0, 0, 0, 0, p),
                    gt::shape(box().im(0), box().im(1), box().im(2), n_comps()),
                    detail::strides(box().im(), n_comps())));
  }

  KG_INLINE const Real& operator()(int m, int i, int j, int k, int p) const
  {
    return storage()(i - ib(0), j - ib(1), k - ib(2), m, p);
  }

  KG_INLINE Real& operator()(int m, int i, int j, int k, int p)
  {
    return storage()(i - ib(0), j - ib(1), k - ib(2), m, p);
  }

  // protected:
  KG_INLINE Storage& storage() { return derived().storageImpl(); }
  KG_INLINE const Storage& storage() const { return derived().storageImpl(); }
  KG_INLINE Derived& derived() { return *static_cast<Derived*>(this); }
  KG_INLINE const Derived& derived() const
  {
    return *static_cast<const Derived*>(this);
  }

private:
  kg::Box3 box_; // size of one patch, including ghost points
  int n_fields_;
  int n_patches_;
};

// ======================================================================
// Mfields

template <typename R>
struct Mfields;

template <typename R>
struct MfieldsCRTPInnerTypes<Mfields<R>>
{
  using Storage = gt::gtensor<R, 5>;
};

template <typename R>
struct Mfields
  : MfieldsBase
  , MfieldsCRTP<Mfields<R>>
{
  using real_t = R;
  using Base = MfieldsCRTP<Mfields<R>>;
  using Storage = typename Base::Storage;
  using space = gt::space::host;

  Mfields(const Grid_t& grid, int n_fields, Int3 ibn)
    : MfieldsBase(grid, n_fields, ibn),
      Base(n_fields, {-ibn, grid.ldims + 2 * ibn}, grid.n_patches()),
      storage_(gt::shape(Base::box().im(0), Base::box().im(1),
                         Base::box().im(2), n_fields, Base::n_patches())),
      grid_{&grid}
  {
    std::fill(storage_.data(), storage_.data() + storage_.size(), real_t{});
  }

  Int3 ldims() const { return grid().ldims; }
  Int3 gdims() const { return grid().domain.gdims; }
  Int3 patchOffset(int p) const { return grid().patches[p].off; }
  const Grid_t& grid() const { return *grid_; }

  auto gt() { return Base::storage().view(); }

  template <typename FUNC>
  void Foreach_3d(int l, int r, FUNC&& F) const
  {
    return grid().Foreach_3d(l, r, std::forward<FUNC>(F));
  }

  virtual void reset(const Grid_t& grid) override
  {
    MfieldsBase::reset(grid);
    Base::reset(grid.n_patches());
    grid_ = &grid;
  }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override;
  const Convert& convert_from() override;

private:
  Storage storage_;
  const Grid_t* grid_;

  Storage& storageImpl() { return storage_; }
  const Storage& storageImpl() const { return storage_; }

  friend class MfieldsCRTP<Mfields<R>>;
};

template <>
const MfieldsBase::Convert Mfields<float>::convert_to_;
extern template const MfieldsBase::Convert Mfields<float>::convert_to_;

template <>
const MfieldsBase::Convert Mfields<float>::convert_from_;
extern template const MfieldsBase::Convert Mfields<float>::convert_from_;

template <>
const MfieldsBase::Convert Mfields<double>::convert_to_;
extern template const MfieldsBase::Convert Mfields<double>::convert_to_;

template <>
const MfieldsBase::Convert Mfields<double>::convert_from_;
extern template const MfieldsBase::Convert Mfields<double>::convert_from_;

template <typename R>
inline auto Mfields<R>::convert_to() -> const Convert&
{
  return convert_to_;
}

template <typename R>
inline auto Mfields<R>::convert_from() -> const Convert&
{
  return convert_from_;
}

// ======================================================================
// MfieldsStateFromMfields

template <typename Mfields>
struct MfieldsStateFromMfields : MfieldsStateBase
{
  using fields_view_t = typename Mfields::fields_view_t;
  using real_t = typename Mfields::real_t;
  using Real = typename Mfields::Real;
  using Storage = typename Mfields::Storage;
  using space = gt::space::host;

  MfieldsStateFromMfields(const Grid_t& grid)
    : MfieldsStateBase{grid, NR_FIELDS,
                       grid.ibn}, // FIXME, still hacky ibn handling...
      mflds_{grid, NR_FIELDS, grid.ibn}
  {}

  const Grid_t& grid() const { return mflds_.grid(); }
  int n_patches() const { return mflds_.n_patches(); }
  Int3 ib() const { return mflds_.ib(); }
  Int3 im() const { return mflds_.im(); }
  int n_comps() const { return mflds_.n_comps(); }

  fields_view_t operator[](int p) { return mflds_[p]; }

  KG_INLINE const real_t& operator()(int m, int i, int j, int k, int p) const
  {
    return mflds_(m, i, j, k, p);
  }

  KG_INLINE real_t& operator()(int m, int i, int j, int k, int p)
  {
    return mflds_(m, i, j, k, p);
  }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  auto& storage() const { return mflds_.storage(); }
  auto& storage() { return mflds_.storage(); }

  Mfields& mflds() { return mflds_; }

  auto gt() { return mflds_.storage(); }

public: // FIXME public so that we can read/write it, friend needs include which
        // gives nvcc issues
  Mfields mflds_;

  //  friend class kg::io::Descr<MfieldsStateFromMfields<Mfields>>;
};

template <>
const MfieldsStateBase::Convert
  MfieldsStateFromMfields<Mfields<float>>::convert_to_;
extern template const MfieldsStateBase::Convert
  MfieldsStateFromMfields<Mfields<float>>::convert_to_;

template <>
const MfieldsStateBase::Convert
  MfieldsStateFromMfields<Mfields<float>>::convert_from_;
extern template const MfieldsStateBase::Convert
  MfieldsStateFromMfields<Mfields<float>>::convert_from_;

template <>
const MfieldsStateBase::Convert
  MfieldsStateFromMfields<Mfields<double>>::convert_to_;
extern template const MfieldsStateBase::Convert
  MfieldsStateFromMfields<Mfields<double>>::convert_to_;

template <>
const MfieldsStateBase::Convert
  MfieldsStateFromMfields<Mfields<double>>::convert_from_;
extern template const MfieldsStateBase::Convert
  MfieldsStateFromMfields<Mfields<double>>::convert_from_;

using namespace gt::placeholders;

template <typename Mfields>
auto adapt(const Mfields& _mflds)
{
  auto& mflds = const_cast<std::remove_const_t<Mfields>&>(_mflds);
  auto ib = mflds.ib(), im = mflds.im(), bnd = -ib;
  return mflds.storage().view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]),
                              _s(bnd[2], -bnd[2]));
}

// ======================================================================

namespace gt
{

template <typename E>
auto host_mirror(const E& e)
{
  // FIXME, empty_like with space would be helpful
  return gt::empty<typename E::value_type>(e.shape());
}

} // namespace gt

template <typename MF>
using hostMirror_t =
  std::decay_t<decltype(hostMirror(std::declval<const MF>()))>;

template <typename MF>
MF& hostMirror(MF& mflds)
{
  return mflds;
}

// FIXME, doesn't actually copy, only for hostMirror use
template <typename MF>
void copy(const MF& from, MF& to)
{
  assert(from.storage().data() == to.storage().data());
}

// ======================================================================
// Mfields_from_gt_t

namespace detail
{
template <typename T, typename S>
struct Mfields_from_type_space
{
  using type = Mfields<T>;
};
} // namespace detail

template <typename E>
using Mfields_from_gt_t =
  typename detail::Mfields_from_type_space<typename E::value_type,
                                           typename E::space>::type;

// ======================================================================
// psc::interior

namespace psc
{
template <typename E>
GT_INLINE auto interior(E&& e, const Int3& ib)
{
  return std::forward<E>(e).view(_s(-ib[0], ib[0]), _s(-ib[1], ib[1]),
                                 _s(-ib[2], ib[2]));
}

namespace mflds
{

template <typename R, typename S = gt::space::host>
GT_INLINE auto gtensor(const Grid_t& grid, int n_comps, Int3 ibn = {})
{
  return gt::empty<R, S>(
    {grid.ldims[0] + 2 * ibn[0], grid.ldims[1] + 2 * ibn[1],
     grid.ldims[2] + 2 * ibn[2], n_comps, grid.n_patches()});
}

} // namespace mflds

} // namespace psc

#endif
