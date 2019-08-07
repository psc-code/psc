
#ifndef FIELDS3D_HXX
#define FIELDS3D_HXX

#include "psc.h"

#include "grid.hxx"
#include <mrc_io.hxx>
#include <kg/SArrayView.h>

#include <mrc_profile.h>

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
// MfieldsBase

struct MfieldsBase
{
  using convert_func_t = void (*)(MfieldsBase&, MfieldsBase&, int, int);
  using Convert = std::unordered_map<std::type_index, convert_func_t>;
  
  MfieldsBase(const Grid_t& grid, int n_fields, Int3 ibn)
    : grid_(&grid),
      n_fields_(n_fields),
      ibn_(ibn)
  {
    instances.push_back(this);
  }

  MfieldsBase(const MfieldsBase&) = delete;
  MfieldsBase& operator=(const MfieldsBase&) = delete;

  MfieldsBase(MfieldsBase&& o) : MfieldsBase(*o.grid_, o.n_fields_, o.ibn_)
  {
  }

  MfieldsBase& operator=(MfieldsBase&& o) = default;
  
  virtual ~MfieldsBase()
  {
    instances.remove(this);
  }

  virtual void reset(const Grid_t& grid) { grid_ = &grid; }
  
  int _n_comps() const { return n_fields_; }
  Int3 ibn() const { return ibn_; }

  virtual void write_as_mrc_fld(mrc_io *io, const std::string& name, const std::vector<std::string>& comp_names) = 0;

  const Grid_t& _grid() const { return *grid_; }
  
  template<typename MF>
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

    // mprintf("get_as %s (%s) %d %d\n", type, psc_mfields_type(mflds_base), mb, me);
    
    auto& mflds = *new MF{_grid(), n_fields_, ibn()};
    
    MfieldsBase::convert(*this, mflds, mb, me);

    prof_stop(pr);
    return mflds;
  }

  template<typename MF>
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

  virtual const Convert& convert_to() { static const Convert convert_to_; return convert_to_; }
  virtual const Convert& convert_from() { static const Convert convert_from_; return convert_from_; }
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
  using convert_func_t = void (*)(MfieldsStateBase&, MfieldsStateBase&, int, int);
  using Convert = std::unordered_map<std::type_index, convert_func_t>;
  
  MfieldsStateBase(const Grid_t& grid, int n_fields, Int3 ibn)
    : grid_(&grid),
      n_fields_(n_fields),
      ibn_(ibn)
  {
    instances.push_back(this);
  }

  virtual ~MfieldsStateBase()
  {
    instances.remove(this);
  }

  virtual void reset(const Grid_t& grid) { grid_ = &grid; }

  int _n_patches() const { return grid_->n_patches(); }
  int _n_comps() const { return n_fields_; }
  Int3 ibn() const { return ibn_; }

  virtual void write_as_mrc_fld(mrc_io *io, const std::string& name, const std::vector<std::string>& comp_names)
  {
    assert(0);
  }

  const Grid_t& _grid() { return *grid_; }

  virtual const Convert& convert_to() { static const Convert convert_to_; return convert_to_; }
  virtual const Convert& convert_from() { static const Convert convert_from_; return convert_from_; }
  static void convert(MfieldsStateBase& mf_from, MfieldsStateBase& mf_to, int mb, int me);

  static std::list<MfieldsStateBase*> instances;

  template<typename MF>
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

    // mprintf("get_as %s (%s) %d %d\n", type, psc_mfields_type(mflds_base), mb, me);
    
    auto& mflds = *new MF{_grid()};
    
    MfieldsStateBase::convert(*this, mflds, mb, me);

    prof_stop(pr);
    return mflds;
  }

  template<typename MF>
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
  using fields_view_t = kg::SArrayView<Real>;

  KG_INLINE const kg::Box3& box() const { return box_; }
  KG_INLINE const Int3& ib() const { return box_.ib(); }
  KG_INLINE const Int3& im() const { return box_.im(); }
  KG_INLINE int ib(int d) const { return box_.ib(d); }
  KG_INLINE int im(int d) const { return box_.im(d); }
  KG_INLINE int n_comps() const { return n_fields_; }
  KG_INLINE int n_patches() const { return n_patches_; }

  typename Storage::iterator begin() { return storage().begin(); }
  typename Storage::iterator end() { return storage().end(); }
  typename Storage::const_iterator begin() const { return storage().begin(); }
  typename Storage::const_iterator end() const { return storage().end(); }

  MfieldsCRTP(int n_fields, const kg::Box3& box, int n_patches)
    : n_fields_(n_fields), box_{box}, n_patches_{n_patches}
  {}

  void reset(int n_patches)
  {
    n_patches_ = n_patches;
    storage().resize(n_patches * n_comps() * box().size());
  }
  
  KG_INLINE fields_view_t operator[](int p)
  {
    size_t stride = n_comps() * box().size();
    return fields_view_t(box(), n_comps(), &storage()[p * stride]);
  }

  void zero_comp(int m)
  {
    for (int p = 0; p < n_patches_; p++) {
      (*this)[p].zero(m);
    }
  }

  void zero()
  {
    for (int m = 0; m < n_fields_; m++) zero_comp(m);
  }
  
  void set_comp(int m, double val)
  {
    for (int p = 0; p < n_patches_; p++) {
      (*this)[p].set(m, val);
    }
  }
  
  void scale_comp(int m, double val)
  {
    for (int p = 0; p < n_patches_; p++) {
      (*this)[p].scale(m, val);
    }
  }

  void scale(double val)
  {
    for (int m = 0; m < n_fields_; m++) scale_comp(m, val);
  }
  
  void copy_comp(int mto, MfieldsBase& from_base, int mfrom)
  {
    // FIXME? dynamic_cast would actually be more appropriate
    Derived& from = static_cast<Derived&>(from_base);
    for (int p = 0; p < n_patches_; p++) {
      (*this)[p].copy_comp(mto, from[p], mfrom);
    }
  }
  
  void axpy_comp(int m_y, double alpha, MfieldsBase& x_base, int m_x)
  {
    // FIXME? dynamic_cast would actually be more appropriate
    Derived& x = static_cast<Derived&>(x_base);
    for (int p = 0; p < n_patches_; p++) {
      (*this)[p].axpy_comp(m_y, alpha, x[p], m_x);
    }
  }

  void axpy(double alpha, MfieldsBase& x)
  {
    for (int m = 0; m < n_fields_; m++) {
      axpy_comp(m, alpha, x, m);
    }
  }

  double max_comp(int m)
  {
    double rv = -std::numeric_limits<double>::max();
    for (int p = 0; p < n_patches_; p++) {
      rv = std::max(rv, double((*this)[p].max_comp(m)));
    }
    return rv;
  }

protected:
  KG_INLINE Storage& storage() { return derived().storageImpl(); }
  KG_INLINE const Storage& storage() const { return derived().storageImpl(); }
  KG_INLINE Derived& derived() { return *static_cast<Derived*>(this); }
  KG_INLINE const Derived& derived() const { return *static_cast<const Derived*>(this); }

private:
  kg::Box3 box_; // size of one patch, including ghost points
  int n_fields_;
  int n_patches_;
};

// ======================================================================
// MfieldsDomain

class MfieldsDomain
{
public:
  MfieldsDomain(const Grid_t& grid)
    : grid_{&grid},
      n_patches_{grid.n_patches()},
      ldims_{grid.ldims},
      gdims_{grid.domain.gdims}
  {
    patch_offsets_.reserve(n_patches());
    for (auto& patch : grid.patches) {
      patch_offsets_.emplace_back(patch.off);
    }
  }

  Int3 ldims() const { return ldims_; }
  Int3 gdims() const { return gdims_; }
  int n_patches() const { return n_patches_; }
  Int3 patchOffset(int p) const { return patch_offsets_[p]; }

  template<typename FUNC>
  void Foreach_3d(int l, int r, FUNC&& F) const
  {
    return grid().Foreach_3d(l, r, std::forward<FUNC>(F));
  }
  
  const Grid_t& grid() const { return *grid_; }
  
private:
  const Grid_t* grid_;
  Int3 ldims_;
  Int3 gdims_;
  int n_patches_;
  std::vector<Int3> patch_offsets_;
};

// ======================================================================
// Mfields

template<typename R>
struct Mfields;

template <typename R>
struct MfieldsCRTPInnerTypes<Mfields<R>>
{
  using Storage = std::vector<R>;
};

template<typename R>
struct Mfields : MfieldsBase, MfieldsCRTP<Mfields<R>>
{
  using real_t = R;
  using Base = MfieldsCRTP<Mfields<R>>;
  using Storage = typename Base::Storage;

  Mfields(const MfieldsDomain& domain, int n_fields, Int3 ibn)
    : MfieldsBase(domain.grid(), n_fields, ibn),
      Base(n_fields, {-ibn, domain.ldims() + 2 * ibn}, domain.n_patches()),
      storage_(size_t(Base::box().size() * n_fields * Base::n_patches())),
      domain_{domain}
  {}

  Int3 ldims() const { return domain_.ldims(); }
  Int3 gdims() const { return domain_.gdims(); }
  Int3 patchOffset(int p) const { return domain_.patchOffset(p); }
  const MfieldsDomain& domain() const { return domain_; }
  const Grid_t& grid() const { return domain_.grid(); }

  template<typename FUNC>
  void Foreach_3d(int l, int r, FUNC&& F) const
  {
    return domain().Foreach_3d(l, r, std::forward<FUNC>(F));
  }
  
  virtual void reset(const Grid_t& grid) override
  {
    MfieldsBase::reset(grid);
    Base::reset(grid.n_patches());
    domain_ = MfieldsDomain(grid);
  }
  
  void write_as_mrc_fld(mrc_io *io, const std::string& name, const std::vector<std::string>& comp_names) override
  {
    MrcIo::write_mflds(io, *this, grid(), name, comp_names);
  }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override;
  const Convert& convert_from() override;

private:
  Storage storage_;
  MfieldsDomain domain_;

  Storage& storageImpl() { return storage_; }
  const Storage& storageImpl() const { return storage_; }

  friend class MfieldsCRTP<Mfields<R>>;
};

template<> const MfieldsBase::Convert Mfields<float>::convert_to_;
extern template const MfieldsBase::Convert Mfields<float>::convert_to_;

template<> const MfieldsBase::Convert Mfields<float>::convert_from_;
extern template const MfieldsBase::Convert Mfields<float>::convert_from_;

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

template<typename Mfields>
struct MfieldsStateFromMfields : MfieldsStateBase
{
  using fields_view_t = typename Mfields::fields_view_t;
  using real_t = typename Mfields::real_t;

  MfieldsStateFromMfields(const Grid_t& grid)
    : MfieldsStateBase{grid, NR_FIELDS, grid.ibn}, // FIXME, still hacky ibn handling...
      mflds_{grid, NR_FIELDS, grid.ibn}
  {}

  const Grid_t& grid() const { return mflds_.grid(); }
  int n_patches() const { return mflds_.n_patches(); }
  
  fields_view_t operator[](int p) { return mflds_[p]; }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  Mfields& mflds() { return mflds_; }

public: // FIXME public so that we can read/write it, friend needs include which gives nvcc issues
  Mfields mflds_;
  
  //  friend class kg::io::Descr<MfieldsStateFromMfields<Mfields>>;
};

#endif
