
#ifndef FIELDS3D_HXX
#define FIELDS3D_HXX

#include "psc.h"

#include "grid.hxx"
#include <mrc_io.hxx>
#include <kg/SArrayContainer.h>

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
// StorageNoOwnership

template <typename T>
class StorageNoOwnership
{
public:
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;

  StorageNoOwnership(pointer data) : data_{data} {}

  const_reference operator[](int offset) const { return data_[offset]; }
  reference operator[](int offset) { return data_[offset]; }

  // FIXME access to underlying storage might better be avoided?
  // use of this makes assumption that storage is contiguous
  const_pointer data() const { return data_; }
  pointer data() { return data_; }

private:
  pointer data_;
};

// ======================================================================
// fields3d_view

template<typename T, typename L=kg::LayoutSOA>
struct fields3d_view;

namespace kg
{

template<typename T, typename L>
struct SArrayContainerInnerTypes<fields3d_view<T, L>>
{
  using Layout = L;
  using Storage = StorageNoOwnership<T>;
};

}

template<typename T, typename L>
struct fields3d_view : kg::SArrayContainer<fields3d_view<T, L>>
{
  using Base = kg::SArrayContainer<fields3d_view<T, L>>;
  using Storage = typename Base::Storage;
  using real_t = typename Base::value_type;

  fields3d_view(Int3 ib, Int3 im, int n_comps, real_t* data)
    : Base{ib, im, n_comps},
      storage_{data}
  {
  }

private:
  Storage storage_;

  Storage& storageImpl() { return storage_; }
  const Storage& storageImpl() const { return storage_; }

  friend class kg::SArrayContainer<fields3d_view<T, L>>;
};

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

  virtual ~MfieldsBase()
  {
    instances.remove(this);
  }

  virtual void reset(const Grid_t& grid) { grid_ = &grid; }
  
  int n_patches() const { return grid_->n_patches(); }
  int n_comps() const { return n_fields_; }
  Int3 ibn() const { return ibn_; }

  virtual void write_as_mrc_fld(mrc_io *io, const std::string& name, const std::vector<std::string>& comp_names) = 0;

  const Grid_t& grid() const { return *grid_; }
  
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
    
    auto& mflds = *new MF{grid(), n_comps(), ibn()};
    
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
  
protected:
  int n_fields_;
  const Grid_t* grid_;
  Int3 ibn_;
};

#if 0

using MfieldsStateBase = MfieldsBase;

#else

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

  int n_patches() const { return grid_->n_patches(); }
  int n_comps() const { return n_fields_; }
  Int3 ibn() const { return ibn_; }

  virtual void write_as_mrc_fld(mrc_io *io, const std::string& name, const std::vector<std::string>& comp_names)
  {
    assert(0);
  }

  const Grid_t& grid() { return *grid_; }

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
    
    auto& mflds = *new MF{grid()};
    
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

#endif

// ======================================================================
// MfieldsCRTP

template <typename D, typename R>
class MfieldsCRTP
{
public:
  using Derived = D;

  using fields_view_t = fields3d_view<R>;

  MfieldsCRTP(int n_fields, Int3 ib, Int3 im, int n_patches)
    : n_fields_(n_fields), ib_{ib}, im_{im}, n_patches_{n_patches}
  {}
  
  fields_view_t operator[](int p)
  {
    return fields_view_t(Int3::fromPointer(ib_), Int3::fromPointer(im_), n_fields_, data_[p].get());
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
  Derived& derived() { return *static_cast<Derived*>(this); }
  const Derived& derived() const { return *static_cast<const Derived*>(this); }

protected:
  std::vector<std::unique_ptr<R[]>> data_;
  Int3 ib_; //> lower left corner for each patch (incl. ghostpoints)
  Int3 im_; //> extent for each patch (incl. ghostpoints)
  int n_fields_;
  int n_patches_;
};

// ======================================================================
// Mfields

template<typename R>
struct Mfields : MfieldsBase, MfieldsCRTP<Mfields<R>, R>
{
  using real_t = R;
  using Base = MfieldsCRTP<Mfields<R>, R>;

  using Base::data_;
  using Base::ib_;
  using Base::im_;

  Mfields(const Grid_t& grid, int n_fields, Int3 ibn)
    : MfieldsBase(grid, n_fields, ibn),
      Base(n_fields, -ibn, grid.ldims + 2 * ibn, grid.n_patches())
  {
    unsigned int size = 1;
    for (int d = 0; d < 3; d++) {
      size *= im_[d];
    }

    data_.reserve(n_patches());
    for (int p = 0; p < n_patches(); p++) {
      data_.emplace_back(new real_t[n_fields * size]{});
    }
  }

  virtual void reset(const Grid_t& grid) override
  {
    MfieldsBase::reset(grid);

    unsigned int size = 1;
    for (int d = 0; d < 3; d++) {
      size *= im_[d];
    }

    data_.clear();
    data_.reserve(n_patches());
    for (int p = 0; p < n_patches(); p++) {
      data_.emplace_back(new real_t[n_comps() * size]);
    }
  }
  
  void write_as_mrc_fld(mrc_io *io, const std::string& name, const std::vector<std::string>& comp_names) override
  {
    MrcIo::write_mflds(io, *this, name, comp_names);
  }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override;
  const Convert& convert_from() override;
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

  fields_view_t operator[](int p) { return mflds_[p]; }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  Mfields& mflds() { return mflds_; }

private:
  Mfields mflds_;
};

#endif
