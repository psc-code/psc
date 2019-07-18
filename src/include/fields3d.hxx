
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

// FIXME, use size_t instead of int, at least for 1d offsets?

// ======================================================================
// fields3d

template<typename T>
class Storage
{
public:
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  
  Storage(pointer data)
    : data_{data} {}

  const_reference operator[](int offset) const { return data_[offset]; }
  reference operator[](int offset) { return data_[offset]; }

  // FIXME access to underlying storage might better be avoided?
  // use of this makes assumption that storage is contiguous
  const_pointer data() const { return data_; }
  pointer data() { return data_; }

  void free()
  {
    ::free(data_);
    data_ = nullptr;
  }
  
private:
  pointer data_;
};

// ======================================================================
// fields3d

template<typename T, typename L=kg::LayoutSOA>
struct fields3d;

template<typename T, typename L>
struct kg::SArrayContainerInnerTypes<fields3d<T, L>>
{
  using Layout = L;
  using Storage = Storage<T>;
};

template<typename T, typename L>
struct fields3d : kg::SArrayContainer<fields3d<T, L>>
{
  using Base = kg::SArrayContainer<fields3d<T, L>>;
  using Storage = typename Base::Storage;
  using real_t = typename Base::value_type;

  fields3d(const Grid_t& grid, Int3 ib, Int3 im, int n_comps)
    : Base{ib, im, n_comps},
      grid_{grid},
      storage_{(real_t *) calloc(Base::size(), sizeof(real_t))}
  {
  }

#if 0
  fields3d(const Grid_t& grid, Int3 ib, Int3 im, int n_comps, real_t* data)
    : Base{ib, im, n_comps},
      grid_{grid},
      storage_{data}
  {
  }
#endif
  
  void dtor()
  {
    storage_.free();
  }

  const Grid_t& grid() const { return grid_; }

private:
  Storage storage_;

  Storage& storageImpl() { return storage_; }
  const Storage& storageImpl() const { return storage_; }

  friend class kg::SArrayContainer<fields3d<T, L>>;
  
  const Grid_t& grid_;
};

// ======================================================================
// fields3d_view

template<typename T, typename L=kg::LayoutSOA>
struct fields3d_view;

template<typename T, typename L>
struct kg::SArrayContainerInnerTypes<fields3d_view<T, L>>
{
  using Layout = L;
  using Storage = Storage<T>;
};

template<typename T, typename L>
struct fields3d_view : kg::SArrayContainer<fields3d<T, L>>
{
  using Base = kg::SArrayContainer<fields3d<T, L>>;
  using Storage = typename Base::Storage;
  using real_t = typename Base::value_type;

  fields3d_view(const Grid_t& grid, Int3 ib, Int3 im, int n_comps, real_t* data)
    : Base{ib, im, n_comps},
      grid_{grid},
      storage_{data}
  {
  }

  const Grid_t& grid() const { return grid_; }

private:
  Storage storage_;

  Storage& storageImpl() { return storage_; }
  const Storage& storageImpl() const { return storage_; }

  friend class kg::SArrayContainer<fields3d<T, L>>;
  
  const Grid_t& grid_;
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

  virtual void zero_comp(int m) = 0;
  virtual void set_comp(int m, double val) = 0;
  virtual void scale_comp(int m, double val) = 0;
  virtual void axpy_comp(int m_y, double alpha, MfieldsBase& x, int m_x) = 0;
  virtual void copy_comp(int mto, MfieldsBase& from, int mfrom) = 0;
  virtual double max_comp(int m) = 0;
  virtual void write_as_mrc_fld(mrc_io *io, const std::string& name, const std::vector<std::string>& comp_names)
  {
    assert(0);
  }

  void zero()            { for (int m = 0; m < n_fields_; m++) zero_comp(m); }
  void scale(double val) { for (int m = 0; m < n_fields_; m++) scale_comp(m, val); }
  void axpy(double alpha, MfieldsBase& x)
  {
    for (int m = 0; m < n_fields_; m++) {
      axpy_comp(m, alpha, x, m);
    }
  }

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
// Mfields

template<typename R>
struct Mfields : MfieldsBase
{
  using real_t = R;
  using fields_view_t = fields3d_view<real_t>;

  Mfields(const Grid_t& grid, int n_fields, Int3 ibn)
    : MfieldsBase(grid, n_fields, ibn)
  {
    unsigned int size = 1;
    for (int d = 0; d < 3; d++) {
      ib[d] = -ibn[d];
      im[d] = grid_->ldims[d] + 2 * ibn[d];
      size *= im[d];
    }

    data.reserve(n_patches());
    for (int p = 0; p < n_patches(); p++) {
      data.emplace_back(new real_t[n_fields * size]{});
    }
  }

  virtual void reset(const Grid_t& grid) override
  {
    MfieldsBase::reset(grid);
    data.clear();

    unsigned int size = 1;
    for (int d = 0; d < 3; d++) {
      size *= im[d];
    }

    data.reserve(n_patches());
    for (int p = 0; p < n_patches(); p++) {
      data.emplace_back(new real_t[n_comps() * size]);
    }
  }
  
  fields_view_t operator[](int p)
  {
    return fields_view_t(grid(), Int3::fromPointer(ib), Int3::fromPointer(im), n_fields_, data[p].get());
  }

  void zero_comp(int m) override
  {
    for (int p = 0; p < n_patches(); p++) {
      (*this)[p].zero(m);
    }
  }

  void set_comp(int m, double val) override
  {
    for (int p = 0; p < n_patches(); p++) {
      (*this)[p].set(m, val);
    }
  }
  
  void scale_comp(int m, double val) override
  {
    for (int p = 0; p < n_patches(); p++) {
      (*this)[p].scale(m, val);
    }
  }

  void copy_comp(int mto, MfieldsBase& from_base, int mfrom) override
  {
    // FIXME? dynamic_cast would actually be more appropriate
    Mfields& from = static_cast<Mfields&>(from_base);
    for (int p = 0; p < n_patches(); p++) {
      (*this)[p].copy_comp(mto, from[p], mfrom);
    }
  }
  
  void axpy_comp(int m_y, double alpha, MfieldsBase& x_base, int m_x) override
  {
    // FIXME? dynamic_cast would actually be more appropriate
    Mfields& x = static_cast<Mfields&>(x_base);
    for (int p = 0; p < n_patches(); p++) {
      (*this)[p].axpy_comp(m_y, alpha, x[p], m_x);
    }
  }

  double max_comp(int m) override
  {
    double rv = -std::numeric_limits<double>::max();
    for (int p = 0; p < n_patches(); p++) {
      rv = std::max(rv, double((*this)[p].max_comp(m)));
    }
    return rv;
  }

  void write_as_mrc_fld(mrc_io *io, const std::string& name, const std::vector<std::string>& comp_names) override
  {
    MrcIo::write_mflds(io, *this, name, comp_names);
  }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  std::vector<std::unique_ptr<real_t[]>> data;
  int ib[3]; //> lower left corner for each patch (incl. ghostpoints)
  int im[3]; //> extent for each patch (incl. ghostpoints)
};

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
