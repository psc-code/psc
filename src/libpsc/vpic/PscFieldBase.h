
#ifndef PSC_FIELD_BASE_H
#define PSC_FIELD_BASE_H

template<class E, class G>
struct PscFieldBase
{
  typedef E Element;
  typedef G Grid;

  PscFieldBase(Grid* grid)
    : g_(grid), is_owner_(true)
  {
    arr_ = new Element[g_->nv]();
  }

  PscFieldBase(Grid* grid, Element* arr)
    : arr_(arr), g_(grid), is_owner_(false)
  {}

  PscFieldBase(const PscFieldBase& other)
    : arr_{other.arr_}, g_{other.g_}, is_owner_{other.is_owner_}
  {
    // FIXME, copying leaves us with two owners of the same data, which
    // can't be good
    assert(!is_owner_);
  }

  ~PscFieldBase()
  {
    if (is_owner_) {
      delete[] arr_;
    }
  }

  Element  operator[](int idx) const { return arr_[idx]; }
  Element& operator[](int idx)       { return arr_[idx]; }

  Element* data() { return arr_; }

  Grid* grid() { return g_; }

protected:
  Element* arr_;
  Grid* g_;

private:
  bool is_owner_; // OPT this could probably somehow be handled by a template...
};


#endif
