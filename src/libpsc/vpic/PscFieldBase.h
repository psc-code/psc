
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
    MALLOC_ALIGNED(arr_, grid->nv, 128);
    CLEAR(arr_, grid->nv);
  }

  PscFieldBase(Grid* grid, Element* arr)
    : arr_(arr), g_(grid), is_owner_(false)
  {
  }

  ~PscFieldBase()
  {
    if (is_owner_) {
      FREE_ALIGNED(arr_);
    }
  }

  Element  operator[](int idx) const { return arr_[idx]; }
  Element& operator[](int idx)       { return arr_[idx]; }

  Element *data() { return arr_; }

  Grid* getGrid() { return g_; }

protected:
  Element* ALIGNED(128) arr_;
  Grid* g_;

private:
  bool is_owner_; // OPT this could probably somehow be handled by a template...
};


#endif
