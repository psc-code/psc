
#ifndef PSC_FIELD_BASE_H
#define PSC_FIELD_BASE_H

template<class E, class G>
struct PscFieldBase
{
  typedef E Element;
  typedef G Grid;

  PscFieldBase(Grid* grid, Element* arr = 0)
    : arr_(arr),
      g(grid)
  {
  }

  ~PscFieldBase()
  {
  }

  Element *data() { return arr_; }

  Grid* getGrid() { return g; }

protected:
  Element* ALIGNED(128) arr_;

public:
  Grid* g;
};


#endif
