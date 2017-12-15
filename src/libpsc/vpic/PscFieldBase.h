
#ifndef PSC_FIELD_BASE_H
#define PSC_FIELD_BASE_H

template<class E, class G>
struct PscFieldBase
{
  typedef E Element;
  typedef G Grid;

  PscFieldBase(Grid* grid) : g(grid)
  {
  }

  ~PscFieldBase()
  {
  }

  Grid* getGrid() { return g; }

  Grid* g;
};


#endif
