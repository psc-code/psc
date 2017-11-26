
#ifndef MATERIAL_H
#define MATERIAL_H

#include "material/material.h"

// ======================================================================
// MaterialList

struct MaterialList {
  material_t* append(material_t* m);
  bool empty();
  
  //private:
  material_t* ml_;
};

inline material_t* MaterialList::append(material_t* m)
{
  return ::append_material(m, &ml_);
}

inline bool MaterialList::empty()
{
  return !ml_;
}


#endif
