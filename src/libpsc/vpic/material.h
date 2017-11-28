
#ifndef MATERIAL_H
#define MATERIAL_H

#include "material/material.h"

// ======================================================================
// VpicMaterialList

struct VpicMaterialList {
  material_t* append(material_t* m)
  {
    return ::append_material(m, &ml_);
  }

  bool empty()
  {
    return !ml_;
  }
  
  //private:
  material_t* ml_;
};


#endif
