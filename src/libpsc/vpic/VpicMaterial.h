
#ifndef VPIC_MATERIAL_H
#define VPIC_MATERIAL_H

#include "VpicListBase.h"

// ======================================================================
// VpicMaterial

struct VpicMaterial : material_t
{
};

// ======================================================================
// VpicMaterialList

struct VpicMaterialList : public VpicListBase<VpicMaterial>
{
  typedef VpicMaterial Material;
  
  static Material* create(const char *name,
			  float epsx, float epsy, float epsz,
			  float mux, float muy, float muz,
			  float sigmax, float sigmay, float sigmaz,
			  float zetax, float zetay, float zetaz)
  {
    return static_cast<Material*>(::material(name, epsx, epsy, epsz, mux, muy, muz,
					     sigmax, sigmay, sigmaz, zetax, zetay, zetaz));
  }

  Material* append(Material* m)
  {
    return static_cast<Material*>(::append_material(m, reinterpret_cast<material_t**>(&head_)));
  }

  operator const material_t * () const
  {
    return head_;
  }
};



#endif
