
#ifndef PSC_MATERIAL_H
#define PSC_MATERIAL_H

#include "psc_vpic_bits.h"
#include "VpicListBase.h"

#include <mrc_common.h>
#include <cassert>

// ======================================================================
// PscMaterial

struct PscMaterial
{
  PscMaterial(const char *name,
	      float epsx, float epsy, float epsz,
	      float mux, float muy, float muz,
	      float sigmax, float sigmay, float sigmaz,
	      float zetax, float zetay, float zetaz)
  {
    int len = name ? strlen(name) : 0;
    if (!len) LOG_ERROR("Cannot create a nameless material");
    this->name = strdup(name);
    this->epsx   = epsx,   this->epsy   = epsy,   this->epsz   = epsz;
    this->mux    = mux,    this->muy    = muy,    this->muz    = muz;
    this->sigmax = sigmax, this->sigmay = sigmay, this->sigmaz = sigmaz;
    this->zetax  = zetax,  this->zetay  = zetay,  this->zetaz  = zetaz;
    this->next = nullptr;
  }
  
  ~PscMaterial()
  {
    free(name);
  }

  char* name;                   // Name of the material
  float epsx, epsy, epsz;       // Relative permittivity along x,y,z axes
  float mux, muy, muz;          // Relative permeability along x,y,z axes
  float sigmax, sigmay, sigmaz; // Electrical conductivity along x,y,z axes
  float zetax,  zetay,  zetaz;  // Magnetic conductivity along x,y,z axes
  MaterialId id;                // Unique identifier for material
  PscMaterial *next;            // Next material in list
};

// ======================================================================
// PscMaterialList

struct PscMaterialList : public VpicListBase<PscMaterial>
{
  typedef PscMaterial Material;

  static Material* create(const char *name,
			  float epsx, float epsy, float epsz,
			  float mux, float muy, float muz,
			  float sigmax, float sigmay, float sigmaz,
			  float zetax, float zetay, float zetaz)
  {
    return new Material(name, epsx, epsy, epsz, mux, muy, muz,
    			sigmax, sigmay, sigmaz, zetax, zetay, zetaz);
  }

  const_iterator find(const char *name) const
  {
    assert(name);
    return std::find_if(cbegin(), cend(),
			[&name](const Material &m) { return strcmp(m.name, name) == 0; });
  }
  
  Material* append(Material* m)
  {
    assert(!m->next);
    if (find(m->name) != cend()) {
      LOG_ERROR("There is already a material named \"%s\" in list", m->name);
    }
    int id = size();
    if (id >= MaterialIdMax) {
      LOG_ERROR("Too many materials in list to append material \"%s\"", m->name);
    }
    m->id = id;
    push_front(*m);
    return m;
  }
};

#endif
