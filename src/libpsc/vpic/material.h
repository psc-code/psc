
#ifndef MATERIAL_H
#define MATERIAL_H

#include "psc_vpic_bits.h"
#include "VpicListBase.h"
#include "material/material.h"

#include <mrc_common.h>

#include <cassert>

static void material_ctor(material_t *m, const char *name,
			  float epsx, float epsy, float epsz,
			  float mux, float muy, float muz,
			  float sigmax, float sigmay, float sigmaz,
			  float zetax, float zetay, float zetaz);
static void material_dtor(material_t *m);

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

// ======================================================================
// PscMaterial

struct PscMaterial : material_t
{
  PscMaterial(const char *name,
	      float epsx, float epsy, float epsz,
	      float mux, float muy, float muz,
	       float sigmax, float sigmay, float sigmaz,
	      float zetax, float zetay, float zetaz)
  {
    material_ctor(this, name, epsx, epsy, epsz, mux, muy, muz,
		  sigmax, sigmay, sigmaz, zetax, zetay, zetaz);
  }
  
  ~PscMaterial()
  {
    material_dtor(this);
  }
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
      ERROR(("There is already a material named \"%s\" in list", m->name ));
    }
    int id = size();
    if (id >= ::max_material) {
      ERROR(("Too many materials in list to append material \"%s\"", m->name));
    }
    m->id   = (material_id)id;
    push_front(*m);
    return m;
  }
};

// ======================================================================

static void material_ctor(material_t *m, const char *name,
	       float epsx, float epsy, float epsz,
	       float mux, float muy, float muz,
	       float sigmax, float sigmay, float sigmaz,
	       float zetax, float zetay, float zetaz)
{
  // copy&paste from material.c, but as ctor not new
  CLEAR(m, 1);

  int len = name ? strlen(name) : 0;
  if (!len) ERROR(( "Cannot create a nameless material" ));
  m->name = strdup(name);
  m->epsx   = epsx,   m->epsy   = epsy,   m->epsz   = epsz;
  m->mux    = mux,    m->muy    = muy,    m->muz    = muz;
  m->sigmax = sigmax, m->sigmay = sigmay, m->sigmaz = sigmaz;
  m->zetax  = zetax,  m->zetay  = zetay,  m->zetaz  = zetaz;
}

static void material_dtor(material_t *m)
{
  FREE(m->name);
}

#endif
