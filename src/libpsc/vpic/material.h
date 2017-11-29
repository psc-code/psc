
#ifndef MATERIAL_H
#define MATERIAL_H

#include "material/material.h"

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
  VpicMaterial(const char *name,
	       float epsx, float epsy, float epsz,
	       float mux, float muy, float muz,
	       float sigmax, float sigmay, float sigmaz,
	       float zetax, float zetay, float zetaz)
  {
    material_ctor(this, name, epsx, epsy, epsz, mux, muy, muz,
		  sigmax, sigmay, sigmaz, zetax, zetay, zetaz);
  }
  
  ~VpicMaterial()
  {
    material_dtor(this);
  }
};

// ======================================================================
// VpicMaterialList

struct VpicMaterialList
{
  typedef VpicMaterial Material;
  
  Material* append(Material* m)
  {
    return static_cast<Material*>(::append_material(m, &ml_));
  }

  bool empty()
  {
    return !ml_;
  }

  operator const material_t * () const
  {
    return ml_;
  }
  
private:
  material_t* ml_;
};

// ======================================================================
// PscMaterialList

struct PscMaterialList
{
  typedef VpicMaterial Material;
  
  Material* append(Material* m)
  {
    return static_cast<Material*>(::append_material(m, reinterpret_cast<material_t**>(&ml_)));
  }

  bool empty()
  {
    return !ml_;
  }

  operator const material_t * () const
  {
    return ml_;
  }
  
private:
  Material* ml_;
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
  MALLOC(m->name, len+1);
  strcpy(m->name, name);
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
