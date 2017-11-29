
#ifndef MATERIAL_H
#define MATERIAL_H

#include "material/material.h"

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
// VpicMaterialConstIter

struct VpicMaterialConstIter {
  typedef VpicMaterialConstIter ConstIter;
  
  VpicMaterialConstIter(const VpicMaterial *node=nullptr) : node_(node)
  {
  }

  bool operator!=(const ConstIter& x) const
  {
    return node_ != x.node_;
  }

  ConstIter& operator++()
  {
    node_ = static_cast<VpicMaterial *>(node_->next);
    return *this;
  }

  const VpicMaterial& operator*() const
  {
    return *node_;
  }
  
  const VpicMaterial *operator->() const
  {
    return node_;
  }
  
private:
  const VpicMaterial *node_;
};

// ======================================================================
// VpicMaterialList

struct VpicMaterialList
{
  typedef VpicMaterial Material;
  typedef VpicMaterialConstIter ConstIter;
  
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
  typedef VpicMaterialConstIter ConstIter;

  ConstIter cbegin() const
  {
    return VpicMaterialConstIter(ml_);
  }
  
  ConstIter cend() const
  {
    return VpicMaterialConstIter(nullptr);
  }

  size_t size() const
  {
    size_t sz = 0;
    for (ConstIter m = cbegin(); m != cend(); ++m) {
      sz++;
    }
    return sz;
  }
  
  ConstIter find(const char *name) const
  {
    assert(name);
    for (ConstIter m = cbegin(); m != cend(); ++m) {
      if (strcmp(name, m->name) == 0) {
	return m;
      }
    }
    return cend();
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
    m->next = ml_;
    ml_ = m;
    return m;
  }

  bool empty()
  {
    return !ml_;
  }

  // operator const material_t * () const
  // {
  //   return ml_;
  // }
  
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
