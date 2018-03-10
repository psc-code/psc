
#include "psc_collision_impl.hxx"

#include "psc_fields_single.h"
#include "psc_particles_single.h"
#include "psc_fields_c.h"
#include "psc_particles_double.h"
#include "fields_item.hxx"

#include <string>

// ======================================================================
// psc_collision: subclass "single"/"double"

psc_collision_ops_<Collision_<MparticlesSingle, MfieldsSingle>> psc_collision_single_ops;
psc_collision_ops_<Collision_<MparticlesDouble, MfieldsC>> psc_collision_double_ops;

template<typename Collision>
struct FieldsItem_coll_stats : FieldsItemCRTP<FieldsItem_coll_stats<Collision>>
{
  using Base = FieldsItemCRTP<FieldsItem_coll_stats<Collision>>;
  using Base::Base;
  
  using mparticles_t = typename Collision::mparticles_t;
  using mfields_t = typename Collision::mfields_t;

  static const char* name()
  {
    return strdup((std::string("coll_stats_") +
		   mparticles_traits<mparticles_t>::name).c_str());
  }
  constexpr static int n_comps = Collision::NR_STATS;
  constexpr static fld_names_t fld_names()
  {
    return { "coll_nudt_min", "coll_nudt_med", "coll_nudt_max", "coll_nudt_large",
	     "coll_ncoll" };
  }
  constexpr static int flags = 0;

  void run(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) override
  {
    PscCollision<Collision> collision(ppsc->collision);
    Collision* coll = collision.sub();
    
    mfields_t mr = this->mres_base_->template get_as<mfields_t>(0, 0);
    
    for (int m = 0; m < coll->NR_STATS; m++) {
      // FIXME, copy could be avoided (?)
      mr->copy_comp(m, *mfields_t(coll->mflds).sub(), m);
    }
    
    mr.put_as(this->mres_base_, 0, coll->NR_STATS);
  }
};

template<typename Collision>
struct FieldsItem_coll_rei : FieldsItemCRTP<FieldsItem_coll_rei<Collision>>
{
  using Base = FieldsItemCRTP<FieldsItem_coll_rei<Collision>>;
  using Base::Base;
  
  using mparticles_t = typename Collision::mparticles_t;
  using mfields_t = typename Collision::mfields_t;

  static const char* name()
  {
    return strdup((std::string("coll_rei_") +
		   mparticles_traits<mparticles_t>::name).c_str());
  }
  constexpr static int n_comps = 3;
  constexpr static fld_names_t fld_names() { return { "coll_rei_x", "coll_rei_y", "coll_rei_z" }; }
  constexpr static int flags = 0;

  void run(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) override
  {
    PscCollision<Collision> collision(ppsc->collision);
    Collision* coll = collision.sub();
    
    mfields_t mr = this->mres_base_->template get_as<mfields_t>(0, 0);
    
    for (int m = 0; m < 3; m++) {
      // FIXME, copy could be avoided (?)
      mr->copy_comp(m, *mfields_t(coll->mflds_rei).sub(), m);
    }
    
    mr.put_as(this->mres_base_, 0, 3);
  }
};
  
// ======================================================================
// psc_output_fields_item: subclass "coll_stats" / "coll_rei"

using CollisionSingle = Collision_<MparticlesSingle, MfieldsSingle>;
FieldsItemOps<FieldsItem_coll_stats<CollisionSingle>> psc_output_fields_item_coll_stats_single_ops;
FieldsItemOps<FieldsItem_coll_rei<CollisionSingle>> psc_output_fields_item_coll_rei_single_ops;

using CollisionDouble = Collision_<MparticlesDouble, MfieldsC>;
FieldsItemOps<FieldsItem_coll_stats<CollisionDouble>> psc_output_fields_item_coll_stats_double_ops;
FieldsItemOps<FieldsItem_coll_rei<CollisionDouble>> psc_output_fields_item_coll_rei_double_ops;

