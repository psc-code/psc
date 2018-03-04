
#include "psc_collision_impl.hxx"

#include "psc_fields_single.h"
#include "psc_particles_single.h"
#include "psc_fields_c.h"
#include "psc_particles_double.h"

// ======================================================================
// psc_collision: subclass "single"/"double"

psc_collision_ops_<Collision_<PscMparticlesSingle, PscMfieldsSingle>> psc_collision_single_ops;
psc_collision_ops_<Collision_<PscMparticlesDouble, PscMfieldsC>> psc_collision_double_ops;


template<typename Collision>
static void copy_stats(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
		       struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  using mparticles_t = typename Collision::mparticles_t;
  using mfields_t = typename Collision::mfields_t;
  PscCollision<Collision> collision(ppsc->collision);
  Collision* coll = collision.sub();
  
  mfields_t mr = mres->get_as<mfields_t>(0, 0);
  
  for (int m = 0; m < coll->NR_STATS; m++) {
    // FIXME, copy could be avoided (?)
    mr->copy_comp(m, *mfields_t(coll->mflds).sub(), m);
  }
  
  mr.put_as(mres, 0, coll->NR_STATS);
}

template<typename Collision>
static void copy_rei(struct psc_output_fields_item *item, struct psc_mfields *mflds_base,
		     struct psc_mparticles *mprts_base, struct psc_mfields *mres)
{
  using mparticles_t = typename Collision::mparticles_t;
  using mfields_t = typename Collision::mfields_t;
  PscCollision<Collision> collision(ppsc->collision);
  Collision* coll = collision.sub();
  
  mfields_t mr = mres->get_as<mfields_t>(0, 0);
  
  for (int m = 0; m < 3; m++) {
    // FIXME, copy could be avoided (?)
    mr->copy_comp(m, *mfields_t(coll->mflds_rei).sub(), m);
  }
  
  mr.put_as(mres, 0, 3);
}
  
// ======================================================================
// psc_output_fields_item: subclass "coll_stats"

struct psc_output_fields_item_ops_coll_single : psc_output_fields_item_ops {
  using Collision = Collision_<PscMparticlesSingle, PscMfieldsSingle>;
  psc_output_fields_item_ops_coll_single() {
    name      = "coll_stats_single";
    nr_comp   = Collision::NR_STATS;
    fld_names[0] = "coll_nudt_min";
    fld_names[1] = "coll_nudt_med";
    fld_names[2] = "coll_nudt_max";
    fld_names[3] = "coll_nudt_large";
    fld_names[4] = "coll_ncoll";
    run_all   = copy_stats<Collision>;
  }
} psc_output_fields_item_coll_stats_single_ops;

struct psc_output_fields_item_ops_coll_double : psc_output_fields_item_ops {
  using Collision = Collision_<PscMparticlesDouble, PscMfieldsC>;
  psc_output_fields_item_ops_coll_double() {
    name      = "coll_stats_double";
    nr_comp   = Collision::NR_STATS;
    fld_names[0] = "coll_nudt_min";
    fld_names[1] = "coll_nudt_med";
    fld_names[2] = "coll_nudt_max";
    fld_names[3] = "coll_nudt_large";
    fld_names[4] = "coll_ncoll";
    run_all   = copy_stats<Collision>;
  }
} psc_output_fields_item_coll_stats_double_ops;

// ======================================================================
// psc_output_fields_item: subclass "coll_rei"

struct psc_output_fields_item_ops_coll_rei_single : psc_output_fields_item_ops {
  using Collision = Collision_<PscMparticlesSingle, PscMfieldsSingle>;
  psc_output_fields_item_ops_coll_rei_single() {
    name      = "coll_rei_single";
    nr_comp   = 3;
    fld_names[0] = "coll_rei_x";
    fld_names[1] = "coll_rei_y";
    fld_names[2] = "coll_rei_z";
    run_all   = copy_rei<Collision>;
  }
} psc_output_fields_item_coll_rei_single_ops;

struct psc_output_fields_item_ops_coll_rei_double : psc_output_fields_item_ops {
  using Collision = Collision_<PscMparticlesDouble, PscMfieldsC>;
  psc_output_fields_item_ops_coll_rei_double() {
    name      = "coll_rei_double";
    nr_comp   = 3;
    fld_names[0] = "coll_rei_x";
    fld_names[1] = "coll_rei_y";
    fld_names[2] = "coll_rei_z";
    run_all   = copy_rei<Collision>;
  }
} psc_output_fields_item_coll_rei_double_ops;
