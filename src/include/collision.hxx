
#pragma once

// ======================================================================
// PscCollision

template<typename S>
struct PscCollision
{
  using sub_t = S;
  
  explicit PscCollision(psc_collision *collision)
    : collision_(collision)
  {}
  
  sub_t* sub() { return mrc_to_subobj(collision_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_collision *collision_;
};

// ======================================================================
// CollisionWrapper

template<typename Collision>
class CollisionWrapper
{
public:
  const static size_t size = sizeof(Collision);

  constexpr static char const* const name = Collision::name;
  
  static void setup(struct psc_collision* _collision)
  {
    PscCollision<Collision> collision(_collision);
    
    new(collision.sub()) Collision(psc_collision_comm(_collision),
				   _collision->every, _collision->nu);
  }

  static void destroy(struct psc_collision* _collision)
  {
    PscCollision<Collision> collision(_collision);
    
    collision->~Collision();
  }

  static void run(struct psc_collision* _collision, struct psc_mparticles* mprts)
  {
    PscCollision<Collision> collision(_collision);

    collision->run(mprts);
  }

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

};

template<typename Collision>
struct psc_collision_ops_ : psc_collision_ops {
  using CollisionWrapper = CollisionWrapper<Collision>;
  psc_collision_ops_() {
    name                  = CollisionWrapper::name;
    size                  = CollisionWrapper::size;
    setup                 = CollisionWrapper::setup;
    destroy               = CollisionWrapper::destroy;
    run                   = CollisionWrapper::run;
  }
};

