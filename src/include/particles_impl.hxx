
#pragma once

template<typename S>
template<typename MP>
inline MP PscMparticles<S>::get_as(uint flags)
{
  auto& mp_from = *sub();
  // If we're already the subtype, nothing to be done
  if (typeid(mp_from) == typeid(typename MP::sub_t)) {
    return MP{mprts_};
  }
  
  const char *type_to = mparticles_traits<MP>::name;
  const char *type_from = psc_mparticles_type(mprts_);

  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_get_as", 1., 0, 0);
  }
  
  prof_start(pr);
  
  //  mprintf("get_as %s -> %s from\n", type_from, type_to);
  //  psc_mparticles_view(mprts_);
  
  struct psc_mparticles *mprts_to_ = psc_mparticles_create(psc_mparticles_comm(mprts_));
  psc_mparticles_set_type(mprts_to_, type_to);
  mprts_to_->grid = mprts_->grid;
  psc_mparticles_setup(mprts_to_);
  auto mprts_to = MP{mprts_to_};
  auto& mp_to = *mprts_to.sub();
  
  if (flags & MP_DONT_COPY) {
    if (!(flags & MP_DONT_RESIZE)) {
      uint n_prts_by_patch[mp_from.n_patches()];
      mp_from.get_size_all(n_prts_by_patch);
      mp_to.reserve_all(n_prts_by_patch);
      mp_to.resize_all(n_prts_by_patch);
    }
  } else {
    assert(!(flags & MP_DONT_RESIZE));
    copy(mp_from, mp_to, flags);
  }
  
  //  mprintf("get_as %s -> %s to\n", type_from, type_to);
  //  psc_mparticles_view(mprts);
  
  prof_stop(pr);
  return mprts_to;
}

template<typename S>
template<typename MP>
inline void PscMparticles<S>::put_as(MP mprts_base, uint flags)
{
  auto& mp_from = *sub();
  auto& mp_to = *mprts_base.sub();
  
  // If we're already the subtype, nothing to be done
  if (typeid(mp_from) == typeid(mp_to)) {
    return;
  }
  
  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_put_as", 1., 0, 0);
  }
  prof_start(pr);
  
  const char *type_from = psc_mparticles_type(mprts_);
  const char *type_to = psc_mparticles_type(mprts_base.mprts());
  
  //  mprintf("put_as %s -> %s from\n", type_from, type_to);
  //  psc_mparticles_view(mprts);
  
  if (flags & MP_DONT_COPY) {
    // let's check that the size of the particle arrays hasn't changed, since
    // it's not obvious what we should do in case it did...
    assert(mp_to.n_patches() == mp_from.n_patches());
    
    uint n_prts_by_patch[mp_from.n_patches()];
    uint n_prts_by_patch_to[mp_to.n_patches()];
    mp_from.get_size_all(n_prts_by_patch);
    mp_to.get_size_all(n_prts_by_patch_to);
    
    for (int p = 0; p < mp_from.n_patches(); p++) {
      assert(n_prts_by_patch[p] == n_prts_by_patch_to[p]);
    }
  } else {
    copy(mp_from, mp_to, flags);
  }
  
  psc_mparticles_destroy(mprts_);
  mprts_ = nullptr;
  
  //  mprintf("put_as %s -> %s to\n", type_from, type_to);
  //  psc_mparticles_view(mprts_to);
  prof_stop(pr);
}

template<typename S>
inline void PscMparticles<S>::copy(MparticlesBase& mp_from, MparticlesBase& mp_to,
				   unsigned int flags)
{
  // FIXME, implementing == wouldn't hurt
  assert(&mp_from.grid() == &mp_to.grid());
  
  auto convert_to = mp_from.convert_to().find(std::type_index(typeid(mp_to)));
  if (convert_to != mp_from.convert_to().cend()) {
    convert_to->second(mp_from, mp_to);
    return;
  }
  
  auto convert_from = mp_to.convert_from().find(std::type_index(typeid(mp_from)));
  if (convert_from != mp_to.convert_from().cend()) {
    convert_from->second(mp_to, mp_from);
    return;
  }

  fprintf(stderr, "ERROR: no conversion known from %s to %s!\n",
	  typeid(mp_from).name(), typeid(mp_to).name());
  assert(0);
}

 

