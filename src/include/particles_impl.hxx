
#pragma once

template<typename S>
template<typename MP>
inline MP PscMparticles<S>::get_as(uint flags)
{
  // If we're already the subtype, nothing to be done
  if (typeid(*sub()) == typeid(typename MP::sub_t)) {
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
  
  copy(*sub(), *mprts_to.sub(), type_from, type_to, flags);
  
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
    
    flags |= MP_DONT_RESIZE;
  }
  
  copy(mp_from, mp_to, type_from, type_to, flags);
  
  psc_mparticles_destroy(mprts_);
  mprts_ = nullptr;
  
  //  mprintf("put_as %s -> %s to\n", type_from, type_to);
  //  psc_mparticles_view(mprts_to);
  prof_stop(pr);
}

template<typename S>
inline void PscMparticles<S>::copy(MparticlesBase& mp_from, MparticlesBase& mp_to,
				   const char *type_from, const char *type_to,
				   unsigned int flags)
{
  // FIXME, could check for equal grid
  assert(mp_from.n_patches() == mp_to.n_patches());
  
  if (flags & MP_DONT_COPY) {
    if (!(flags & MP_DONT_RESIZE)) {
      uint n_prts_by_patch[mp_from.n_patches()];
      mp_from.get_size_all(n_prts_by_patch);
      mp_to.reserve_all(n_prts_by_patch);
      mp_to.resize_all(n_prts_by_patch);
    }
    return;
  }
  
  assert(!(flags & MP_DONT_RESIZE));
  
  auto it = mp_from.convert_to().find(std::string{"copy_to_"} + type_to);
  auto copy_to = it != mp_from.convert_to().cend() ? it->second : nullptr;
  if (copy_to) {
    copy_to(mp_from, mp_to);
  } else {
    auto it = mp_to.convert_from().find(std::string{"copy_from_"} + type_from);
    auto copy_from = it != mp_to.convert_from().cend() ? it->second : nullptr;
    if (copy_from) {
      copy_from(mp_to, mp_from);
    } else {
      fprintf(stderr, "ERROR: no 'copy_to_%s' in psc_mparticles '%s' and "
	      "no 'copy_from_%s' in '%s'!\n",
	      type_to, type_from, type_from, type_to);
      assert(0);
    }
  }
}

 

