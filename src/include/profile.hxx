
#include "mrc_profile.h"

// ======================================================================
// Profile
//
// This class behaves just as the underlying class T, but does perform
// profiling as the methods are called.
// The fact that is has wrappers for all potentially used methods isn't
// really pretty :(

template<typename T>
class Profile
{
public:
  template<typename... Args>
  Profile(Args&&... args)
    : obj_(std::forward<Args>(args)...)
  {
    if (!pr_) {
      pr_ = prof_register(typeid(T).name(), 1., 0, 0);
    }
  }

#define FWD_METHOD(name)				\
  template<typename ...Args>				\
  void name(Args&&... args)				\
  {							\
    prof_start(pr_);					\
    obj_.name(std::forward<Args>(args)...);		\
    prof_stop(pr_);					\
  }

  FWD_METHOD(operator())
  FWD_METHOD(push_mprts)
  FWD_METHOD(push_H)
  FWD_METHOD(push_E)
  FWD_METHOD(fill_ghosts)
  FWD_METHOD(add_ghosts)
  FWD_METHOD(fill_ghosts_E)
  FWD_METHOD(fill_ghosts_H)
  FWD_METHOD(add_ghosts_J)
  FWD_METHOD(gauss)
  FWD_METHOD(continuity_before_particle_push)
  FWD_METHOD(continuity_after_particle_push)

  template<typename ...Args>
  std::vector<uint> initial(Args&&... args)
  {
    prof_start(pr_);
    auto rv = obj_.initial(std::forward<Args>(args)...);
    prof_stop(pr_);
    return rv;
  }

private:
  T obj_;
  static int pr_;
};

template<typename T>
int Profile<T>::pr_;

