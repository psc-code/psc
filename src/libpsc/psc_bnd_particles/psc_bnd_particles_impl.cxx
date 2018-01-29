
#include "psc_bnd_particles_private.h"
#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "psc_particles_fortran.h"

//#include "psc_bnd_particles_open.cxx"

#include "bnd_particles_impl.hxx"
#include "bnd_particles_ordered_impl.hxx"

// ======================================================================
// psc_bnd_particles: subclass "single"

struct psc_bnd_particles_ops_single : psc_bnd_particles_ops {
  using sub_t = psc_bnd_particles_sub<mparticles_single_t>;
  psc_bnd_particles_ops_single() {
    name                    = "single";
    size                    = sizeof(sub_t);
    destroy                 = sub_t::destroy;
    setup                   = sub_t::setup;
    unsetup                 = sub_t::unsetup;
    exchange_particles      = sub_t::exchange_particles;
    //open_calc_moments       = psc_bnd_particles_sub_open_calc_moments;
  }
} psc_bnd_particles_single_ops;

// ======================================================================
// psc_bnd_particles: subclass "double"

struct psc_bnd_particles_ops_double : psc_bnd_particles_ops {
  using sub_t = psc_bnd_particles_sub<mparticles_double_t>;
  psc_bnd_particles_ops_double() {
    name                    = "double";
    size                    = sizeof(sub_t);
    destroy                 = sub_t::destroy;
    setup                   = sub_t::setup;
    unsetup                 = sub_t::unsetup;
    exchange_particles      = sub_t::exchange_particles;
    //open_calc_moments       = psc_bnd_particles_sub_open_calc_moments;
  }
} psc_bnd_particles_double_ops;

// ======================================================================
// psc_bnd_particles: subclass "fortran"

struct psc_bnd_particles_ops_fortran : psc_bnd_particles_ops {
  using sub_t = psc_bnd_particles_sub<mparticles_fortran_t>;
  psc_bnd_particles_ops_fortran() {
    name                    = "fortran";
    size                    = sizeof(sub_t);
    destroy                 = sub_t::destroy;
    setup                   = sub_t::setup;
    unsetup                 = sub_t::unsetup;
    exchange_particles      = sub_t::exchange_particles;
    //open_calc_moments       = psc_bnd_particles_sub_open_calc_moments;
  }
} psc_bnd_particles_fortran_ops;

// ======================================================================
// psc_bnd_particles: subclass "single2"

struct psc_bnd_particles_ops_single2 : psc_bnd_particles_ops {
  using sub_t = psc_bnd_particles_sub<mparticles_single_t,
				      bnd_particles_policy_ordered<mparticles_single_t>>;
  psc_bnd_particles_ops_single2() {
    name                    = "single2";
    size                    = sizeof(sub_t);
    destroy                 = sub_t::destroy;
    setup                   = sub_t::setup;
    unsetup                 = sub_t::unsetup;
    exchange_particles      = sub_t::exchange_particles;
  }
} psc_bnd_particles_single2_ops;
