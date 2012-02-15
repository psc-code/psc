
#include "psc_moments_private.h"

#include "psc_bnd.h"

// ======================================================================
// forward to subclass

void
psc_moments_calc_densities(struct psc_moments *moments,
			   mfields_base_t *flds, mparticles_base_t *particles,
			   mfields_c_t *res)
{
  struct psc_moments_ops *ops = psc_moments_ops(moments);
  assert(ops && ops->calc_densities);
  ops->calc_densities(moments, flds, particles, res);
}

void
psc_moments_calc_v(struct psc_moments *moments,
		   mfields_base_t *flds, mparticles_base_t *particles,
		   mfields_c_t *res)
{
  struct psc_moments_ops *ops = psc_moments_ops(moments);
  assert(ops && ops->calc_v);
  ops->calc_v(moments, flds, particles, res);
}

void
psc_moments_calc_vv(struct psc_moments *moments,
		    mfields_base_t *flds, mparticles_base_t *particles,
		    mfields_c_t *res)
{
  struct psc_moments_ops *ops = psc_moments_ops(moments);
  assert(ops && ops->calc_vv);
  ops->calc_vv(moments, flds, particles, res);
}

void
psc_moments_calc_photon_n(struct psc_moments *moments,
			  mphotons_t *photons, mfields_c_t *res)
{
  struct psc_moments_ops *ops = psc_moments_ops(moments);
  assert(ops && ops->calc_photon_n);
  ops->calc_photon_n(moments, photons, res);
}

// ======================================================================
// psc_moments_create

static void
_psc_moments_create(struct psc_moments *moments)
{
  moments->bnd = psc_bnd_create(psc_moments_comm(moments));
  psc_bnd_set_name(moments->bnd, "psc_moments_bnd");
  psc_bnd_set_type(moments->bnd, "c");
  psc_bnd_set_psc(moments->bnd, ppsc);
  psc_moments_add_child(moments, (struct mrc_obj *) moments->bnd);
}

// ======================================================================
// psc_moments_init

static void
psc_moments_init()
{
  mrc_class_register_subclass(&mrc_class_psc_moments, &psc_moments_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_moments, &psc_moments_1st_ops);
  mrc_class_register_subclass(&mrc_class_psc_moments, &psc_moments_1st_cc_ops);
  mrc_class_register_subclass(&mrc_class_psc_moments, &psc_moments_fortran_ops);
}

// ======================================================================
// psc_moments class

struct mrc_class_psc_moments mrc_class_psc_moments = {
  .name             = "psc_moments",
  .size             = sizeof(struct psc_moments),
  .init             = psc_moments_init,
  .create           = _psc_moments_create,
};

