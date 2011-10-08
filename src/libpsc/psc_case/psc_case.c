
#include "psc_case_private.h"

#include <psc_push_fields.h>
#include <psc_bnd_fields.h>
#include <psc_balance.h>
#include <mrc_params.h>
#include <stdlib.h>

// ----------------------------------------------------------------------
// psc_init_field_pml helper

static void
psc_init_field_pml(struct psc_case *_case, mfields_base_t *flds)
{
  psc_mfields_copy_comp(flds, DX, flds, EX);
  psc_mfields_copy_comp(flds, DY, flds, EY);
  psc_mfields_copy_comp(flds, DZ, flds, EZ);
  psc_mfields_copy_comp(flds, BX, flds, HX);
  psc_mfields_copy_comp(flds, BY, flds, HY);
  psc_mfields_copy_comp(flds, BZ, flds, HZ);
  psc_mfields_set_comp(flds, EPS, 1.);
  psc_mfields_set_comp(flds, MU, 1.);
}

// ----------------------------------------------------------------------
// _psc_case_create

static void
_psc_case_create(struct psc_case *_case)
{
  _case->psc = psc_create(psc_case_comm(_case));
}

// ----------------------------------------------------------------------
// _psc_case_destroy

static void
_psc_case_destroy(struct psc_case *_case)
{
  psc_destroy(_case->psc);
  ppsc = NULL;
}

// ----------------------------------------------------------------------
// _psc_case_set_from_options

static void
_psc_case_set_from_options(struct psc_case *_case)
{
  // subclass set_from_options gets called first automatically 
  psc_set_from_options(_case->psc);
}

// ----------------------------------------------------------------------
// _psc_case_setup

static void
_psc_case_setup(struct psc_case *_case)
{
  struct psc *psc = _case->psc;

  psc_setup_coeff(psc);
  psc_setup_domain(psc); // needs to be done before setting up psc_bnd

  //TODO: Set this somewhere less hackish
  psc->patchmanager.currentcase = _case;
  if (psc_case_ops(_case)->setup != NULL) psc_case_ops(_case)->setup(_case);
  
  // alloc / initialize particles
  int particle_label_offset;
  int *nr_particles_by_patch = calloc(psc->nr_patches, sizeof(*nr_particles_by_patch));
  psc_case_init_partition(_case, nr_particles_by_patch, &particle_label_offset);
  psc_balance_initial(psc->balance, psc, &nr_particles_by_patch);

  psc->particles = 
    psc_mparticles_base_create(mrc_domain_comm(psc->mrc_domain));
  psc_mparticles_base_set_type(psc->particles, s_particles_base);
  psc_mparticles_base_set_name(psc->particles, "mparticles");
  psc_mparticles_base_set_domain_nr_particles(psc->particles, psc->mrc_domain,
					  nr_particles_by_patch);
  psc_mparticles_base_setup(psc->particles);

  psc_case_init_particles(_case, nr_particles_by_patch, particle_label_offset);
  free(nr_particles_by_patch);

  // alloc / initialize fields
  psc->flds = psc_mfields_create(mrc_domain_comm(psc->mrc_domain));
  psc_mfields_list_add(&psc_mfields_base_list, &psc->flds);
  psc_mfields_set_type(psc->flds, psc->prm.fields_base);
  psc_mfields_set_name(psc->flds, "mfields");
  psc_mfields_set_domain(psc->flds, psc->mrc_domain);
  psc_mfields_set_param_int(psc->flds, "nr_fields", NR_FIELDS);
  psc_mfields_set_param_int3(psc->flds, "ibn", psc->ibn);
  psc_mfields_setup(psc->flds);
  psc_case_init_field(_case, psc->flds);
  if (psc->domain.use_pml) {
    psc_init_field_pml(_case, psc->flds);
  }

  // alloc / initialize photons
  psc_case_init_photons(_case);

  // this sets up everything except allocating fields and particles,
  // and intializing them
  psc_setup(psc);

  psc_setup_fortran(psc);
}

// ----------------------------------------------------------------------
// _psc_case_view

static void
_psc_case_view(struct psc_case *_case)
{
  psc_view(_case->psc);
}

// ----------------------------------------------------------------------
// psc_case_init_npt

void
psc_case_init_npt(struct psc_case *_case, int kind, double x[3],
		  struct psc_particle_npt *npt)
{
  if (!psc_case_ops(_case)->init_npt)
    return;

  psc_case_ops(_case)->init_npt(_case, kind, x, npt);
}

// ----------------------------------------------------------------------
// psc_case_init_field

void
psc_case_init_field(struct psc_case *_case, mfields_base_t *flds_base)
{
  mfields_c_t *flds = psc_mfields_get_c(flds_base, JXI, HX + 3);

  // case-specific other initial condition
  if (psc_case_ops(_case)->init_field) {
    psc_case_ops(_case)->init_field(_case, flds);
  }

  psc_mfields_put_c(flds, flds_base, JXI, HX + 3);
}

// ----------------------------------------------------------------------
// psc_case_integrate

void
psc_case_integrate(struct psc_case *_case)
{
  if (psc_case_ops(_case)->integrate) {
    psc_case_ops(_case)->integrate(_case);
  } else {
    psc_integrate(_case->psc);
  }
}

// ----------------------------------------------------------------------
// psc_case_get_psc

struct psc *
psc_case_get_psc(struct psc_case *_case)
{
  return _case->psc;
}

// ======================================================================
// psc_case_init

static struct psc_case_ops psc_case_default_ops = {
  .name                  = "default",
};

static void
psc_case_init()
{
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_default_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_test_xy_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_test_yz_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_test_xz_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_test_z_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_harris_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_harris_xy_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_langmuir_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_wakefield_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_thinfoil_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_foils_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_curvedfoil_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_singlepart_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_collisions_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_cone_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_microsphere_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_photon_test_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_bubble_ops);
  mrc_class_register_subclass(&mrc_class_psc_case, &psc_case_dynamic_ops);
}

// ======================================================================
// psc_case class

#define VAR(x) (void *)offsetof(struct psc_case, x)
static struct param psc_case_descr[] = {
  { "seed_by_time"          , VAR(seed_by_time)     , PARAM_BOOL(false)   },
  { "nr_kinds"              , VAR(nr_kinds)         , PARAM_INT(2)        },
  {},
};
#undef VAR

struct mrc_class_psc_case mrc_class_psc_case = {
  .name             = "psc_case",
  .size             = sizeof(struct psc_case),
  .init             = psc_case_init,
  .param_descr      = psc_case_descr,
  .create           = _psc_case_create,
  .destroy          = _psc_case_destroy,
  .set_from_options = _psc_case_set_from_options,
  .setup            = _psc_case_setup,
  .view             = _psc_case_view,
};

