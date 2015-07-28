
#include "psc.h"
#include "psc_particles_private.h"

#include <mrc_profile.h>
#include <string.h>

// ======================================================================
// psc_particles_reorder

void
psc_particles_reorder(struct psc_particles *prts)
{
  if (psc_particles_ops(prts)->reorder) {
    psc_particles_ops(prts)->reorder(prts);
  }
}

// ======================================================================
// psc_particles_get_as

struct psc_particles *
psc_particles_get_as(struct psc_particles *prts_base, const char *type,
		     unsigned int flags)
{
  const char *type_base = psc_particles_type(prts_base);
  // If we're already the subtype, nothing to be done
  if (strcmp(type_base, type) == 0)
    return prts_base;

  static int pr;
  if (!pr) {
    pr = prof_register("particles_get_as", 1., 0, 0);
  }
  prof_start(pr);

  struct psc_particles *prts = psc_particles_create(psc_particles_comm(prts_base));
  psc_particles_set_type(prts, type);
  prts->n_part = prts_base->n_part;
  prts->flags = flags;
  prts->p = prts_base->p;
  psc_particles_setup(prts);

  if (!(flags & MP_DONT_COPY)) {
    char s[strlen(type) + 12]; sprintf(s, "copy_to_%s", type);
    psc_particles_copy_to_func_t copy_to = (psc_particles_copy_to_func_t)
      psc_particles_get_method(prts_base, s);
    if (copy_to) {
      copy_to(prts_base, prts, flags);
    } else {
      sprintf(s, "copy_from_%s", type_base);
      psc_particles_copy_to_func_t copy_from = (psc_particles_copy_from_func_t)
	psc_particles_get_method(prts, s);
      if (copy_from) {
	copy_from(prts, prts_base, flags);
      } else {
	fprintf(stderr, "ERROR: no 'copy_to_%s' in psc_particles '%s' and "
		"no 'copy_from_%s' in '%s'!\n",
		type, psc_particles_type(prts_base), type_base, psc_particles_type(prts));
	assert(0);
      }
    }
  }

  prof_stop(pr);
  return prts;
}

void
psc_particles_put_as(struct psc_particles *prts, struct psc_particles *prts_base,
		     unsigned int flags)
{
  // If we're already the subtype, nothing to be done
  const char *type = psc_particles_type(prts);
  const char *type_base = psc_particles_type(prts_base);
  if (strcmp(type_base, type) == 0)
    return;

  static int pr;
  if (!pr) {
    pr = prof_register("particles_put_as", 1., 0, 0);
  }
  prof_start(pr);

  if (!(flags & MP_DONT_COPY)) {
    char s[strlen(type) + 12]; sprintf(s, "copy_from_%s", type);
    psc_particles_copy_from_func_t copy_from = (psc_particles_copy_from_func_t)
      psc_particles_get_method(prts_base, s);
    if (copy_from) {
      copy_from(prts_base, prts, MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS);
    } else {
      sprintf(s, "copy_to_%s", type_base);
      psc_particles_copy_from_func_t copy_to = (psc_particles_copy_from_func_t)
	psc_particles_get_method(prts, s);
      if (copy_to) {
	copy_to(prts, prts_base, MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS);
      } else {
	fprintf(stderr, "ERROR: no 'copy_from_%s' in psc_particles '%s' and "
		"no 'copy_to_%s' in '%s'!\n",
		type, psc_particles_type(prts_base), type_base, psc_particles_type(prts));
	assert(0);
      }
    }
  }
  psc_particles_destroy(prts);

  prof_stop(pr);
}

// ======================================================================
// psc_particles_init

static void
psc_particles_init()
{
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_fortran_ops);
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_single_by_block_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_cuda_ops);
#endif
#ifdef USE_CUDA2
  mrc_class_register_subclass(&mrc_class_psc_particles, &psc_particles_cuda2_ops);
#endif
}

// ======================================================================
// psc_particles class

struct mrc_class_psc_particles mrc_class_psc_particles = {
  .name             = "psc_particles",
  .size             = sizeof(struct psc_particles),
  .init             = psc_particles_init,
};

