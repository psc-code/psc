
#include "psc.h"
#include "psc_particles_single.h"
#include "psc_particles_c.h"
#include "psc_particles_double.h"

#include <mrc_io.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// ======================================================================
// psc_particles "single"

static void
psc_particles_single_setup(struct psc_particles *prts)
{
  struct psc_particles_single *sub = psc_particles_single(prts);

  sub->n_alloced = psc_particles_size(prts) * 1.2;
  sub->particles = calloc(sub->n_alloced, sizeof(*sub->particles));
  sub->particles_alt = calloc(sub->n_alloced, sizeof(*sub->particles_alt));
  sub->b_idx = calloc(sub->n_alloced, sizeof(*sub->b_idx));
  sub->b_ids = calloc(sub->n_alloced, sizeof(*sub->b_ids));

  for (int d = 0; d < 3; d++) {
    sub->b_mx[d] = ppsc->patch[prts->p].ldims[d];
    sub->b_dxi[d] = 1.f / ppsc->patch[prts->p].dx[d];
  }
  sub->nr_blocks = sub->b_mx[0] * sub->b_mx[1] * sub->b_mx[2];
  sub->b_cnt = calloc(sub->nr_blocks + 1, sizeof(*sub->b_cnt));
}

static void
psc_particles_single_destroy(struct psc_particles *prts)
{
  struct psc_particles_single *sub = psc_particles_single(prts);

  free(sub->particles);
  free(sub->particles_alt);
  free(sub->b_idx);
  free(sub->b_ids);
  free(sub->b_cnt);
}

static void
psc_particles_single_reorder(struct psc_particles *prts)
{
  struct psc_particles_single *sub = psc_particles_single(prts);

  if (!sub->need_reorder) {
    return;
  }

  for (int n = 0; n < psc_particles_size(prts); n++) {
    sub->particles_alt[n] = sub->particles[sub->b_ids[n]];
  }
  
  // swap in alt array
  particle_single_t *tmp = sub->particles;
  sub->particles = sub->particles_alt;
  sub->particles_alt = tmp;
  sub->need_reorder = false;
}

void
particles_single_realloc(struct psc_particles *prts, int new_n_part)
{
  struct psc_particles_single *sub = psc_particles_single(prts);

  if (new_n_part <= sub->n_alloced)
    return;

  sub->n_alloced = new_n_part * 1.2;
  sub->particles = realloc(sub->particles, sub->n_alloced * sizeof(*sub->particles));
  sub->b_idx = realloc(sub->b_idx, sub->n_alloced * sizeof(*sub->b_idx));
  sub->b_ids = realloc(sub->b_ids, sub->n_alloced * sizeof(*sub->b_ids));
  free(sub->particles_alt);
  sub->particles_alt = malloc(sub->n_alloced * sizeof(*sub->particles_alt));
}

// ======================================================================

#ifdef HAVE_LIBHDF5_HL

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------
// psc_particles_single_write

static void
psc_particles_single_write(struct psc_particles *prts, struct mrc_io *io)
{
  int ierr;
  assert(sizeof(particle_single_t) / sizeof(particle_single_real_t) == 8);
  assert(sizeof(particle_single_real_t) == sizeof(float));

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, prts), H5P_DEFAULT); H5_CHK(group);
  // save/restore n_alloced, too?
  ierr = H5LTset_attribute_int(group, ".", "p", &prts->p, 1); CE;
  int n_prts = psc_particles_size(prts);
  ierr = H5LTset_attribute_int(group, ".", "n_part", &n_prts, 1); CE;
  ierr = H5LTset_attribute_uint(group, ".", "flags", &prts->flags, 1); CE;
  if (n_prts > 0) {
    // in a rather ugly way, we write the int "kind" member as a float
    hsize_t hdims[2] = { n_prts, 8 };
    ierr = H5LTmake_dataset_float(group, "particles_single", 2, hdims,
				  (float *) particles_single_get_one(prts, 0)); CE;
  }
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_particles_single_read

static void
psc_particles_single_read(struct psc_particles *prts, struct mrc_io *io)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, prts), H5P_DEFAULT); H5_CHK(group);
  ierr = H5LTget_attribute_int(group, ".", "p", &prts->p); CE;
  int n_prts;
  ierr = H5LTget_attribute_int(group, ".", "n_part", &n_prts); CE;
  psc_particles_resize(prts, n_prts);
  ierr = H5LTget_attribute_uint(group, ".", "flags", &prts->flags); CE;
  psc_particles_setup(prts);
  if (n_prts > 0) {
    ierr = H5LTread_dataset_float(group, "particles_single",
				  (float *) particles_single_get_one(prts, 0)); CE;
  }
  ierr = H5Gclose(group); CE;
}

#endif

// ======================================================================

static inline void
calc_vxi(particle_single_real_t vxi[3], particle_single_t *part)
{
  particle_single_real_t root =
    1.f / sqrtf(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

static void
psc_particles_single_copy_to_c(struct psc_particles *prts_base,
			       struct psc_particles *prts_c, unsigned int flags)
{
  particle_single_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }

  int n_prts = psc_particles_size(prts_base);
  psc_particles_resize(prts_c, n_prts);
  assert(n_prts <= psc_particles_c(prts_c)->n_alloced);
  for (int n = 0; n < n_prts; n++) {
    particle_single_t *part_base = particles_single_get_one(prts_base, n);
    particle_c_t *part = particles_c_get_one(prts_c, n);
    
    particle_c_real_t qni = ppsc->kinds[part_base->kind].q;
    particle_c_real_t mni = ppsc->kinds[part_base->kind].m;
    particle_c_real_t wni = part_base->qni_wni / qni;
    
    particle_single_real_t vxi[3];
    calc_vxi(vxi, part_base);
    part->xi  = part_base->xi - dth[0] * vxi[0];
    part->yi  = part_base->yi - dth[1] * vxi[1];
    part->zi  = part_base->zi - dth[2] * vxi[2];
    part->pxi = part_base->pxi;
    part->pyi = part_base->pyi;
    part->pzi = part_base->pzi;
    part->qni = qni;
    part->mni = mni;
    part->wni = wni;
    part->kind = part_base->kind;
  }
}

static void
psc_particles_single_copy_from_c(struct psc_particles *prts_base,
				 struct psc_particles *prts_c, unsigned int flags)
{
  particle_single_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }

  struct psc_particles_single *sub = psc_particles_single(prts_base);
  int n_prts = psc_particles_size(prts_c);
  psc_particles_resize(prts_base, n_prts);
  assert(n_prts <= sub->n_alloced);
  for (int n = 0; n < n_prts; n++) {
    particle_single_t *part_base = particles_single_get_one(prts_base, n);
    particle_c_t *part = particles_c_get_one(prts_c, n);
    
    particle_single_real_t qni_wni;
    if (part->qni != 0.) {
      qni_wni = part->qni * part->wni;
    } else {
      qni_wni = part->wni;
    }
    
    part_base->xi          = part->xi;
    part_base->yi          = part->yi;
    part_base->zi          = part->zi;
    part_base->pxi         = part->pxi;
    part_base->pyi         = part->pyi;
    part_base->pzi         = part->pzi;
    part_base->qni_wni     = qni_wni;
    part_base->kind        = part->kind;

    particle_single_real_t vxi[3];
    calc_vxi(vxi, part_base);
    part_base->xi += dth[0] * vxi[0];
    part_base->yi += dth[1] * vxi[1];
    part_base->zi += dth[2] * vxi[2];
  }
}

// ----------------------------------------------------------------------
// psc_particles_single_copy_to_double

static void
psc_particles_single_copy_to_double(struct psc_particles *prts_base,
				    struct psc_particles *prts, unsigned int flags)
{
  int n_prts = psc_particles_size(prts_base);
  psc_particles_resize(prts, n_prts);
  assert(n_prts <= psc_particles_double(prts)->n_alloced);
  for (int n = 0; n < n_prts; n++) {
    particle_single_t *part_base = particles_single_get_one(prts_base, n);
    particle_double_t *part = particles_double_get_one(prts, n);
    
    part->xi      = part_base->xi;
    part->yi      = part_base->yi;
    part->zi      = part_base->zi;
    part->pxi     = part_base->pxi;
    part->pyi     = part_base->pyi;
    part->pzi     = part_base->pzi;
    part->qni_wni = part_base->qni_wni;
    part->kind    = part_base->kind;
  }
}

// ----------------------------------------------------------------------
// psc_particles_single_copy_from_double

static void
psc_particles_single_copy_from_double(struct psc_particles *prts_base,
				      struct psc_particles *prts, unsigned int flags)
{
  struct psc_particles_single *sub = psc_particles_single(prts_base);
  int n_prts = psc_particles_size(prts);
  psc_particles_resize(prts_base, n_prts);
  assert(n_prts <= sub->n_alloced);
  for (int n = 0; n < n_prts; n++) {
    particle_single_t *part_base = particles_single_get_one(prts_base, n);
    particle_double_t *part = particles_double_get_one(prts, n);

    part_base->xi      = part->xi;
    part_base->yi      = part->yi;
    part_base->zi      = part->zi;
    part_base->pxi     = part->pxi;
    part_base->pyi     = part->pyi;
    part_base->pzi     = part->pzi;
    part_base->qni_wni = part->qni_wni;
    part_base->kind    = part->kind;
  }
}

// ======================================================================
// psc_particles: subclass "single"

static struct mrc_obj_method psc_particles_single_methods[] = {
  MRC_OBJ_METHOD("copy_to_c"       , psc_particles_single_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c"     , psc_particles_single_copy_from_c),
  MRC_OBJ_METHOD("copy_to_double"  , psc_particles_single_copy_to_double),
  MRC_OBJ_METHOD("copy_from_double", psc_particles_single_copy_from_double),
  {}
};

struct psc_particles_ops psc_particles_single_ops = {
  .name                    = "single",
  .size                    = sizeof(struct psc_particles_single),
  .setup                   = psc_particles_single_setup,
  .destroy                 = psc_particles_single_destroy,
#ifdef HAVE_LIBHDF5_HL
  .read                    = psc_particles_single_read,
  .write                   = psc_particles_single_write,
#endif
  .reorder                 = psc_particles_single_reorder,
};

// ======================================================================
// psc_mparticles: subclass "single"
  
struct psc_mparticles_ops psc_mparticles_single_ops = {
  .name                    = "single",
  .methods                 = psc_particles_single_methods,
};

