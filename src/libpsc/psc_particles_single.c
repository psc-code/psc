
#include "psc.h"
#include "psc_particles_as_single.h"
#include "psc_particles_inc.h"
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

  prts->n_alloced = psc_particles_size(prts) * 1.2;
  sub->particles = calloc(prts->n_alloced, sizeof(*sub->particles));
  sub->particles_alt = calloc(prts->n_alloced, sizeof(*sub->particles_alt));
  sub->b_idx = calloc(prts->n_alloced, sizeof(*sub->b_idx));
  sub->b_ids = calloc(prts->n_alloced, sizeof(*sub->b_ids));

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

  if (new_n_part <= prts->n_alloced)
    return;

  prts->n_alloced = new_n_part * 1.2;
  sub->particles = realloc(sub->particles, prts->n_alloced * sizeof(*sub->particles));
  sub->b_idx = realloc(sub->b_idx, prts->n_alloced * sizeof(*sub->b_idx));
  sub->b_ids = realloc(sub->b_ids, prts->n_alloced * sizeof(*sub->b_ids));
  free(sub->particles_alt);
  sub->particles_alt = malloc(prts->n_alloced * sizeof(*sub->particles_alt));
}

// ======================================================================
// conversion to/from "c"

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
get_particle_c(particle_single_t *prt, int n, struct psc_particles *prts_c)
{
  particle_single_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  
  particle_c_t *prt_c = particles_c_get_one(prts_c, n);

  prt->xi      = prt_c->xi;
  prt->yi      = prt_c->yi;
  prt->zi      = prt_c->zi;
  prt->pxi     = prt_c->pxi;
  prt->pyi     = prt_c->pyi;
  prt->pzi     = prt_c->pzi;
  prt->kind    = prt_c->kind;
  prt->qni_wni = prt_c->qni * prt_c->wni;

  particle_single_real_t vxi[3];
  calc_vxi(vxi, prt);
  prt->xi += dth[0] * vxi[0];
  prt->yi += dth[1] * vxi[1];
  prt->zi += dth[2] * vxi[2];
}

static void
put_particle_c(particle_single_t *prt, int n, struct psc_particles *prts_c)
{
  particle_single_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  
  particle_single_real_t vxi[3];
  calc_vxi(vxi, prt);

  particle_c_t *prt_c = particles_c_get_one(prts_c, n);

  particle_c_real_t qni = ppsc->kinds[prt->kind].q;
  particle_c_real_t mni = ppsc->kinds[prt->kind].m;
  particle_c_real_t wni = prt->qni_wni / qni;

  prt_c->xi      = prt->xi - dth[0] * vxi[0];
  prt_c->yi      = prt->yi - dth[1] * vxi[1];
  prt_c->zi      = prt->zi - dth[2] * vxi[2];
  prt_c->pxi     = prt->pxi;
  prt_c->pyi     = prt->pyi;
  prt_c->pzi     = prt->pzi;
  prt_c->kind    = prt->kind;
  prt_c->qni     = qni;
  prt_c->wni     = wni;
  prt_c->mni     = mni;
}

static void
psc_mparticles_single_copy_to_c(int p, struct psc_mparticles *mprts,
				struct psc_mparticles *mprts_c, unsigned int flags)
{
  psc_mparticles_copy_to(p, mprts, mprts_c, flags, put_particle_c);
}

static void
psc_mparticles_single_copy_from_c(int p, struct psc_mparticles *mprts,
				  struct psc_mparticles *mprts_c, unsigned int flags)
{
  psc_mparticles_copy_from(p, mprts, mprts_c, flags, get_particle_c);
}

// ======================================================================
// conversion to/from "double"

static void
put_particle_double(particle_single_t *prt, int n, struct psc_particles *prts_dbl)
{
  particle_double_t *prt_dbl = particles_double_get_one(prts_dbl, n);
  
  prt_dbl->xi      = prt->xi;
  prt_dbl->yi      = prt->yi;
  prt_dbl->zi      = prt->zi;
  prt_dbl->pxi     = prt->pxi;
  prt_dbl->pyi     = prt->pyi;
  prt_dbl->pzi     = prt->pzi;
  prt_dbl->qni_wni = prt->qni_wni;
  prt_dbl->kind    = prt->kind;
}

static void
get_particle_double(particle_single_t *prt, int n, struct psc_particles *prts_dbl)
{
  particle_double_t *prt_dbl = particles_double_get_one(prts_dbl, n);

  prt->xi      = prt_dbl->xi;
  prt->yi      = prt_dbl->yi;
  prt->zi      = prt_dbl->zi;
  prt->pxi     = prt_dbl->pxi;
  prt->pyi     = prt_dbl->pyi;
  prt->pzi     = prt_dbl->pzi;
  prt->qni_wni = prt_dbl->qni_wni;
  prt->kind    = prt_dbl->kind;
}

static void
psc_mparticles_single_copy_to_double(int p, struct psc_mparticles *mprts,
				    struct psc_mparticles *mprts_dbl, unsigned int flags)
{
  psc_mparticles_copy_to(p, mprts, mprts_dbl, flags, put_particle_double);
}

static void
psc_mparticles_single_copy_from_double(int p, struct psc_mparticles *mprts,
				       struct psc_mparticles *mprts_dbl, unsigned int flags)
{
  psc_mparticles_copy_from(p, mprts, mprts_dbl, flags, get_particle_double);
}

// ======================================================================
// psc_particles: subclass "single"

struct psc_particles_ops psc_particles_single_ops = {
  .name                    = "single",
  .size                    = sizeof(struct psc_particles_single),
  .setup                   = psc_particles_single_setup,
  .destroy                 = psc_particles_single_destroy,
  .reorder                 = psc_particles_single_reorder,
};

// ======================================================================
// psc_mparticles: subclass "single"

#ifdef HAVE_LIBHDF5_HL

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------
// psc_mparticles_single_write

static void
psc_mparticles_single_write(struct psc_mparticles *mprts, struct mrc_io *io)
{
  int ierr;
  assert(sizeof(particle_t) / sizeof(particle_real_t) == 8);
  assert(sizeof(particle_real_t) == sizeof(float));

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *_prts = psc_mparticles_get_patch(mprts, p);
    particle_range_t prts = particle_range_mprts(mprts, p);
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gcreate(group, pname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(pgroup);
    // save/restore n_alloced, too?
    ierr = H5LTset_attribute_uint(pgroup, ".", "flags", &_prts->flags, 1); CE;
    int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);
    ierr = H5LTset_attribute_int(pgroup, ".", "n_prts", &n_prts, 1); CE;
    if (n_prts > 0) {
      // in a rather ugly way, we write the int "kind" member as a float
      hsize_t hdims[2] = { n_prts, 8 };
      ierr = H5LTmake_dataset_float(pgroup, "data", 2, hdims,
				    (float *) particle_iter_deref(prts.begin)); CE;
    }
    ierr = H5Gclose(pgroup); CE;
  }
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_mparticles_single_read

static void
psc_mparticles_single_read(struct psc_mparticles *mprts, struct mrc_io *io)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);
  // FIXME those should be superclass bits
  mprts->domain = mrc_io_read_ref(io, mprts, "domain", mrc_domain);
  mrc_domain_get_patches(mprts->domain, &mprts->nr_patches);
  mrc_io_read_int(io, mprts, "flags", (int *) &mprts->flags);

  mprts->prts = calloc(mprts->nr_patches, sizeof(*mprts->prts));
  /* mprts->nr_particles_by_patch = */
  /*   calloc(mprts->nr_patches, sizeof(*mprts->nr_particles_by_patch)); */

  for (int p = 0; p < mprts->nr_patches; p++) {
    char name[20]; sprintf(name, "prts%d", p);
    struct psc_particles *_prts = psc_particles_create(MPI_COMM_NULL);
    psc_particles_set_type(_prts, "single");
    psc_particles_set_name(_prts, name);
    _prts->p = p;
    mprts->prts[p] = _prts;

    particle_range_t prts = particle_range_mprts(mprts, p);
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    ierr = H5LTget_attribute_uint(pgroup, ".", "flags", &_prts->flags); CE;
    int n_prts;
    ierr = H5LTget_attribute_int(pgroup, ".", "n_prts", &n_prts); CE;
    psc_particles_set_n_prts(_prts, n_prts);
    psc_particles_setup(_prts);
    
    if (n_prts > 0) {
      ierr = H5LTread_dataset_float(pgroup, "data",
				    (float *) particle_iter_deref(prts.begin)); CE;
    }
    ierr = H5Gclose(pgroup); CE;
  }
  ierr = H5Gclose(group); CE;
  psc_mparticles_view(mprts);
}

#endif

// ----------------------------------------------------------------------
// psc_mparticles subclass "single"
  
static struct mrc_obj_method psc_mparticles_single_methods[] = {
  MRC_OBJ_METHOD("copy_to_c"       , psc_mparticles_single_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c"     , psc_mparticles_single_copy_from_c),
  MRC_OBJ_METHOD("copy_to_double"  , psc_mparticles_single_copy_to_double),
  MRC_OBJ_METHOD("copy_from_double", psc_mparticles_single_copy_from_double),
  {}
};

struct psc_mparticles_ops psc_mparticles_single_ops = {
  .name                    = "single",
  .methods                 = psc_mparticles_single_methods,
#ifdef HAVE_LIBHDF5_HL
  .write                   = psc_mparticles_single_write,
  .read                    = psc_mparticles_single_read,
#endif
};

