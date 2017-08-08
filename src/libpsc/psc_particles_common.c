
#include <mrc_io.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// ======================================================================
// psc_mparticles "single" / "double" / "c"

// ----------------------------------------------------------------------
// psc_mparticles_sub_setup_patch

static void
MPFX(setup_patch)(struct psc_mparticles *mprts, int p, int n_prts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);
  struct MPFX(patch) *patch = &sub->patch[p];

  patch->n_prts = n_prts;
  int n_alloced = psc_mparticles_n_prts_by_patch(mprts, p) * 1.2;
  patch->n_alloced = n_alloced;
  patch->prt_array = calloc(n_alloced, sizeof(*patch->prt_array));

#if PSC_PARTICLES_AS_SINGLE
  patch->prt_array_alt = calloc(n_alloced, sizeof(*patch->prt_array_alt));
  patch->b_idx = calloc(n_alloced, sizeof(*patch->b_idx));
  patch->b_ids = calloc(n_alloced, sizeof(*patch->b_ids));

  int *b_mx = patch->b_mx;
  for (int d = 0; d < 3; d++) {
    b_mx[d] = ppsc->patch[p].ldims[d];
    patch->b_dxi[d] = 1.f / ppsc->patch[p].dx[d];
  }
  patch->nr_blocks = b_mx[0] * b_mx[1] * b_mx[2];
  patch->b_cnt = calloc(patch->nr_blocks + 1, sizeof(*patch->b_cnt));
#endif

#if PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  patch->prt_array_alt = calloc(n_alloced, sizeof(*patch->prt_array_alt));
  patch->b_idx = calloc(n_alloced, sizeof(*patch->b_idx));
  patch->b_ids = calloc(n_alloced, sizeof(*patch->b_ids));

  for (int d = 0; d < 3; d++) {
    patch->b_mx[d] = ppsc->patch[p].ldims[d];
    patch->b_dxi[d] = 1.f / ppsc->patch[p].dx[d];
  }
  patch->nr_blocks = patch->b_mx[0] * patch->b_mx[1] * patch->b_mx[2];
  patch->b_cnt = calloc(patch->nr_blocks + 1, sizeof(*patch->b_cnt));
  patch->b_off = calloc(patch->nr_blocks + 2, sizeof(*patch->b_off));
#endif
}

// ----------------------------------------------------------------------
// psc_mparticles_sub_destroy_patch

static void
MPFX(destroy_patch)(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);
  struct MPFX(patch) *patch = &sub->patch[p];

  free(patch->prt_array);

#if PSC_PARTICLES_AS_SINGLE
  free(patch->prt_array_alt);
  free(patch->b_idx);
  free(patch->b_ids);
  free(patch->b_cnt);
#endif

#if PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  free(patch->prt_array_alt);
  free(patch->b_idx);
  free(patch->b_ids);
  free(patch->b_cnt);
  free(patch->b_off);
#endif
}

// ----------------------------------------------------------------------
// psc_mparticles_sub_realloc_patch

static void
MPFX(realloc_patch)(struct psc_mparticles *mprts, int p, int new_n_prts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);
  struct MPFX(patch) *patch = &sub->patch[p];

  if (new_n_prts <= patch->n_alloced)
    return;

  int n_alloced = new_n_prts * 1.2;
  patch->n_alloced = n_alloced;

  patch->prt_array = realloc(patch->prt_array, n_alloced * sizeof(*patch->prt_array));

#if PSC_PARTICLES_AS_SINGLE
  free(patch->prt_array_alt);
  patch->prt_array_alt = malloc(n_alloced * sizeof(*patch->prt_array_alt));
  patch->b_idx = realloc(patch->b_idx, n_alloced * sizeof(*patch->b_idx));
  patch->b_ids = realloc(patch->b_ids, n_alloced * sizeof(*patch->b_ids));
#endif

#if PSC_PARTICLES_AS_SINGLE_BY_BLOCK
  free(patch->prt_array_alt);
  patch->prt_array_alt = malloc(n_alloced * sizeof(*patch->prt_array_alt));
  patch->b_idx = realloc(patch->b_idx, n_alloced * sizeof(*patch->b_idx));
  patch->b_ids = realloc(patch->b_ids, n_alloced * sizeof(*patch->b_ids));
#endif
}

// ----------------------------------------------------------------------
// psc_mparticles_sub_setup
//
// FIXME does deliberately not call old-style setup_super(), rather duplicates
// that code

static void
MPFX(setup)(struct psc_mparticles *mprts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);
  
  assert(mprts->nr_particles_by_patch);

  mprts->mpatch = calloc(mprts->nr_patches, sizeof(*mprts->mpatch));
  sub->patch = calloc(mprts->nr_patches, sizeof(*sub->patch));
  for (int p = 0; p < mprts->nr_patches; p++) {
    MPFX(setup_patch)(mprts, p, mprts->nr_particles_by_patch[p]);
  }

  free(mprts->nr_particles_by_patch);
  mprts->nr_particles_by_patch = NULL;
}

// ----------------------------------------------------------------------
// psc_mparticls_sub_write/read
  
#if (PSC_PARTICLES_AS_DOUBLE || PSC_PARTICLES_AS_SINGLE) && HAVE_LIBHDF5_HL

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------
// psc_mparticles_sub_write

static void
MPFX(write)(struct psc_mparticles *mprts, struct mrc_io *io)
{
  int ierr;
  assert(sizeof(particle_t) / sizeof(particle_real_t) == 8);

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);
  for (int p = 0; p < mprts->nr_patches; p++) {
    particle_range_t prts = particle_range_mprts(mprts, p);
    char pname[10]; sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gcreate(group, pname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);
    ierr = H5LTset_attribute_int(pgroup, ".", "n_prts", &n_prts, 1); CE;
    if (n_prts > 0) {
      // in a rather ugly way, we write the int "kind" member as a float / double
      hsize_t hdims[2] = { n_prts, 8 };
#if PSC_PARTICLES_AS_DOUBLE
      ierr = H5LTmake_dataset_double(pgroup, "data", 2, hdims,
				    (double *) particle_iter_deref(prts.begin)); CE;
#elif PSC_PARTICLES_AS_SINGLE
      ierr = H5LTmake_dataset_float(pgroup, "data", 2, hdims,
				    (float *) particle_iter_deref(prts.begin)); CE;
#else
      assert(0);
#endif
    }
    ierr = H5Gclose(pgroup); CE;
  }
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_mparticles_sub_read

static void
MPFX(read)(struct psc_mparticles *mprts, struct mrc_io *io)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);
  // FIXME those should be superclass bits
  mprts->domain = mrc_io_read_ref(io, mprts, "domain", mrc_domain);
  mrc_domain_get_patches(mprts->domain, &mprts->nr_patches);
  mrc_io_read_int(io, mprts, "flags", (int *) &mprts->flags);

  mprts->mpatch = calloc(mprts->nr_patches, sizeof(*mprts->mpatch));

  for (int p = 0; p < mprts->nr_patches; p++) {
    particle_range_t prts = particle_range_mprts(mprts, p);
    char pname[10]; sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts;
    ierr = H5LTget_attribute_int(pgroup, ".", "n_prts", &n_prts); CE;
    MPFX(setup_patch)(mprts, p, n_prts);
    
    if (n_prts > 0) {
#if PSC_PARTICLES_AS_SINGLE
      ierr = H5LTread_dataset_float(pgroup, "data",
				    (float *) particle_iter_deref(prts.begin)); CE;
#elif PSC_PARTICLES_AS_DOUBLE
      ierr = H5LTread_dataset_double(pgroup, "data",
				    (double *) particle_iter_deref(prts.begin)); CE;
#else
      assert(0);
#endif
    }
    ierr = H5Gclose(pgroup); CE;
  }
  ierr = H5Gclose(group); CE;
}

#else

static void
MPFX(write)(struct psc_mparticles *mprts, struct mrc_io *io)
{
  assert(0);
}

static void
MPFX(read)(struct psc_mparticles *mprts, struct mrc_io *io)
{
  assert(0);
}

#endif

static void
MPFX(destroy)(struct psc_mparticles *mprts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  for (int p = 0; p < mprts->nr_patches; p++) {
    MPFX(destroy_patch)(mprts, p);
  }
  free(sub->patch);
}

static int
MPFX(get_n_alloced)(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  return sub->patch[p].n_alloced;
}

static void
MPFX(resize_patch)(struct psc_mparticles *mprts, int p, int n_prts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  assert(n_prts <= sub->patch[p].n_alloced);
  sub->patch[p].n_prts = n_prts;
}

static void
MPFX(set_n_prts)(struct psc_mparticles *mprts, int p, int n_prts)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  sub->patch[p].n_prts = n_prts;
}

static int
MPFX(get_n_prts)(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  return sub->patch[p].n_prts;
}

static void
MPFX(get_n_prts_all)(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct psc_mparticles_sub *sub = psc_mparticles_sub(mprts);

  for (int p = 0; p < mprts->nr_patches; p++) {
    n_prts_by_patch[p] = sub->patch[p].n_prts;
  }
}

// ----------------------------------------------------------------------
// psc_mparticles_ops

struct psc_mparticles_ops MPFX(ops) = {
  .name                    = PARTICLE_TYPE,
  .size                    = sizeof(struct psc_mparticles_sub),
  .methods                 = MPFX(methods),
  .setup                   = MPFX(setup),
  .destroy                 = MPFX(destroy),
  .write                   = MPFX(write),
  .read                    = MPFX(read),
  .realloc                 = MPFX(realloc_patch),
  .get_n_alloced           = MPFX(get_n_alloced),
  .resize_patch            = MPFX(resize_patch),
  .set_n_prts              = MPFX(set_n_prts),
  .get_n_prts              = MPFX(get_n_prts),
  .get_n_prts_all          = MPFX(get_n_prts_all),
};

