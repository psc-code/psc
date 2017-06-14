
#include "psc.h"
#include "psc_cuda.h"
#include "psc_bnd_cuda.h"
#include "psc_particles_cuda.h"
#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "psc_push_particles.h"
#include "particles_cuda.h"

#include "cuda_mparticles.h"

#include <mrc_io.h>

EXTERN_C void cuda_init(int rank);

// ======================================================================
// psc_particles "cuda"

static void
psc_particles_cuda_setup(struct psc_particles *prts)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

  struct psc_patch *patch = &ppsc->patch[prts->p];

  if (!prts->flags) {
    // FIXME, they get set too early, so auto-dispatch "1vb" doesn't work
    prts->flags = MP_NEED_BLOCK_OFFSETS | MP_BLOCKSIZE_4X4X4 | MP_NO_CHECKERBOARD;
  }

  int bs[3];
  for (int d = 0; d < 3; d++) {
    switch (prts->flags & MP_BLOCKSIZE_MASK) {
    case MP_BLOCKSIZE_1X1X1: bs[d] = 1; break;
    case MP_BLOCKSIZE_2X2X2: bs[d] = 2; break;
    case MP_BLOCKSIZE_4X4X4: bs[d] = 4; break;
    case MP_BLOCKSIZE_8X8X8: bs[d] = 8; break;
    default: assert(0);
    }
    if (ppsc->domain.gdims[d] == 1) {
      bs[d] = 1;
    }
    assert(patch->ldims[d] % bs[d] == 0); // not sure what breaks if not
    cuda->b_mx[d] = (patch->ldims[d] + bs[d] - 1) / bs[d];
    cuda->b_dxi[d] = 1.f / (bs[d] * ppsc->patch[prts->p].dx[d]);
  }
  cuda->nr_blocks = cuda->b_mx[0] * cuda->b_mx[1] * cuda->b_mx[2];

  for (int d = 0; d < 3; d++) {
    if (prts->flags & MP_NO_CHECKERBOARD) {
      bs[d] = 1;
    } else {
      bs[d] = (patch->ldims[d] == 1) ? 1 : 2;
    }
  }
  cell_map_init(&cuda->map, cuda->b_mx, bs);
}

static void
psc_particles_cuda_destroy(struct psc_particles *prts)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

  cell_map_free(&cuda->map);
}

#ifdef HAVE_LIBHDF5_HL

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------
// psc_particles_cuda_write

static void
psc_particles_cuda_write(struct psc_particles *prts, struct mrc_io *io)
{
  int ierr;
  assert(sizeof(particle_cuda_real_t) == sizeof(float));

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, prts), H5P_DEFAULT); H5_CHK(group);
  // save/restore n_alloced, too?
  ierr = H5LTset_attribute_int(group, ".", "p", &prts->p, 1); CE;
  ierr = H5LTset_attribute_int(group, ".", "n_part", &prts->n_part, 1); CE;
  ierr = H5LTset_attribute_uint(group, ".", "flags", &prts->flags, 1); CE;
  if (prts->n_part > 0) {
    float4 *xi4  = calloc(prts->n_part, sizeof(float4));
    float4 *pxi4 = calloc(prts->n_part, sizeof(float4));
  
    __particles_cuda_from_device(prts, xi4, pxi4);
  
    hsize_t hdims[2] = { prts->n_part, 4 };
    ierr = H5LTmake_dataset_float(group, "xi4", 2, hdims, (float *) xi4); CE;
    ierr = H5LTmake_dataset_float(group, "pxi4", 2, hdims, (float *) pxi4); CE;

    free(xi4);
    free(pxi4);
  }
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_particles_cuda_read

static void
psc_particles_cuda_read(struct psc_particles *prts, struct mrc_io *io)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, prts), H5P_DEFAULT); H5_CHK(group);
  ierr = H5LTget_attribute_int(group, ".", "p", &prts->p); CE;
  ierr = H5LTget_attribute_int(group, ".", "n_part", &prts->n_part); CE;
  ierr = H5LTget_attribute_uint(group, ".", "flags", &prts->flags); CE;
  psc_particles_setup(prts);
  ierr = H5Gclose(group); CE;
}

#endif

// FIXME, should go away and always be done within cuda for consistency

static inline int
find_cellIdx(struct psc_patch *patch, struct cell_map *map,
	     struct psc_particles *pp, int n)
{
  particle_c_t *p = particles_c_get_one(pp, n);
  particle_c_real_t dxi = 1.f / patch->dx[0];
  particle_c_real_t dyi = 1.f / patch->dx[1];
  particle_c_real_t dzi = 1.f / patch->dx[2];
  particle_c_real_t xi[3] = { p->xi * dxi, p->yi * dyi, p->zi * dzi };
  int pos[3];
  for (int d = 0; d < 3; d++) {
    pos[d] = particle_c_real_fint(xi[d]);
  }
  
  return cell_map_3to1(map, pos);
}

static inline int
find_blockIdx(struct psc_patch *patch, struct cell_map *map,
	      struct psc_particles *pp, int n, int blocksize[3])
{
  int cell_idx = find_cellIdx(patch, map, pp, n);
  return cell_idx / (blocksize[0] * blocksize[1] * blocksize[2]);
}

static inline void
blockIdx_to_blockCrd(struct psc_patch *patch, struct cell_map *map,
		     int bidx, int bi[3], int blocksize[3])
{
  int cidx = bidx * (blocksize[0] * blocksize[1] * blocksize[2]);
  cell_map_1to3(map, cidx, bi);
  for (int d = 0; d < 3; d++) {
    bi[d] /= blocksize[d];
  }
}

static void
copy_from(struct psc_particles *prts_cuda, struct psc_particles *prts,
	  void (*get_particle)(struct cuda_mparticles_prt *prt, int n, struct psc_particles *prts))
{
  struct psc_mparticles *mprts = psc_particles_cuda(prts_cuda)->mprts;
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int off = 0;
  for (int p = 0; p < prts->p; p++) {
    off += psc_mparticles_get_patch(mprts, p)->n_part;
  }

  cuda_mparticles_set_particles(cmprts, prts->n_part, off,
				(void (*)(struct cuda_mparticles_prt *, int, void *)) get_particle,
				prts);
}

// ======================================================================
// conversion to "c"

static inline void
calc_vxi(particle_c_real_t vxi[3], particle_c_t *part)
{
  particle_c_real_t root =
    1.f / particle_c_real_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

static void
get_particle_c(struct cuda_mparticles_prt *prt, int n, struct psc_particles *prts)
{
  particle_single_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  
  particle_c_t *part = particles_c_get_one(prts, n);

  particle_c_real_t vxi[3];
  calc_vxi(vxi, part);

  prt->xi[0]   = part->xi + dth[0] * vxi[0];
  prt->xi[1]   = part->yi + dth[1] * vxi[1];
  prt->xi[2]   = part->zi + dth[2] * vxi[2];
  prt->pxi[0]  = part->pxi;
  prt->pxi[1]  = part->pyi;
  prt->pxi[2]  = part->pzi;
  prt->kind    = part->kind;
  prt->qni_wni = part->qni * part->wni;
}

static void
psc_particles_cuda_copy_from_c(struct psc_particles *prts_cuda,
			       struct psc_particles *prts, unsigned int flags)
{
  assert(prts_cuda->n_part == prts->n_part);
  
  copy_from(prts_cuda, prts, get_particle_c);
}

static void
psc_particles_cuda_copy_to_c(struct psc_particles *prts_cuda,
			     struct psc_particles *prts_c, unsigned int flags)
{
  struct psc_particles_c *c = psc_particles_c(prts_c);
  prts_c->n_part = prts_cuda->n_part;
  assert(prts_c->n_part <= c->n_alloced);
  
  particle_single_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  
  float4 *xi4  = calloc(prts_cuda->n_part, sizeof(float4));
  float4 *pxi4 = calloc(prts_cuda->n_part, sizeof(float4));
  
  __particles_cuda_from_device(prts_cuda, xi4, pxi4);
  
  for (int n = 0; n < prts_c->n_part; n++) {
    particle_c_real_t qni_wni = pxi4[n].w;
    unsigned int kind = cuda_float_as_int(xi4[n].w);
    
    particle_c_t *part_base = particles_c_get_one(prts_c, n);
    part_base->xi  = xi4[n].x;
    part_base->yi  = xi4[n].y;
    part_base->zi  = xi4[n].z;
    part_base->pxi = pxi4[n].x;
    part_base->pyi = pxi4[n].y;
    part_base->pzi = pxi4[n].z;
    part_base->qni = ppsc->kinds[kind].q;
    part_base->mni = ppsc->kinds[kind].m;
    part_base->wni = qni_wni / part_base->qni;
    part_base->kind = kind;

    particle_c_real_t vxi[3];
    calc_vxi(vxi, part_base);
    part_base->xi -= dth[0] * vxi[0];
    part_base->yi -= dth[1] * vxi[1];
    part_base->zi -= dth[2] * vxi[2];
  }

  free(xi4);
  free(pxi4);
}

// ======================================================================
// conversion to "single"

static void
get_particle_single(struct cuda_mparticles_prt *prt, int n, struct psc_particles *prts)
{
  particle_single_t *part = particles_single_get_one(prts, n);

  prt->xi[0]   = part->xi;
  prt->xi[1]   = part->yi;
  prt->xi[2]   = part->zi;
  prt->pxi[0]  = part->pxi;
  prt->pxi[1]  = part->pyi;
  prt->pxi[2]  = part->pzi;
  prt->kind    = part->kind;
  prt->qni_wni = part->qni_wni;
}

static void
psc_particles_cuda_copy_from_single(struct psc_particles *prts_cuda,
				    struct psc_particles *prts, unsigned int flags)
{
  assert(prts_cuda->n_part == prts->n_part);
  
  copy_from(prts_cuda, prts, get_particle_single);
}

static void
psc_particles_cuda_copy_to_single(struct psc_particles *prts_cuda,
				  struct psc_particles *prts, unsigned int flags)
{
  int p = prts_cuda->p;
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts_cuda);
  struct psc_particles_single *sngl = psc_particles_single(prts);
  prts->n_part = prts_cuda->n_part;
  assert(prts->n_part <= sngl->n_alloced);
  
  float4 *xi4  = calloc(prts_cuda->n_part, sizeof(float4));
  float4 *pxi4 = calloc(prts_cuda->n_part, sizeof(float4));
  
  __particles_cuda_from_device(prts_cuda, xi4, pxi4);
  
  for (int n = 0; n < prts->n_part; n++) {
    particle_single_t *part_base = particles_single_get_one(prts, n);

    part_base->xi  = xi4[n].x;
    part_base->yi  = xi4[n].y;
    part_base->zi  = xi4[n].z;
    part_base->kind = cuda_float_as_int(xi4[n].w);
    part_base->pxi = pxi4[n].x;
    part_base->pyi = pxi4[n].y;
    part_base->pzi = pxi4[n].z;
    part_base->qni_wni = pxi4[n].w;

    for (int d = 0; d < 3; d++) {
      int bi = particle_single_real_fint((&part_base->xi)[d] * cuda->b_dxi[d]);
      if (bi < 0 || bi >= cuda->b_mx[d]) {
	MHERE;
	mprintf("XXX p %d xi %.10g %.10g %.10g\n", p, part_base->xi, part_base->yi, part_base->zi);
	mprintf("XXX p %d n %d d %d xi %.10g b_dxi %.10g bi %d // %d\n",
		p, n, d, (&part_base->xi)[d] * cuda->b_dxi[d], cuda->b_dxi[d], bi, cuda->b_mx[d]);
      }
    }

  }

  free(xi4);
  free(pxi4);
}

// ======================================================================
// conversion to "double"

static void
get_particle_double(struct cuda_mparticles_prt *prt, int n, struct psc_particles *prts)
{
  particle_double_t *part = particles_double_get_one(prts, n);

  prt->xi[0]   = part->xi;
  prt->xi[1]   = part->yi;
  prt->xi[2]   = part->zi;
  prt->pxi[0]  = part->pxi;
  prt->pxi[1]  = part->pyi;
  prt->pxi[2]  = part->pzi;
  prt->kind    = part->kind;
  prt->qni_wni = part->qni_wni;
}

static void
psc_particles_cuda_copy_from_double(struct psc_particles *prts_cuda,
				    struct psc_particles *prts, unsigned int flags)
{
  assert(prts_cuda->n_part == prts->n_part);
  
  copy_from(prts_cuda, prts, get_particle_double);
}

static void
psc_particles_cuda_copy_to_double(struct psc_particles *prts_cuda,
				  struct psc_particles *prts, unsigned int flags)
{
  struct psc_particles_double *sngl = psc_particles_double(prts);
  prts->n_part = prts_cuda->n_part;
  assert(prts->n_part <= sngl->n_alloced);
  
  float4 *xi4  = calloc(prts_cuda->n_part, sizeof(float4));
  float4 *pxi4 = calloc(prts_cuda->n_part, sizeof(float4));
  
  __particles_cuda_from_device(prts_cuda, xi4, pxi4);
  
  for (int n = 0; n < prts->n_part; n++) {
    particle_double_t *part_base = particles_double_get_one(prts, n);

    part_base->xi  = xi4[n].x;
    part_base->yi  = xi4[n].y;
    part_base->zi  = xi4[n].z;
    part_base->kind = cuda_float_as_int(xi4[n].w);
    part_base->pxi = pxi4[n].x;
    part_base->pyi = pxi4[n].y;
    part_base->pzi = pxi4[n].z;
    part_base->qni_wni = pxi4[n].w;
  }

  free(xi4);
  free(pxi4);
}

// ======================================================================
// psc_particles: subclass "cuda"

static struct mrc_obj_method psc_particles_cuda_methods[] = {
  MRC_OBJ_METHOD("copy_to_c"       , psc_particles_cuda_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c"     , psc_particles_cuda_copy_from_c),
  MRC_OBJ_METHOD("copy_to_single"  , psc_particles_cuda_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_particles_cuda_copy_from_single),
  MRC_OBJ_METHOD("copy_to_double"  , psc_particles_cuda_copy_to_double),
  MRC_OBJ_METHOD("copy_from_double", psc_particles_cuda_copy_from_double),
  {}
};

struct psc_particles_ops psc_particles_cuda_ops = {
  .name                    = "cuda",
  .size                    = sizeof(struct psc_particles_cuda),
  .methods                 = psc_particles_cuda_methods,
  .setup                   = psc_particles_cuda_setup,
  .destroy                 = psc_particles_cuda_destroy,
#ifdef HAVE_LIBHDF5_HL
  .read                    = psc_particles_cuda_read,
  .write                   = psc_particles_cuda_write,
#endif
};

// ======================================================================
// psc_mparticles "cuda"

// ----------------------------------------------------------------------
// psc_mparticles_cuda_setup

static void
psc_mparticles_cuda_setup(struct psc_mparticles *mprts)
{
  psc_mparticles_setup_super(mprts);
  __psc_mparticles_cuda_setup(mprts);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_destroy

static void
psc_mparticles_cuda_destroy(struct psc_mparticles *mprts)
{
  __psc_mparticles_cuda_free(mprts);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_read

static void
psc_mparticles_cuda_read(struct psc_mparticles *mprts, struct mrc_io *io)
{
  mprts->domain = mrc_io_read_ref(io, mprts, "domain", mrc_domain);
  mrc_domain_get_patches(mprts->domain, &mprts->nr_patches);
  mrc_io_read_int(io, mprts, "flags", (int *) &mprts->flags);
  
  mprts->nr_particles_by_patch =
    calloc(mprts->nr_patches, sizeof(*mprts->nr_particles_by_patch));
  mprts->prts = calloc(mprts->nr_patches, sizeof(*mprts->prts));
  for (int p = 0; p < mprts->nr_patches; p++) {
    char name[20]; sprintf(name, "prts%d", p);
    mprts->prts[p] = mrc_io_read_ref(io, mprts, name, psc_particles); // doesn't actually read data
    mprts->nr_particles_by_patch[p] = mprts->prts[p]->n_part;
  }

  psc_mparticles_setup(mprts);

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = mprts->prts[p];
    if (prts->n_part > 0) {
      int ierr;
      char name[20]; sprintf(name, "prts%d", p);
      char *path;
      mrc_io_read_attr_string(io, mrc_io_obj_path(io, mprts), name, &path);
      hid_t group = H5Gopen(h5_file, path, H5P_DEFAULT); H5_CHK(group);
      free(path);
      float4 *xi4  = calloc(prts->n_part, sizeof(float4));
      float4 *pxi4 = calloc(prts->n_part, sizeof(float4));
      
      ierr = H5LTread_dataset_float(group, "xi4", (float *) xi4); CE;
      ierr = H5LTread_dataset_float(group, "pxi4", (float *) pxi4); CE;
      
      struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
      struct cuda_mparticles *cmprts = mprts_cuda->cmprts;
      
      cuda_mparticles_to_device(cmprts, xi4, pxi4, prts->n_part, off);
      
      free(xi4);
      free(pxi4);
      ierr = H5Gclose(group); CE;
    }
    off += prts->n_part;
  }

  psc_mparticles_setup_internals(mprts);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_setup_internals

static void
psc_mparticles_cuda_setup_internals(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  assert((mprts->flags & MP_NEED_BLOCK_OFFSETS) &&
	 !(mprts->flags & MP_NEED_CELL_OFFSETS));

  unsigned int n_prts_by_patch[mprts->nr_patches];
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    n_prts_by_patch[p] = prts->n_part;
  }

  cuda_mparticles_sort_initial(cmprts, n_prts_by_patch);
}

// ======================================================================
// psc_mparticles: subclass "cuda"
  
struct psc_mparticles_ops psc_mparticles_cuda_ops = {
  .name                    = "cuda",
  .size                    = sizeof(struct psc_mparticles_cuda),
  .setup                   = psc_mparticles_cuda_setup,
  .destroy                 = psc_mparticles_cuda_destroy,
  .read                    = psc_mparticles_cuda_read,
  //  .write                   = psc_mparticles_cuda_write,
  .setup_internals         = psc_mparticles_cuda_setup_internals,
};

void
psc_mparticles_cuda_reorder(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  if (mprts_cuda->need_reorder) {
    cuda_mprts_reorder(mprts);
    mprts_cuda->need_reorder = false;
  }
}
