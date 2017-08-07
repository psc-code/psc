
#include "psc.h"
#include "cuda_iface.h"
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

// FIXME, should go away and always be done within cuda for consistency

#if 0
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
#endif

static void
copy_from(int p, struct psc_mparticles *mprts_cuda, struct psc_mparticles *mprts,
	  void (*get_particle)(struct cuda_mparticles_prt *prt, int n, void *ctx))
{
  struct psc_particles *prts_cuda = psc_mparticles_get_patch(mprts_cuda, p);
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts_cuda)->cmprts;

  psc_particles_set_n_prts(prts_cuda, psc_particles_size(prts));

  unsigned int off = 0;
  for (int pp = 0; pp < p; pp++) {
    off += psc_mparticles_n_prts_by_patch(mprts_cuda, pp);
  }

  cuda_mparticles_set_particles(cmprts, psc_mparticles_n_prts_by_patch(mprts, p), off,
				get_particle, prts);
}

static void
copy_to(int p, struct psc_mparticles *mprts_cuda, struct psc_mparticles *mprts,
	void (*put_particle)(struct cuda_mparticles_prt *prt, int n, void *ctx))
{
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts_cuda)->cmprts;

  psc_mparticles_resize_patch(mprts, p, psc_mparticles_n_prts_by_patch(mprts_cuda, p));

  unsigned int off = 0;
  for (int pp = 0; pp < p; pp++) {
    off += psc_mparticles_n_prts_by_patch(mprts_cuda, pp);
  }

  cuda_mparticles_get_particles(cmprts, psc_mparticles_n_prts_by_patch(mprts_cuda, p), off,
				put_particle, prts);
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
get_particle_c(struct cuda_mparticles_prt *prt, int n, void *ctx)
{
  struct psc_particles *prts = ctx;
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
put_particle_c(struct cuda_mparticles_prt *prt, int n, void *ctx)
{
  struct psc_particles *prts = ctx;
  particle_single_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  
  particle_c_real_t qni_wni = prt->qni_wni;
  unsigned int kind = prt->kind;
  
  particle_c_t *part = particles_c_get_one(prts, n);
  part->xi  = prt->xi[0];
  part->yi  = prt->xi[1];
  part->zi  = prt->xi[2];
  part->pxi = prt->pxi[0];
  part->pyi = prt->pxi[1];
  part->pzi = prt->pxi[2];
  part->qni = ppsc->kinds[kind].q;
  part->mni = ppsc->kinds[kind].m;
  part->wni = qni_wni / part->qni;
  part->kind = kind;
  
  particle_c_real_t vxi[3];
  calc_vxi(vxi, part);
  part->xi -= dth[0] * vxi[0];
  part->yi -= dth[1] * vxi[1];
  part->zi -= dth[2] * vxi[2];
}

static void
psc_particles_cuda_copy_from_c(int p, struct psc_mparticles *mprts_cuda,
			       struct psc_mparticles *mprts, unsigned int flags)
{
  copy_from(p, mprts_cuda, mprts, get_particle_c);
}

static void
psc_particles_cuda_copy_to_c(int p, struct psc_mparticles *mprts_cuda,
			     struct psc_mparticles *mprts, unsigned int flags)
{
  copy_to(p, mprts_cuda, mprts, put_particle_c);
}

// ======================================================================
// conversion to "single"

static void
get_particle_single(struct cuda_mparticles_prt *prt, int n, void *ctx)
{
  struct psc_particles *prts = ctx;
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
put_particle_single(struct cuda_mparticles_prt *prt, int n, void *ctx)
{
  struct psc_particles *prts = ctx;
  particle_single_t *part = particles_single_get_one(prts, n);
  
  part->xi      = prt->xi[0];
  part->yi      = prt->xi[1];
  part->zi      = prt->xi[2];
  part->kind    = prt->kind;
  part->pxi     = prt->pxi[0];
  part->pyi     = prt->pxi[1];
  part->pzi     = prt->pxi[2];
  part->qni_wni = prt->qni_wni;
}

static void
psc_particles_cuda_copy_from_single(int p, struct psc_mparticles *mprts_cuda,
				    struct psc_mparticles *mprts, unsigned int flags)
{
  copy_from(p, mprts_cuda, mprts, get_particle_single);
}

static void
psc_particles_cuda_copy_to_single(int p, struct psc_mparticles *mprts_cuda,
				  struct psc_mparticles *mprts, unsigned int flags)
{
  copy_to(p, mprts_cuda, mprts, put_particle_single);
}

// ======================================================================
// conversion to "double"

static void
get_particle_double(struct cuda_mparticles_prt *prt, int n, void *ctx)
{
  struct psc_particles *prts = ctx;
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
put_particle_double(struct cuda_mparticles_prt *prt, int n, void *ctx)
{
  struct psc_particles *prts = ctx;
  particle_double_t *part = particles_double_get_one(prts, n);
  
  part->xi      = prt->xi[0];
  part->yi      = prt->xi[1];
  part->zi      = prt->xi[2];
  part->kind    = prt->kind;
  part->pxi     = prt->pxi[0];
  part->pyi     = prt->pxi[1];
  part->pzi     = prt->pxi[2];
  part->qni_wni = prt->qni_wni;
}

static void
psc_particles_cuda_copy_from_double(int p, struct psc_mparticles *mprts_cuda,
				    struct psc_mparticles *mprts, unsigned int flags)
{
  copy_from(p, mprts_cuda, mprts, get_particle_double);
}

static void
psc_particles_cuda_copy_to_double(int p, struct psc_mparticles *mprts_cuda,
				  struct psc_mparticles *mprts, unsigned int flags)
{
  copy_to(p, mprts_cuda, mprts, put_particle_double);
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
};

// ======================================================================
// psc_mparticles "cuda"

// ----------------------------------------------------------------------
// psc_mparticles_cuda_setup

static void
psc_mparticles_cuda_setup(struct psc_mparticles *mprts)
{
  psc_mparticles_setup_super(mprts);

  cuda_base_init();

  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = cuda_mparticles_create();
  mprts_cuda->cmprts = cmprts;

  if (mprts->nr_patches == 0) {
    return;
  }
  
  if (!mprts->flags) {
    // FIXME, they get set too late, so auto-dispatch "1vb" doesn't work
    mprts->flags = MP_NEED_BLOCK_OFFSETS | MP_BLOCKSIZE_4X4X4 | MP_NO_CHECKERBOARD;
  }

  int *ldims = ppsc->patch[0].ldims;
  double *dx = ppsc->patch[0].dx;

  // check that all patches≈ have equal size
  for (int p = 1; p < mprts->nr_patches; p++) {
    for (int d = 0; d < 3; d++) {
      assert(ppsc->patch[p].ldims[d] == ldims[d]);
      assert(ppsc->patch[p].dx[d] == dx[d]);
    }
  }

  struct cuda_domain_info domain_info = {};

  domain_info.n_patches = mprts->nr_patches;
  domain_info.xb_by_patch = calloc(domain_info.n_patches, sizeof(*domain_info.xb_by_patch));

  int bs[3];
  for (int d = 0; d < 3; d++) {
    switch (mprts->flags & MP_BLOCKSIZE_MASK) {
    case MP_BLOCKSIZE_1X1X1: bs[d] = 1; break;
    case MP_BLOCKSIZE_2X2X2: bs[d] = 2; break;
    case MP_BLOCKSIZE_4X4X4: bs[d] = 4; break;
    case MP_BLOCKSIZE_8X8X8: bs[d] = 8; break;
    default: assert(0);
    }
    if (ppsc->domain.gdims[d] == 1) {
      bs[d] = 1;
    }
    assert(ldims[d] % bs[d] == 0); // FIXME not sure what breaks if not
    mprts_cuda->b_mx[d]  = ldims[d] / bs[d];
    mprts_cuda->b_dxi[d] = 1.f / (bs[d] * dx[d]);
    domain_info.ldims[d] = ldims[d];
    domain_info.bs[d]    = bs[d];
    domain_info.dx[d]    = dx[d];

    for (int p = 0; p < domain_info.n_patches; p++) {
      domain_info.xb_by_patch[p][0] = ppsc->patch[p].xb[0];
      domain_info.xb_by_patch[p][1] = ppsc->patch[p].xb[1];
      domain_info.xb_by_patch[p][2] = ppsc->patch[p].xb[2];
    }
  }

  cuda_mparticles_set_domain_info(cmprts, &domain_info);

  free(domain_info.xb_by_patch);

  unsigned int n_prts_by_patch[cmprts->n_patches];
  for (int p = 0; p < cmprts->n_patches; p++) {
    n_prts_by_patch[p] = psc_mparticles_n_prts_by_patch(mprts, p);
  }

  cuda_mparticles_alloc(cmprts, n_prts_by_patch);

  __psc_mparticles_cuda_setup(mprts);

  mprts_cuda->bnd = calloc(mprts->nr_patches, sizeof(*mprts_cuda->bnd));
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_destroy

static void
psc_mparticles_cuda_destroy(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;
  assert(cmprts);
  
  __psc_mparticles_cuda_free(mprts);

  cuda_mparticles_dealloc(cmprts);
  cuda_mparticles_destroy(cmprts);
  mprts_cuda->cmprts = NULL;
}

// ----------------------------------------------------------------------
// psc_mparticls_cuda_reorder

void
psc_mparticles_cuda_reorder(struct psc_mparticles *mprts)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  if (cmprts->need_reorder) {
    cuda_mprts_reorder(mprts);
    cmprts->need_reorder = false;
  }
}
#ifdef HAVE_LIBHDF5_HL

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------
// psc_mparticles_cuda_write

static void
psc_mparticles_cuda_write(struct psc_mparticles *mprts, struct mrc_io *io)
{
  int ierr;
  
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);
  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gcreate(group, pname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);
    ierr = H5LTset_attribute_int(pgroup, ".", "n_prts", &n_prts, 1); CE;
    if (n_prts > 0) {
      float4 *xi4  = calloc(n_prts, sizeof(float4));
      float4 *pxi4 = calloc(n_prts, sizeof(float4));
      
      __particles_cuda_from_device(mprts, xi4, pxi4, off, n_prts);
      
      hsize_t hdims[2] = { n_prts, 4 };
      ierr = H5LTmake_dataset_float(pgroup, "xi4", 2, hdims, (float *) xi4); CE;
      ierr = H5LTmake_dataset_float(pgroup, "pxi4", 2, hdims, (float *) pxi4); CE;
      
      free(xi4);
      free(pxi4);
    }
    ierr = H5Gclose(pgroup); CE;
    off += n_prts;
  }
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_read

static void
psc_mparticles_cuda_read(struct psc_mparticles *mprts, struct mrc_io *io)
{
  mprts->domain = mrc_io_read_ref(io, mprts, "domain", mrc_domain);
  mrc_io_read_int(io, mprts, "flags", (int *) &mprts->flags);
  mrc_domain_get_patches(mprts->domain, &mprts->nr_patches);

  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);
  mprts->nr_particles_by_patch =
    calloc(mprts->nr_patches, sizeof(*mprts->nr_particles_by_patch));
  mprts->prts = calloc(mprts->nr_patches, sizeof(*mprts->prts));

  for (int p = 0; p < mprts->nr_patches; p++) {
    char name[20]; sprintf(name, "prts%d", p);
    mprts->prts[p] = psc_particles_create(MPI_COMM_NULL);
    psc_particles_set_type(mprts->prts[p], "cuda");
    psc_particles_set_name(mprts->prts[p], name);
    mprts->prts[p]->p = p;

    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts;
    ierr = H5LTget_attribute_int(pgroup, ".", "n_prts", &n_prts); CE;
    mprts->nr_particles_by_patch[p] = n_prts;

    ierr = H5Gclose(pgroup); CE;
  }

  psc_mparticles_setup(mprts);

  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);
    if (n_prts > 0) {
      float4 *xi4  = calloc(n_prts, sizeof(float4));
      float4 *pxi4 = calloc(n_prts, sizeof(float4));
      
      ierr = H5LTread_dataset_float(pgroup, "xi4", (float *) xi4); CE;
      ierr = H5LTread_dataset_float(pgroup, "pxi4", (float *) pxi4); CE;
      
      struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
      struct cuda_mparticles *cmprts = mprts_cuda->cmprts;
      
      cuda_mparticles_to_device(cmprts, xi4, pxi4, n_prts, off);
      
      free(xi4);
      free(pxi4);
    }

    ierr = H5Gclose(pgroup); CE;
    off += n_prts;
  }

  ierr = H5Gclose(group); CE;
  psc_mparticles_setup_internals(mprts);
}

#endif

// ----------------------------------------------------------------------
// psc_mparticles_cuda_update_n_part

static void
psc_mparticles_cuda_update_n_part(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int n_prts_by_patch[cmprts->n_patches];
  cuda_mparticles_get_n_prts_by_patch(cmprts, n_prts_by_patch);

  unsigned int n_prts = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    psc_particles_set_n_prts(prts, n_prts_by_patch[p]);
    n_prts += n_prts_by_patch[p];
  }
  assert(cmprts->n_prts == n_prts);
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

  unsigned int n_prts = 0;
  unsigned int n_prts_by_patch[mprts->nr_patches];
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    n_prts_by_patch[p] = psc_particles_size(prts);
    n_prts += n_prts_by_patch[p];
  }
  cmprts->n_prts = n_prts; // another hack to deal with copy_from having
  // updated n_part for the patches

  cuda_mparticles_sort_initial(cmprts, n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_get_nr_particles

static unsigned int
psc_mparticles_cuda_get_nr_particles(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  return cmprts->n_prts;
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_inject

#include <psc_particles_as_single.h> // FIXME

void
psc_mparticles_cuda_inject(struct psc_mparticles *mprts_base, struct cuda_mparticles_prt *buf,
			   unsigned int *buf_n_by_patch)
{
  assert(strcmp(psc_mparticles_type(mprts_base), "cuda") == 0);
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts_base);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

#if 1
  cuda_mparticles_inject(cmprts, buf, buf_n_by_patch);
#else
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);

  unsigned buf_n = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    
    int i = psc_particles_size(prts);
    particles_realloc(prts, i + buf_n_by_patch[p]);
    for (int cnt = 0; cnt < buf_n_by_patch[p]; cnt++) {
      particle_t *prt = particles_get_one(prts, i++);
      
      prt->xi = buf[buf_n + cnt].xi[0];
      prt->yi = buf[buf_n + cnt].xi[1];
      prt->zi = buf[buf_n + cnt].xi[2];
      prt->pxi = buf[buf_n + cnt].pxi[0];
      prt->pyi = buf[buf_n + cnt].pxi[1];
      prt->pzi = buf[buf_n + cnt].pxi[2];
      prt->kind = buf[buf_n + cnt].kind;
      prt->qni_wni = buf[buf_n + cnt].qni_wni;
    }
    buf_n += buf_n_by_patch[p];
    psc_mparticles_resize_patch(mprts, p, i);
  }

  psc_mparticles_put_as(mprts, mprts_base, 0);
#endif
}

// ----------------------------------------------------------------------
// psc_mparticles: subclass "cuda"
  
struct psc_mparticles_ops psc_mparticles_cuda_ops = {
  .name                    = "cuda",
  .size                    = sizeof(struct psc_mparticles_cuda),
  .methods                 = psc_particles_cuda_methods,
  .setup                   = psc_mparticles_cuda_setup,
  .destroy                 = psc_mparticles_cuda_destroy,
  .read                    = psc_mparticles_cuda_read,
  .write                   = psc_mparticles_cuda_write,
  .update_n_part           = psc_mparticles_cuda_update_n_part,
  .setup_internals         = psc_mparticles_cuda_setup_internals,
  .get_nr_particles        = psc_mparticles_cuda_get_nr_particles,
};

