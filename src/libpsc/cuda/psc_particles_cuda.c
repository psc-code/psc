
#include "psc.h"
#include "cuda_iface.h"
#include "psc_cuda.h"
#include "psc_bnd_cuda.h"
#include "psc_particles_cuda.h"
#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "psc_push_particles.h"
#include "particles_cuda.h"

#include "json-builder.h"

#include <mrc_io.h>

EXTERN_C void cuda_init(int rank);

// ======================================================================
// psc_mparticles "cuda"

// FIXME, should go away and always be done within cuda for consistency

#if 0
static inline int
find_cellIdx(struct psc_patch *patch, struct cell_map *map,
	     struct psc_mparticles *mprts, int p, int n)
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
	      struct psc_mparticles *mprts, int p, int n, int blocksize[3])
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

struct copy_ctx {
  struct psc_mparticles *mprts;
  int p;
};

static void
copy_from(struct psc_mparticles *mprts, struct psc_mparticles *mprts_from,
	  void (*get_particle)(struct cuda_mparticles_prt *prt, int n, void *ctx))
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  unsigned int n_prts_by_patch[mprts->nr_patches];
  cuda_mparticles_get_size_all(cmprts, n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    int n_prts = n_prts_by_patch[p];
    struct copy_ctx ctx = { .mprts = mprts_from, .p = p };
    cuda_mparticles_set_particles(cmprts, n_prts, off, get_particle, &ctx);

    off += n_prts;
  }
}

static void
copy_to(struct psc_mparticles *mprts, struct psc_mparticles *mprts_to,
	void (*put_particle)(struct cuda_mparticles_prt *prt, int n, void *ctx))
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;
  
  unsigned int n_prts_by_patch[mprts->nr_patches];
  cuda_mparticles_get_size_all(cmprts, n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    int n_prts = n_prts_by_patch[p];
    struct copy_ctx ctx = { .mprts = mprts_to, .p = p };
    cuda_mparticles_get_particles(cmprts, n_prts, off, put_particle, &ctx);

    off += n_prts;
  }
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
get_particle_c(struct cuda_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy_ctx *ctx = _ctx;
  particle_c_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  
  particle_c_t *prt_c = psc_mparticles_c_get_one(ctx->mprts, ctx->p, n);

  particle_c_real_t vxi[3];
  calc_vxi(vxi, prt_c);

  prt->xi[0]   = prt_c->xi + dth[0] * vxi[0];
  prt->xi[1]   = prt_c->yi + dth[1] * vxi[1];
  prt->xi[2]   = prt_c->zi + dth[2] * vxi[2];
  prt->pxi[0]  = prt_c->pxi;
  prt->pxi[1]  = prt_c->pyi;
  prt->pxi[2]  = prt_c->pzi;
  prt->kind    = prt_c->kind;
  prt->qni_wni = prt_c->qni * prt_c->wni;
}

static void
put_particle_c(struct cuda_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy_ctx *ctx = _ctx;
  particle_c_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  
  particle_c_real_t qni_wni = prt->qni_wni;
  unsigned int kind = prt->kind;
  
  particle_c_t *prt_c = psc_mparticles_c_get_one(ctx->mprts, ctx->p, n);
  prt_c->xi  = prt->xi[0];
  prt_c->yi  = prt->xi[1];
  prt_c->zi  = prt->xi[2];
  prt_c->pxi = prt->pxi[0];
  prt_c->pyi = prt->pxi[1];
  prt_c->pzi = prt->pxi[2];
  prt_c->qni = ppsc->kinds[kind].q;
  prt_c->mni = ppsc->kinds[kind].m;
  prt_c->wni = qni_wni / prt_c->qni;
  prt_c->kind = kind;
  
  particle_c_real_t vxi[3];
  calc_vxi(vxi, prt_c);
  prt_c->xi -= dth[0] * vxi[0];
  prt_c->yi -= dth[1] * vxi[1];
  prt_c->zi -= dth[2] * vxi[2];
}

static void
psc_mparticles_cuda_copy_from_c(struct psc_mparticles *mprts_cuda,
			       struct psc_mparticles *mprts, unsigned int flags)
{
  copy_from(mprts_cuda, mprts, get_particle_c);
}

static void
psc_mparticles_cuda_copy_to_c(struct psc_mparticles *mprts_cuda,
			     struct psc_mparticles *mprts, unsigned int flags)
{
  copy_to(mprts_cuda, mprts, put_particle_c);
}

// ======================================================================
// conversion to "single"

static void
get_particle_single(struct cuda_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy_ctx *ctx = _ctx;
  particle_single_t *part = psc_mparticles_single_get_one(ctx->mprts, ctx->p, n);

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
put_particle_single(struct cuda_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy_ctx *ctx = _ctx;
  particle_single_t *part = psc_mparticles_single_get_one(ctx->mprts, ctx->p, n);
  
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
psc_mparticles_cuda_copy_from_single(struct psc_mparticles *mprts_cuda,
				    struct psc_mparticles *mprts, unsigned int flags)
{
  copy_from(mprts_cuda, mprts, get_particle_single);
}

static void
psc_mparticles_cuda_copy_to_single(struct psc_mparticles *mprts_cuda,
				  struct psc_mparticles *mprts, unsigned int flags)
{
  copy_to(mprts_cuda, mprts, put_particle_single);
}

// ======================================================================
// conversion to "double"

static void
get_particle_double(struct cuda_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy_ctx *ctx = _ctx;
  particle_double_t *part = psc_mparticles_double_get_one(ctx->mprts, ctx->p, n);

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
put_particle_double(struct cuda_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy_ctx *ctx = _ctx;
  particle_double_t *part = psc_mparticles_double_get_one(ctx->mprts, ctx->p, n);
  
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
psc_mparticles_cuda_copy_from_double(struct psc_mparticles *mprts_cuda,
				    struct psc_mparticles *mprts, unsigned int flags)
{
  copy_from(mprts_cuda, mprts, get_particle_double);
}

static void
psc_mparticles_cuda_copy_to_double(struct psc_mparticles *mprts_cuda,
				  struct psc_mparticles *mprts, unsigned int flags)
{
  copy_to(mprts_cuda, mprts, put_particle_double);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_methods

static struct mrc_obj_method psc_mparticles_cuda_methods[] = {
  MRC_OBJ_METHOD("copy_to_c"       , psc_mparticles_cuda_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c"     , psc_mparticles_cuda_copy_from_c),
  MRC_OBJ_METHOD("copy_to_single"  , psc_mparticles_cuda_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_mparticles_cuda_copy_from_single),
  MRC_OBJ_METHOD("copy_to_double"  , psc_mparticles_cuda_copy_to_double),
  MRC_OBJ_METHOD("copy_from_double", psc_mparticles_cuda_copy_from_double),
  {}
};

// ----------------------------------------------------------------------
// psc_mparticles_cuda_setup

static void
psc_mparticles_cuda_setup(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  cuda_base_init();

  struct cuda_mparticles *cmprts = cuda_mparticles_create();
  mprts_cuda->cmprts = cmprts;

  psc_mparticles_setup_super(mprts);

  int n_patches = mprts->nr_patches;
  assert(n_patches != 0);
  
  if (!mprts->flags) {
    // FIXME, they get set too late, so auto-dispatch "1vb" doesn't work
    mprts->flags = MP_NEED_BLOCK_OFFSETS | MP_BLOCKSIZE_4X4X4 | MP_NO_CHECKERBOARD;
  }

  int *ldims = ppsc->patch[0].ldims;
  double *dx = ppsc->patch[0].dx;

    // check that all patches have equal size
  for (int p = 1; p < n_patches; p++) {
    for (int d = 0; d < 3; d++) {
      assert(ppsc->patch[p].ldims[d] == ldims[d]);
      assert(ppsc->patch[p].dx[d] == dx[d]);
    }
  }

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
  }

  json_value *info = json_object_new(0);
  
  json_object_push(info, "n_patches", json_integer_new(n_patches));
  
  json_value *arr_ldims = json_array_new(3);
  json_object_push(info, "ldims", arr_ldims);
  
  json_value *arr_bs = json_array_new(3);
  json_object_push(info, "bs", arr_bs);
  
  json_value *arr_dx = json_array_new(3);
  json_object_push(info, "dx", arr_dx);
  
  for (int d = 0; d < 3; d++) {
    json_array_push(arr_ldims, json_integer_new(ldims[d]));
    json_array_push(arr_bs, json_integer_new(bs[d]));
    json_array_push(arr_dx, json_double_new(dx[d]));
  }
  
  json_value *arr_xb_by_patch = json_array_new(n_patches);
  json_object_push(info, "xb_by_patch", arr_xb_by_patch);
  for (int p = 0; p < n_patches; p++) {
    json_value *arr_xb = json_array_new(3);
    json_array_push(arr_xb_by_patch, arr_xb);
    for (int d = 0; d < 3; d++) {
      json_array_push(arr_xb, json_double_new(ppsc->patch[p].xb[d]));
    }
  }

  json_value *obj = json_object_new(0);
  json_object_push(obj, "info", info);

  cuda_mparticles_ctor(cmprts, mrc_json_from_json_parser(obj));

  json_builder_free(obj);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_destroy

static void
psc_mparticles_cuda_destroy(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;
  assert(cmprts);
  
  cuda_mparticles_destroy(cmprts);
  mprts_cuda->cmprts = NULL;
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_reserve_all

static void
psc_mparticles_cuda_reserve_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  cuda_mparticles_reserve_all(cmprts, (unsigned int *) n_prts_by_patch);
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
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;
  int ierr;
  
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  unsigned int n_prts_by_patch[mprts->nr_patches];
  cuda_mparticles_get_size_all(cmprts, n_prts_by_patch);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);
  unsigned int off = 0;
  // FIXME, reorder first if necessary
  for (int p = 0; p < mprts->nr_patches; p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gcreate(group, pname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = n_prts_by_patch[p];
    ierr = H5LTset_attribute_int(pgroup, ".", "n_prts", &n_prts, 1); CE;
    if (n_prts > 0) {
      float_4 *xi4  = calloc(n_prts, sizeof(*xi4));
      float_4 *pxi4 = calloc(n_prts, sizeof(*pxi4));
      
      cuda_mparticles_from_device(cmprts, xi4, pxi4, n_prts, off);
      
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
  psc_mparticles_read_super(mprts, io);

  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);

  int n_prts_by_patch[mprts->nr_patches];

  for (int p = 0; p < mprts->nr_patches; p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts;
    ierr = H5LTget_attribute_int(pgroup, ".", "n_prts", &n_prts); CE;
    n_prts_by_patch[p] = n_prts;
    ierr = H5Gclose(pgroup); CE;
  }

  psc_mparticles_setup(mprts);
  psc_mparticles_reserve_all(mprts, n_prts_by_patch);

  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = n_prts_by_patch[p];
    if (n_prts > 0) {
      float_4 *xi4  = calloc(n_prts, sizeof(float4));
      float_4 *pxi4 = calloc(n_prts, sizeof(float4));
      
      ierr = H5LTread_dataset_float(pgroup, "xi4", (float *) xi4); CE;
      ierr = H5LTread_dataset_float(pgroup, "pxi4", (float *) pxi4); CE;
      
      struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;
      
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
// psc_mparticles_cuda_setup_internals

static void
psc_mparticles_cuda_setup_internals(struct psc_mparticles *mprts)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  cuda_mparticles_setup_internals(cmprts);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_get_nr_particles

static unsigned int
psc_mparticles_cuda_get_nr_particles(struct psc_mparticles *mprts)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  return cuda_mparticles_get_n_prts(cmprts);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_get_size_all

static void
psc_mparticles_cuda_get_size_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  cuda_mparticles_get_size_all(cmprts, (unsigned int *) n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_resize_all

static void
psc_mparticles_cuda_resize_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  cuda_mparticles_resize_all(cmprts, (unsigned int *) n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_inject

void
psc_mparticles_cuda_inject(struct psc_mparticles *mprts, struct cuda_mparticles_prt *buf,
			   unsigned int *buf_n_by_patch)
{
  assert(strcmp(psc_mparticles_type(mprts), "cuda") == 0);
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  cuda_mparticles_inject(cmprts, buf, buf_n_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_patch_get_b_dxi

const particle_cuda_real_t *
psc_mparticles_cuda_patch_get_b_dxi(struct psc_mparticles *mprts, int p)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  return cuda_mparticles_patch_get_b_dxi(cmprts, p);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_patch_get_b_mx

const int *
psc_mparticles_cuda_patch_get_b_mx(struct psc_mparticles *mprts, int p)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  return cuda_mparticles_patch_get_b_mx(cmprts, p);
}

// ----------------------------------------------------------------------
// psc_mparticles: subclass "cuda"
  
struct psc_mparticles_ops psc_mparticles_cuda_ops = {
  .name                    = "cuda",
  .size                    = sizeof(struct psc_mparticles_cuda),
  .methods                 = psc_mparticles_cuda_methods,
  .setup                   = psc_mparticles_cuda_setup,
  .destroy                 = psc_mparticles_cuda_destroy,
  .read                    = psc_mparticles_cuda_read,
  .write                   = psc_mparticles_cuda_write,
  .setup_internals         = psc_mparticles_cuda_setup_internals,
  .reserve_all             = psc_mparticles_cuda_reserve_all,
  .get_nr_particles        = psc_mparticles_cuda_get_nr_particles,
  .resize_all              = psc_mparticles_cuda_resize_all,
  .get_size_all            = psc_mparticles_cuda_get_size_all,
};

