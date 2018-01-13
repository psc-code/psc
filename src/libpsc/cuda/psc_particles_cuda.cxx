
#include "psc.h"
#include "cuda_iface.h"
#include "psc_particles_cuda.h"
#include "psc_particles_single.h"
#include "psc_particles_double.h"
#include "psc_push_particles.h"

#include <mrc_io.h>

// ======================================================================
// psc_mparticles "cuda"

// FIXME, should go away and always be done within cuda for consistency

struct copy_ctx {
  struct psc_mparticles *mprts;
  int p;
};

static void
copy_from(mparticles_cuda_t mprts, struct psc_mparticles *mprts_from,
	  void (*get_particle)(struct cuda_mparticles_prt *prt, int n, void *ctx))
{
  struct cuda_mparticles *cmprts = mprts.sub()->cmprts;

  unsigned int n_prts_by_patch[mprts.n_patches()];
  cmprts->get_size_all(n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    int n_prts = n_prts_by_patch[p];
    struct copy_ctx ctx = { .mprts = mprts_from, .p = p };
    cmprts->set_particles(n_prts, off, get_particle, &ctx);

    off += n_prts;
  }
}

static void
copy_to(mparticles_cuda_t mprts, struct psc_mparticles *mprts_to,
	void (*put_particle)(struct cuda_mparticles_prt *prt, int n, void *ctx))
{
  struct cuda_mparticles *cmprts = mprts.sub()->cmprts;
  
  unsigned int n_prts_by_patch[mprts.n_patches()];
  cmprts->get_size_all(n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    int n_prts = n_prts_by_patch[p];
    struct copy_ctx ctx = { .mprts = mprts_to, .p = p };
    cmprts->get_particles(n_prts, off, put_particle, &ctx);

    off += n_prts;
  }
}

// ======================================================================
// conversion to "single"

static void
get_particle_single(struct cuda_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy_ctx *ctx = (struct copy_ctx *) _ctx;
  particle_single_t *part = &mparticles_single_t(ctx->mprts)[ctx->p][n];

  prt->xi[0]   = part->xi;
  prt->xi[1]   = part->yi;
  prt->xi[2]   = part->zi;
  prt->pxi[0]  = part->pxi;
  prt->pxi[1]  = part->pyi;
  prt->pxi[2]  = part->pzi;
  prt->kind    = part->kind_;
  prt->qni_wni = part->qni_wni;
}

static void
put_particle_single(struct cuda_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy_ctx *ctx = (struct copy_ctx *) _ctx;
  particle_single_t *part = &mparticles_single_t(ctx->mprts)[ctx->p][n];
  
  part->xi      = prt->xi[0];
  part->yi      = prt->xi[1];
  part->zi      = prt->xi[2];
  part->kind_   = prt->kind;
  part->pxi     = prt->pxi[0];
  part->pyi     = prt->pxi[1];
  part->pzi     = prt->pxi[2];
  part->qni_wni = prt->qni_wni;
}

static void
psc_mparticles_cuda_copy_from_single(struct psc_mparticles *mprts_cuda,
				    struct psc_mparticles *mprts, unsigned int flags)
{
  copy_from(mparticles_cuda_t(mprts_cuda), mprts, get_particle_single);
}

static void
psc_mparticles_cuda_copy_to_single(struct psc_mparticles *mprts_cuda,
				  struct psc_mparticles *mprts, unsigned int flags)
{
  copy_to(mparticles_cuda_t(mprts_cuda), mprts, put_particle_single);
}

// ======================================================================
// conversion to "double"

static void
get_particle_double(struct cuda_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy_ctx *ctx = (struct copy_ctx *) _ctx;
  particle_double_t *part = &mparticles_double_t(ctx->mprts)[ctx->p][n];

  prt->xi[0]   = part->xi;
  prt->xi[1]   = part->yi;
  prt->xi[2]   = part->zi;
  prt->pxi[0]  = part->pxi;
  prt->pxi[1]  = part->pyi;
  prt->pxi[2]  = part->pzi;
  prt->kind    = part->kind_;
  prt->qni_wni = part->qni_wni;
}

static void
put_particle_double(struct cuda_mparticles_prt *prt, int n, void *_ctx)
{
  struct copy_ctx *ctx = (struct copy_ctx *) _ctx;
  particle_double_t *part = &mparticles_double_t(ctx->mprts)[ctx->p][n];
  
  part->xi      = prt->xi[0];
  part->yi      = prt->xi[1];
  part->zi      = prt->xi[2];
  part->kind_   = prt->kind;
  part->pxi     = prt->pxi[0];
  part->pyi     = prt->pxi[1];
  part->pzi     = prt->pxi[2];
  part->qni_wni = prt->qni_wni;
}

static void
psc_mparticles_cuda_copy_from_double(struct psc_mparticles *mprts_cuda,
				    struct psc_mparticles *mprts, unsigned int flags)
{
  copy_from(mparticles_cuda_t(mprts_cuda), mprts, get_particle_double);
}

static void
psc_mparticles_cuda_copy_to_double(struct psc_mparticles *mprts_cuda,
				  struct psc_mparticles *mprts, unsigned int flags)
{
  copy_to(mparticles_cuda_t(mprts_cuda), mprts, put_particle_double);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_methods

static struct mrc_obj_method psc_mparticles_cuda_methods[] = {
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

  psc_mparticles_setup_super(mprts);

  int n_patches = mprts->nr_patches;
  assert(n_patches != 0);
  
  if (!mprts->flags) {
    // FIXME, they get set too late, so auto-dispatch "1vb" doesn't work
    mprts->flags = MP_BLOCKSIZE_4X4X4;
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

  mrc_json_t json = mrc_json_object_new(0);

  mrc_json_t info = mrc_json_object_new(0);
  mrc_json_object_push(json, "info", info);
  mrc_json_object_push_integer(info, "n_patches", n_patches);
  mrc_json_object_push_integer_array(info, "ldims", 3, ldims);
  mrc_json_object_push_integer_array(info, "bs", 3, bs);
  mrc_json_object_push_double_array(info, "dx", 3, dx);
  
  mrc_json_t json_xb_by_patch = mrc_json_array_new(n_patches);
  mrc_json_object_push(info, "xb_by_patch", json_xb_by_patch);
  for (int p = 0; p < n_patches; p++) {
    mrc_json_array_push_double_array(json_xb_by_patch, 3, ppsc->patch[p].xb);
  }

  double fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  mrc_json_object_push_double(info, "fnqs", fnqs);
  mrc_json_object_push_double(info, "eta", ppsc->coeff.eta);
  mrc_json_object_push_double(info, "dt", ppsc->dt);

  double kind_q[ppsc->nr_kinds], kind_m[ppsc->nr_kinds];
  for (int k = 0; k < ppsc->nr_kinds; k++) {
    kind_q[k] = ppsc->kinds[k].q;
    kind_m[k] = ppsc->kinds[k].m;
  }
  mrc_json_object_push_double_array(info, "kind_q", ppsc->nr_kinds, kind_q);
  mrc_json_object_push_double_array(info, "kind_m", ppsc->nr_kinds, kind_m);

  cuda_mparticles *cmprts = new cuda_mparticles(json);
  mprts_cuda->cmprts = cmprts;

  // FIXME json_builder_free(obj);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_destroy

static void
psc_mparticles_cuda_destroy(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  assert(mprts_cuda->cmprts);
  
  delete mprts_cuda->cmprts;
  mprts_cuda->cmprts = NULL;
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_reserve_all

static void
psc_mparticles_cuda_reserve_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  cmprts->reserve_all((unsigned int *) n_prts_by_patch);
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
  cmprts->get_size_all(n_prts_by_patch);

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
      float_4 *xi4  = (float_4 *) calloc(n_prts, sizeof(*xi4));
      float_4 *pxi4 = (float_4 *) calloc(n_prts, sizeof(*pxi4));
      
      cmprts->from_device(xi4, pxi4, n_prts, off);
      
      hsize_t hdims[2];
      hdims[0] = n_prts; hdims[1] = 4;
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
      float_4 *xi4  = (float_4*) calloc(n_prts, sizeof(float_4));
      float_4 *pxi4 = (float_4*) calloc(n_prts, sizeof(float_4));
      
      ierr = H5LTread_dataset_float(pgroup, "xi4", (float *) xi4); CE;
      ierr = H5LTread_dataset_float(pgroup, "pxi4", (float *) pxi4); CE;
      
      struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;
      
      cmprts->to_device(xi4, pxi4, n_prts, off);
      
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

  cmprts->setup_internals();
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_get_nr_particles

static unsigned int
psc_mparticles_cuda_get_nr_particles(struct psc_mparticles *mprts)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  return cmprts->get_n_prts();
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_get_size_all

static void
psc_mparticles_cuda_get_size_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  cmprts->get_size_all((unsigned int *) n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_resize_all

static void
psc_mparticles_cuda_resize_all(struct psc_mparticles *mprts, int *n_prts_by_patch)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  cmprts->resize_all((unsigned int *) n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_inject

void
psc_mparticles_cuda_inject(struct psc_mparticles *mprts, struct cuda_mparticles_prt *buf,
			   unsigned int *buf_n_by_patch)
{
  assert(strcmp(psc_mparticles_type(mprts), "cuda") == 0);
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  cmprts->inject(buf, buf_n_by_patch);
}

// ----------------------------------------------------------------------
// mparticles_cuda_t::patch_t::get_b_mx

const int* mparticles_cuda_t::patch_t::get_b_mx() const
{
  return mp_.sub()->cmprts->patch_get_b_mx(p_);
}

// ----------------------------------------------------------------------
// mparticles_cuda_t::patch_t::get_b_dxi

const mparticles_cuda_t::real_t* mparticles_cuda_t::patch_t::get_b_dxi() const
{
  return mp_.sub()->cmprts->patch_get_b_dxi(p_);
}

// ----------------------------------------------------------------------
// psc_mparticles: subclass "cuda"
  
struct psc_mparticles_ops_cuda : psc_mparticles_ops {
  psc_mparticles_ops_cuda() {
    name                    = "cuda";
    size                    = sizeof(struct psc_mparticles_cuda);
    methods                 = psc_mparticles_cuda_methods;
    setup                   = psc_mparticles_cuda_setup;
    destroy                 = psc_mparticles_cuda_destroy;
    read                    = psc_mparticles_cuda_read;
    write                   = psc_mparticles_cuda_write;
    setup_internals         = psc_mparticles_cuda_setup_internals;
    reserve_all             = psc_mparticles_cuda_reserve_all;
    get_nr_particles        = psc_mparticles_cuda_get_nr_particles;
    resize_all              = psc_mparticles_cuda_resize_all;
    get_size_all            = psc_mparticles_cuda_get_size_all;
  }
} psc_mparticles_cuda_ops;

