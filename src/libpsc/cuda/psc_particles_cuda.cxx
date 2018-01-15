
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
  unsigned int n_prts_by_patch[mprts->n_patches];
  mprts->get_size_all(n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    int n_prts = n_prts_by_patch[p];
    struct copy_ctx ctx = { .mprts = mprts_from, .p = p };
    mprts->set_particles(n_prts, off, get_particle, &ctx);

    off += n_prts;
  }
}

static void
copy_to(mparticles_cuda_t mprts, struct psc_mparticles *mprts_to,
	void (*put_particle)(struct cuda_mparticles_prt *prt, int n, void *ctx))
{
  unsigned int n_prts_by_patch[mprts->n_patches];
  mprts->get_size_all(n_prts_by_patch);
  
  unsigned int off = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    int n_prts = n_prts_by_patch[p];
    struct copy_ctx ctx = { .mprts = mprts_to, .p = p };
    mprts->get_particles(n_prts, off, put_particle, &ctx);

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
psc_mparticles_cuda_setup(struct psc_mparticles *_mprts)
{
  mparticles_cuda_t mprts(_mprts);

  cuda_base_init();

  psc_mparticles_setup_super(_mprts);

  int n_patches = mprts.n_patches();
  assert(n_patches != 0);
  
  if (!_mprts->flags) {
    // FIXME, they get set too late, so auto-dispatch "1vb" doesn't work
    _mprts->flags = MP_BLOCKSIZE_4X4X4;
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
    switch (_mprts->flags & MP_BLOCKSIZE_MASK) {
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

  Grid_t& grid = ppsc->grid;
  grid.gdims = ppsc->domain.gdims;
  grid.ldims = ldims;
  grid.dx = dx;
  grid.fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  grid.eta = ppsc->coeff.eta;
  grid.dt = ppsc->dt;
  
  grid.patches.resize(n_patches);
  for (int p = 0; p < n_patches; p++) {
    grid.patches[p].xb = ppsc->patch[p].xb;
  }

  for (int k = 0; k < ppsc->nr_kinds; k++) {
    grid.kinds.push_back(Grid_t::Kind(ppsc->kinds[k].q, ppsc->kinds[k].m, ppsc->kinds[k].name));
  }

  new(mprts.sub_) cuda_mparticles(grid, bs);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_destroy

static void
psc_mparticles_cuda_destroy(struct psc_mparticles *_mprts)
{
  mparticles_cuda_t mprts(_mprts);
  
  mprts.sub_->~cuda_mparticles();
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_reserve_all

static void
psc_mparticles_cuda_reserve_all(struct psc_mparticles *_mprts, int *n_prts_by_patch)
{
  mparticles_cuda_t mprts(_mprts);

  mprts->reserve_all((unsigned int *) n_prts_by_patch);
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
psc_mparticles_cuda_write(struct psc_mparticles *_mprts, struct mrc_io *io)
{
  mparticles_cuda_t mprts(_mprts);
  int ierr;
  
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  unsigned int n_prts_by_patch[mprts.n_patches()];
  mprts->get_size_all(n_prts_by_patch);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, _mprts), H5P_DEFAULT); H5_CHK(group);
  unsigned int off = 0;
  // FIXME, reorder first if necessary
  for (int p = 0; p < mprts.n_patches(); p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gcreate(group, pname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = n_prts_by_patch[p];
    ierr = H5LTset_attribute_int(pgroup, ".", "n_prts", &n_prts, 1); CE;
    if (n_prts > 0) {
      float_4 *xi4  = (float_4 *) calloc(n_prts, sizeof(*xi4));
      float_4 *pxi4 = (float_4 *) calloc(n_prts, sizeof(*pxi4));
      
      mprts->from_device(xi4, pxi4, n_prts, off);
      
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
psc_mparticles_cuda_read(struct psc_mparticles *_mprts, struct mrc_io *io)
{
  mparticles_cuda_t mprts(_mprts);

  psc_mparticles_read_super(_mprts, io);

  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, _mprts), H5P_DEFAULT); H5_CHK(group);

  int n_prts_by_patch[mprts.n_patches()];

  for (int p = 0; p < mprts.n_patches(); p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts;
    ierr = H5LTget_attribute_int(pgroup, ".", "n_prts", &n_prts); CE;
    n_prts_by_patch[p] = n_prts;
    ierr = H5Gclose(pgroup); CE;
  }

  psc_mparticles_setup(_mprts);
  psc_mparticles_reserve_all(_mprts, n_prts_by_patch);

  unsigned int off = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = n_prts_by_patch[p];
    if (n_prts > 0) {
      float_4 *xi4  = (float_4*) calloc(n_prts, sizeof(float_4));
      float_4 *pxi4 = (float_4*) calloc(n_prts, sizeof(float_4));
      
      ierr = H5LTread_dataset_float(pgroup, "xi4", (float *) xi4); CE;
      ierr = H5LTread_dataset_float(pgroup, "pxi4", (float *) pxi4); CE;
      
      mparticles_cuda_t mprts = mparticles_cuda_t(_mprts);
      mprts->to_device(xi4, pxi4, n_prts, off);
      
      free(xi4);
      free(pxi4);
    }

    ierr = H5Gclose(pgroup); CE;
    off += n_prts;
  }

  ierr = H5Gclose(group); CE;
  psc_mparticles_setup_internals(_mprts);
}

#endif

// ----------------------------------------------------------------------
// psc_mparticles_cuda_setup_internals

static void
psc_mparticles_cuda_setup_internals(struct psc_mparticles *_mprts)
{
  mparticles_cuda_t mprts(_mprts);

  mprts->setup_internals();
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_get_nr_particles

static unsigned int
psc_mparticles_cuda_get_nr_particles(struct psc_mparticles *_mprts)
{
  mparticles_cuda_t mprts(_mprts);

  return mprts->get_n_prts();
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_get_size_all

static void
psc_mparticles_cuda_get_size_all(struct psc_mparticles *_mprts, int *n_prts_by_patch)
{
  mparticles_cuda_t mprts(_mprts);

  mprts->get_size_all((unsigned int *) n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_resize_all

static void
psc_mparticles_cuda_resize_all(struct psc_mparticles *_mprts, int *n_prts_by_patch)
{
  mparticles_cuda_t mprts(_mprts);

  mprts->resize_all((unsigned int *) n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_inject

void
psc_mparticles_cuda_inject(struct psc_mparticles *_mprts, struct cuda_mparticles_prt *buf,
			   unsigned int *buf_n_by_patch)
{
  mparticles_cuda_t mprts(_mprts);

  mprts->inject(buf, buf_n_by_patch);
}

// ----------------------------------------------------------------------
// mparticles_cuda_t::patch_t::get_b_mx

const int* mparticles_cuda_t::patch_t::get_b_mx() const
{
  return mp_->patch_get_b_mx(p_);
}

// ----------------------------------------------------------------------
// mparticles_cuda_t::patch_t::get_b_dxi

const mparticles_cuda_t::real_t* mparticles_cuda_t::patch_t::get_b_dxi() const
{
  return mp_->patch_get_b_dxi(p_);
}

// ----------------------------------------------------------------------
// psc_mparticles: subclass "cuda"
  
struct psc_mparticles_ops_cuda : psc_mparticles_ops {
  psc_mparticles_ops_cuda() {
    name                    = "cuda";
    size                    = sizeof(struct cuda_mparticles);
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

