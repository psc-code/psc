
#include "psc.h"
#include "cuda_mparticles.h"
#include "psc_particles_cuda.h"
#include "psc_push_particles.h"
#include "psc_particles_single.h"
#include "psc_particles_double.h"

#include <mrc_io.h>

#if 0
#define dprintf(...) mprintf(__VA_ARGS__)
#else
#define dprintf(...) do {} while (0)
#endif

// ======================================================================
// MparticleCuda implementation

MparticlesCuda::MparticlesCuda(const Grid_t& grid)
  : MparticlesBase(grid),
    pi_(grid)
{
  dprintf("CMPRTS: ctor\n");
  cmprts_ = new cuda_mparticles<BS144>(grid);
}

MparticlesCuda::~MparticlesCuda()
{
  dprintf("CMPRTS: dtor\n");
  delete cmprts_;
}

void MparticlesCuda::reserve_all(const uint *n_prts_by_patch)
{
  dprintf("CMPRTS: reserve_all\n");
  for (int p = 0; p < cmprts_->n_patches; p++) {
    dprintf("  p %d: %d\n", p, n_prts_by_patch[p]);
  }
  cmprts_->reserve_all(n_prts_by_patch);
}

void MparticlesCuda::get_size_all(uint *n_prts_by_patch) const
{
  dprintf("CMPRTS: get_size_all\n");
  cmprts_->get_size_all(n_prts_by_patch);
  for (int p = 0; p < cmprts_->n_patches; p++) {
    dprintf("  p %d: %d\n", p, n_prts_by_patch[p]);
  }
}

void MparticlesCuda::resize_all(const uint *n_prts_by_patch)
{
  dprintf("CMPRTS: resize_all\n");
  cmprts_->resize_all(n_prts_by_patch);
}

int MparticlesCuda::get_n_prts() const
{
  dprintf("CMPRTS: get_n_prts\n");
  return cmprts_->get_n_prts();
}

void MparticlesCuda::reset(const Grid_t& grid)
{
  this->~MparticlesCuda();
  new(this) MparticlesCuda(grid);
}

void MparticlesCuda::inject_buf(cuda_mparticles_prt *buf, uint *buf_n_by_patch)
{
  dprintf("CMPRTS: inject\n");
  cmprts_->inject(buf, buf_n_by_patch);
}

void MparticlesCuda::dump(const std::string& filename)
{
  cmprts_->dump(filename);
}

// ======================================================================
// conversion

template<typename MP>
struct ConvertToCuda
{
  using particle_t = typename MP::particle_t;

  ConvertToCuda(MP& mprts_other, int p)
    : mprts_other_(mprts_other), p_(p)
  {}

  cuda_mparticles_prt operator()(int n)
  {
    const particle_t& prt_other = mprts_other_[p_][n];

    cuda_mparticles_prt prt;
    prt.xi[0]   = prt_other.xi;
    prt.xi[1]   = prt_other.yi;
    prt.xi[2]   = prt_other.zi;
    prt.pxi[0]  = prt_other.pxi;
    prt.pxi[1]  = prt_other.pyi;
    prt.pxi[2]  = prt_other.pzi;
    prt.kind    = prt_other.kind_;
    prt.qni_wni = prt_other.qni_wni_;

    return prt;
  }

private:
  MP& mprts_other_;
  int p_;
};

template<typename MP>
struct ConvertFromCuda
{
  using particle_t = typename MP::particle_t;

  ConvertFromCuda(MP& mprts_other, int p)
    : mprts_other_(mprts_other), p_(p)
  {}

  void operator()(int n, const cuda_mparticles_prt &prt)
  {
    particle_t& prt_other = mprts_other_[p_][n];

    prt_other.xi      = prt.xi[0];
    prt_other.yi      = prt.xi[1];
    prt_other.zi      = prt.xi[2];
    prt_other.kind_   = prt.kind;
    prt_other.pxi     = prt.pxi[0];
    prt_other.pyi     = prt.pxi[1];
    prt_other.pzi     = prt.pxi[2];
    prt_other.qni_wni_ = prt.qni_wni;
  }

private:
  MP& mprts_other_;
  int p_;
};

template<typename MP>
static void copy_from(MparticlesCuda& mp, MP& mp_other)
{
  int n_patches = mp.n_patches();
  uint n_prts_by_patch[n_patches];
  mp_other.get_size_all(n_prts_by_patch);
  mp.reserve_all(n_prts_by_patch);
  mp.resize_all(n_prts_by_patch);

  for (int p = 0; p < n_patches; p++) {
    ConvertToCuda<MP> convert_to_cuda(mp_other, p);
    mp.cmprts()->set_particles(p, convert_to_cuda);
  }

  mp.cmprts()->setup_internals();
}

template<typename MP>
static void copy_to(MparticlesCuda& mp, MP& mp_other)
{
  int n_patches = mp_other.n_patches();
  uint n_prts_by_patch[n_patches];
  mp.get_size_all(n_prts_by_patch);
  mp_other.reserve_all(n_prts_by_patch);
  mp_other.resize_all(n_prts_by_patch);

  if (mp.cmprts()->need_reorder) {
    mp.cmprts()->reorder();
  }
  for (int p = 0; p < n_patches; p++) {
    ConvertFromCuda<MP> convert_from_cuda(mp_other, p);
    mp.cmprts()->get_particles(p, convert_from_cuda);
  }
}

// ======================================================================
// conversion to "single"/"double"

template<typename MP>
static void copy_from(MparticlesBase& mp, MparticlesBase& mp_other)
{
  copy_from(dynamic_cast<MparticlesCuda&>(mp), dynamic_cast<MP&>(mp_other));
}

template<typename MP>
static void copy_to(MparticlesBase& mp, MparticlesBase& mp_other)
{
  copy_to(dynamic_cast<MparticlesCuda&>(mp), dynamic_cast<MP&>(mp_other));
}

const MparticlesCuda::Convert MparticlesCuda::convert_to_ = {
  { std::type_index(typeid(MparticlesSingle)), copy_to<MparticlesSingle>   },
  { std::type_index(typeid(MparticlesDouble)), copy_to<MparticlesDouble>   },
};

const MparticlesCuda::Convert MparticlesCuda::convert_from_ = {
  { std::type_index(typeid(MparticlesSingle)), copy_from<MparticlesSingle>   },
  { std::type_index(typeid(MparticlesDouble)), copy_from<MparticlesDouble>   },
};

// ======================================================================
// psc_mparticles "cuda"

#if 0
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
  PscMparticlesCuda mprts(_mprts);
  int ierr;
  
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  uint n_prts_by_patch[mprts.n_patches()];
  mprts->get_size_all(n_prts_by_patch);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, _mprts), H5P_DEFAULT); H5_CHK(group);
  uint off = 0;
  // FIXME, reorder first if necessary
  for (int p = 0; p < mprts.n_patches(); p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gcreate(group, pname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = n_prts_by_patch[p];
    ierr = H5LTset_attribute_int(pgroup, ".", "n_prts", &n_prts, 1); CE;
    if (n_prts > 0) {
      float4 *xi4  = (float4 *) calloc(n_prts, sizeof(*xi4));
      float4 *pxi4 = (float4 *) calloc(n_prts, sizeof(*pxi4));
      
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
  PscMparticlesCuda mprts(_mprts);

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
  mprts->reserve_all(n_prts_by_patch);

  uint off = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = n_prts_by_patch[p];
    if (n_prts > 0) {
      float4 *xi4  = (float4*) calloc(n_prts, sizeof(float4));
      float4 *pxi4 = (float4*) calloc(n_prts, sizeof(float4));
      
      ierr = H5LTread_dataset_float(pgroup, "xi4", (float *) xi4); CE;
      ierr = H5LTread_dataset_float(pgroup, "pxi4", (float *) pxi4); CE;
      
      PscMparticlesCuda mprts = PscMparticlesCuda(_mprts);
      mprts->to_device(xi4, pxi4, n_prts, off);
      
      free(xi4);
      free(pxi4);
    }

    ierr = H5Gclose(pgroup); CE;
    off += n_prts;
  }

  ierr = H5Gclose(group); CE;
  mprts->setup_internals();
}

#endif
#endif

// ----------------------------------------------------------------------
// psc_mparticles: subclass "cuda"
  
psc_mparticles_ops_<MparticlesCuda> psc_mparticles_cuda_ops;
