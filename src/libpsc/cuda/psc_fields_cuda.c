
#include "psc.h"
#include "cuda_iface.h"
#include "psc_fields_cuda.h"
#include "psc_fields_c.h"
#include "psc_fields_single.h"
#include "psc_cuda.h"

#include <mrc_params.h>

// OPT, CUDA fields have too many ghostpoints, and 7 points in the invar direction!

// ----------------------------------------------------------------------
// macros to access C (host) versions of the fields

#define F3_OFF_CUDA(pf, fldnr, jx,jy,jz)				\
  ((((((fldnr)								\
       * (pf)->im[2] + ((jz)-(pf)->ib[2]))				\
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))				\
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))

#ifndef BOUNDS_CHECK

#define F3_CUDA(pf, fldnr, jx,jy,jz)		\
  (h_flds[F3_OFF_CUDA(pf, fldnr, jx,jy,jz)])

#else

#define F3_CUDA(pf, fldnr, jx,jy,jz)				\
  (*({int off = F3_OFF_CUDA(pf, fldnr, jx,jy,jz);			\
      assert(fldnr >= 0 && fldnr < (pf)->nr_comp);			\
      assert(jx >= (pf)->ib[0] && jx < (pf)->ib[0] + (pf)->im[0]);	\
      assert(jy >= (pf)->ib[1] && jy < (pf)->ib[1] + (pf)->im[1]);	\
      assert(jz >= (pf)->ib[2] && jz < (pf)->ib[2] + (pf)->im[2]);	\
      &(h_flds[off]);						\
    }))

#endif

// ======================================================================

// ----------------------------------------------------------------------
// psc_fields_cuda_axpy_comp

static void
psc_fields_cuda_axpy_comp(struct psc_fields *y, int ym, double a, struct psc_fields *x, int xm)
{
  assert(ppsc->domain.gdims[0] == 1);
  cuda_axpy_comp_yz(y, ym, a, x, xm);
}

// ----------------------------------------------------------------------
// psc_fields_cuda_zero_comp

static void
psc_fields_cuda_zero_comp(struct psc_fields *x, int xm)
{
  assert(ppsc->domain.gdims[0] == 1);
  cuda_zero_comp_yz(x, xm);
}

#ifdef HAVE_LIBHDF5_HL

#include <mrc_io.h>

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------
// psc_fields_cuda_write

static void
psc_fields_cuda_write(struct psc_fields *flds, struct mrc_io *io)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, flds), H5P_DEFAULT); H5_CHK(group);
  ierr = H5LTset_attribute_int(group, ".", "p", &flds->p, 1); CE;
  ierr = H5LTset_attribute_int(group, ".", "ib", flds->ib, 3); CE;
  ierr = H5LTset_attribute_int(group, ".", "im", flds->im, 3); CE;
  ierr = H5LTset_attribute_int(group, ".", "nr_comp", &flds->nr_comp, 1); CE;
  // write components separately instead?
  hsize_t hdims[4] = { flds->nr_comp, flds->im[2], flds->im[1], flds->im[0] };
  float *h_flds = malloc(flds->nr_comp * psc_fields_size(flds) * sizeof(*h_flds));
  __fields_cuda_from_device(flds, h_flds, 0, flds->nr_comp);
  ierr = H5LTmake_dataset_float(group, "fields_cuda", 4, hdims, h_flds); CE;
  free(h_flds);
  ierr = H5Gclose(group); CE;
}

#endif

// ======================================================================
// convert from/to "c"

static void
psc_fields_cuda_copy_from_c(struct psc_fields *flds_cuda, struct psc_fields *flds_c,
			    int mb, int me)
{
  float *h_flds = malloc(flds_cuda->nr_comp * psc_fields_size(flds_cuda) * sizeof(*h_flds));

  for (int m = mb; m < me; m++) {
    for (int jz = flds_cuda->ib[2]; jz < flds_cuda->ib[2] + flds_cuda->im[2]; jz++) {
      for (int jy = flds_cuda->ib[1]; jy < flds_cuda->ib[1] + flds_cuda->im[1]; jy++) {
	for (int jx = flds_cuda->ib[0]; jx < flds_cuda->ib[0] + flds_cuda->im[0]; jx++) {
	  F3_CUDA(flds_cuda, m, jx,jy,jz) = F3_C(flds_c, m, jx,jy,jz);
	}
      }
    }
  }

  __fields_cuda_to_device(flds_cuda, h_flds, mb, me);

  free(h_flds);
}

static void
psc_fields_cuda_copy_to_c(struct psc_fields *flds_cuda, struct psc_fields *flds_c,
			  int mb, int me)
{
  float *h_flds = malloc(flds_cuda->nr_comp * psc_fields_size(flds_cuda) * sizeof(*h_flds));

  __fields_cuda_from_device(flds_cuda, h_flds, mb, me);
  
  for (int m = mb; m < me; m++) {
    for (int jz = flds_cuda->ib[2]; jz < flds_cuda->ib[2] + flds_cuda->im[2]; jz++) {
      for (int jy = flds_cuda->ib[1]; jy < flds_cuda->ib[1] + flds_cuda->im[1]; jy++) {
	for (int jx = flds_cuda->ib[0]; jx < flds_cuda->ib[0] + flds_cuda->im[0]; jx++) {
	  F3_C(flds_c, m, jx,jy,jz) = F3_CUDA(flds_cuda, m, jx,jy,jz);
	}
      }
    }
  }

  free(h_flds);
}

// ======================================================================
// convert from/to "single"

static void
psc_fields_cuda_copy_from_single(struct psc_fields *flds_cuda, struct psc_fields *flds_single,
				 int mb, int me)
{
  float *h_flds = malloc(flds_cuda->nr_comp * psc_fields_size(flds_cuda) * sizeof(*h_flds));

  for (int m = mb; m < me; m++) {
    for (int jz = flds_cuda->ib[2]; jz < flds_cuda->ib[2] + flds_cuda->im[2]; jz++) {
      for (int jy = flds_cuda->ib[1]; jy < flds_cuda->ib[1] + flds_cuda->im[1]; jy++) {
	for (int jx = flds_cuda->ib[0]; jx < flds_cuda->ib[0] + flds_cuda->im[0]; jx++) {
	  F3_CUDA(flds_cuda, m, jx,jy,jz) = F3_S(flds_single, m, jx,jy,jz);
	}
      }
    }
  }

  __fields_cuda_to_device(flds_cuda, h_flds, mb, me);

  free(h_flds);
}

static void
psc_fields_cuda_copy_to_single(struct psc_fields *flds_cuda, struct psc_fields *flds_single,
			       int mb, int me)
{
  float *h_flds = malloc(flds_cuda->nr_comp * psc_fields_size(flds_cuda) * sizeof(*h_flds));

  __fields_cuda_from_device(flds_cuda, h_flds, mb, me);
  
  for (int m = mb; m < me; m++) {
    for (int jz = flds_cuda->ib[2]; jz < flds_cuda->ib[2] + flds_cuda->im[2]; jz++) {
      for (int jy = flds_cuda->ib[1]; jy < flds_cuda->ib[1] + flds_cuda->im[1]; jy++) {
	for (int jx = flds_cuda->ib[0]; jx < flds_cuda->ib[0] + flds_cuda->im[0]; jx++) {
	  F3_S(flds_single, m, jx,jy,jz) = F3_CUDA(flds_cuda, m, jx,jy,jz);
	}
      }
    }
  }

  free(h_flds);
}

// ======================================================================
// psc_fields: subclass "cuda"
  
static struct mrc_obj_method psc_fields_cuda_methods[] = {
  MRC_OBJ_METHOD("copy_to_c"       , psc_fields_cuda_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c"     , psc_fields_cuda_copy_from_c),
  MRC_OBJ_METHOD("copy_to_single"  , psc_fields_cuda_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_fields_cuda_copy_from_single),
  {}
};

struct psc_fields_ops psc_fields_cuda_ops = {
  .name                  = "cuda",
  .size                  = sizeof(struct psc_fields_cuda),
  .methods               = psc_fields_cuda_methods,
#ifdef HAVE_LIBHDF5_HL
  .write                 = psc_fields_cuda_write,
#endif
  .axpy_comp             = psc_fields_cuda_axpy_comp,
  .zero_comp             = psc_fields_cuda_zero_comp,
};

// ======================================================================

// ----------------------------------------------------------------------
// psc_mfields_cuda_setup

static void
psc_mfields_cuda_setup(struct psc_mfields *mflds)
{
  psc_mfields_setup_super(mflds);

  cuda_base_init();

  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  struct cuda_mfields *cmflds = cuda_mfields_create();
  mflds_cuda->cmflds = cmflds;

  __psc_mfields_cuda_setup(mflds);
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_destroy

static void
psc_mfields_cuda_destroy(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  struct cuda_mfields *cmflds = mflds_cuda->cmflds;

  __psc_mfields_cuda_destroy(mflds);

  cuda_mfields_destroy(cmflds);
  mflds_cuda->cmflds = NULL;
}

#ifdef HAVE_LIBHDF5_HL

// FIXME, it'd be nicer to have this a psc_fields::cuda::read() method,
// but it would need to know the parent to ask for the allocated memory,
// or something...

// ----------------------------------------------------------------------
// psc_fields_cuda_read

static void
psc_fields_cuda_read(struct psc_fields *flds, struct mrc_io *io, const char *path)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, path, H5P_DEFAULT); H5_CHK(group);
  int ib[3], im[3], nr_comp;
  ierr = H5LTget_attribute_int(group, ".", "p", &flds->p); CE;
  ierr = H5LTget_attribute_int(group, ".", "ib", ib); CE;
  ierr = H5LTget_attribute_int(group, ".", "im", im); CE;
  ierr = H5LTget_attribute_int(group, ".", "nr_comp", &nr_comp); CE;
  for (int d = 0; d < 3; d++) {
    assert(ib[d] == flds->ib[d]);
    assert(im[d] == flds->im[d]);
  }
  assert(nr_comp == flds->nr_comp);
  //  psc_fields_setup(flds);
  float *h_flds = malloc(flds->nr_comp * psc_fields_size(flds) * sizeof(*h_flds));
  ierr = H5LTread_dataset_float(group, "fields_cuda", h_flds); CE;
  __fields_cuda_to_device(flds, h_flds, 0, flds->nr_comp);
  free(h_flds);
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_read

static void
psc_mfields_cuda_read(struct psc_mfields *mflds, struct mrc_io *io)
{
  mflds->domain = mrc_io_read_ref(io, mflds, "domain", mrc_domain);
  mrc_domain_get_patches(mflds->domain, &mflds->nr_patches);

  mflds->comp_name = calloc(mflds->nr_fields, sizeof(*mflds->comp_name));
  for (int m = 0; m < mflds->nr_fields; m++) {
    char name[20]; sprintf(name, "comp_name_%d", m);
    char *s;
    mrc_io_read_string(io, mflds, name, &s);
    if (s) {
      psc_mfields_set_comp_name(mflds, m, s);
    }
  }

  psc_mfields_cuda_setup(mflds);

  for (int p = 0; p < mflds->nr_patches; p++) {
    char name[20]; sprintf(name, "flds%d", p);
    char *path;
    mrc_io_read_attr_string(io, mrc_io_obj_path(io, mflds), name, &path);
    psc_fields_cuda_read(mflds->flds[p], io, path);
  }
}

#endif

// ======================================================================
// psc_mfields: subclass "cuda"
  
struct psc_mfields_ops psc_mfields_cuda_ops = {
  .name                  = "cuda",
  .size                  = sizeof(struct psc_mfields_cuda),
  .setup                 = psc_mfields_cuda_setup,
  .destroy               = psc_mfields_cuda_destroy,
#ifdef HAVE_LIBHDF5_HL
  .read                  = psc_mfields_cuda_read,
#endif
};

