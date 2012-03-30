
#include "psc.h"
#include "psc_fields_c.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

void
fields_c_alloc(fields_c_t *pf, int ib[3], int ie[3], int nr_comp, int first_comp)
{
  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    pf->ib[d] = ib[d];
    pf->im[d] = ie[d] - ib[d];
    size *= pf->im[d];
  }
  pf->nr_comp = nr_comp;
  pf->first_comp = first_comp;
#ifdef USE_CBE
  // The Cell processor translation can use the C fields with one modification:
  // the data needs to be 128 byte aligned (to speed off-loading to spes). This
  // change is roughly put in below.
  void *m;
  int ierr = posix_memalign(&m, 128, nr_comp * size * sizeof(*pf->flds));
  pf->flds =  m; 
  assert(ierr == 0);
#else
  pf->flds = calloc(nr_comp * size, sizeof(*pf->flds));
#endif
}

void
fields_c_free(fields_c_t *pf)
{
  free(pf->flds);
}

void
fields_c_zero(fields_c_t *pf, int m)
{
  memset(&F3_C(pf, m, pf->ib[0], pf->ib[1], pf->ib[2]), 0,
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_c_real_t));
}

void
fields_c_set(fields_c_t *pf, int m, fields_c_real_t val)
{
  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3_C(pf, m, jx, jy, jz) = val;
      }
    }
  }
}

void
fields_c_copy_comp(fields_c_t *pto, int m_to, fields_c_t *pfrom, int m_from)
{
  for (int jz = pto->ib[2]; jz < pto->ib[2] + pto->im[2]; jz++) {
    for (int jy = pto->ib[1]; jy < pto->ib[1] + pto->im[1]; jy++) {
      for (int jx = pto->ib[0]; jx < pto->ib[0] + pto->im[0]; jx++) {
	F3_C(pto, m_to, jx, jy, jz) = F3_C(pfrom, m_from, jx, jy, jz);
      }
    }
  }
}

void
fields_c_axpy(fields_c_t *y, fields_c_real_t a, fields_c_t *x)
{
  assert(y->nr_comp == x->nr_comp);
  for (int m = 0; m < y->nr_comp; m++) {
    for (int jz = y->ib[2]; jz < y->ib[2] + y->im[2]; jz++) {
      for (int jy = y->ib[1]; jy < y->ib[1] + y->im[1]; jy++) {
	for (int jx = y->ib[0]; jx < y->ib[0] + y->im[0]; jx++) {
	  F3_C(y, m, jx, jy, jz) += a * F3_C(x, m, jx, jy, jz);
	}
      }
    }
  }
}

void
fields_c_axpy_comp(fields_c_t *y, int ym, fields_c_real_t a, fields_c_t *x, int xm)
{
  for (int jz = y->ib[2]; jz < y->ib[2] + y->im[2]; jz++) {
    for (int jy = y->ib[1]; jy < y->ib[1] + y->im[1]; jy++) {
      for (int jx = y->ib[0]; jx < y->ib[0] + y->im[0]; jx++) {
	F3_C(y, ym, jx, jy, jz) += a * F3_C(x, xm, jx, jy, jz);
      }
    }
  }
}

void
fields_c_scale(fields_c_t *pf, fields_c_real_t val)
{
  for (int m = 0; m < pf->nr_comp; m++) {
    for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
      for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
	for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	  F3_C(pf, m, jx, jy, jz) *= val;
	}
      }
    }
  }
}

static void
_psc_mfields_c_zero_comp(mfields_c_t *flds, int m)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    fields_c_zero(psc_mfields_get_patch_c(flds, p), m);
  }
}

static void
_psc_mfields_c_axpy(mfields_c_t *yf, fields_c_real_t alpha, mfields_c_t *xf)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    fields_c_axpy(psc_mfields_get_patch_c(yf, p), alpha, psc_mfields_get_patch_c(xf, p));
  }
}

static void
_psc_mfields_c_axpy_comp(mfields_c_t *yf, int ym, fields_c_real_t alpha,
			 mfields_c_t *xf, int xm)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    fields_c_axpy_comp(psc_mfields_get_patch_c(yf, p), ym, alpha,
		       psc_mfields_get_patch_c(xf, p), xm);
  }
}

static void
_psc_mfields_c_scale(mfields_c_t *yf, double alpha)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    fields_c_scale(psc_mfields_get_patch_c(yf, p), alpha);
  }
}

static void
_psc_mfields_c_set_comp(mfields_c_t *yf, int m, fields_c_real_t alpha)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    fields_c_set(psc_mfields_get_patch_c(yf, p), m, alpha);
  }
}

static void
_psc_mfields_c_copy_comp(mfields_c_t *to, int mto, mfields_c_t *from, int mfrom)
{
  for (int p = 0; p < to->nr_patches; p++) {
    fields_c_copy_comp(psc_mfields_get_patch_c(to, p), mto,
		       psc_mfields_get_patch_c(from, p), mfrom);
  }
}

// ======================================================================
// psc_mfields_c

static void
_psc_mfields_c_setup(mfields_c_t *flds)
{
  psc_mfields_setup_super(flds);

  struct mrc_patch *patches = mrc_domain_get_patches(flds->domain,
						     &flds->nr_patches);
  flds->data = calloc(flds->nr_patches, sizeof(fields_c_t));
  for (int p = 0; p < flds->nr_patches; p++) {
    int ilg[3] = { -flds->ibn[0], -flds->ibn[1], -flds->ibn[2] };
    int ihg[3] = { patches[p].ldims[0] + flds->ibn[0],
		   patches[p].ldims[1] + flds->ibn[1],
		   patches[p].ldims[2] + flds->ibn[2] };
    fields_c_alloc(psc_mfields_get_patch_c(flds, p), ilg, ihg,
		   flds->nr_fields, flds->first_comp);
  }
}

static void
_psc_mfields_c_destroy(mfields_c_t *flds)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    fields_c_free(psc_mfields_get_patch_c(flds, p));
  }
  free(flds->data);
}

#ifdef HAVE_LIBHDF5_HL

#include <mrc_io.h>

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

static void
_psc_mfields_c_write(mfields_c_t *mfields, struct mrc_io *io)
{
  int ierr;
  const char *path = psc_mfields_name(mfields);
  mrc_io_write_obj_ref(io, path, "domain", (struct mrc_obj *) mfields->domain);
  mrc_io_write_attr_int(io, path, "nr_patches", mfields->nr_patches);

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, path, H5P_DEFAULT); H5_CHK(group);
  for (int m = 0; m < mfields->nr_fields; m++) {
    char namec[10]; sprintf(namec, "m%d", m);
    const char *s = psc_mfields_comp_name(mfields, m);
    if (!s) {
      s = "(null)";
    }
    ierr = H5LTset_attribute_string(group, ".", namec, s); CE;
  }
  for (int p = 0; p < mfields->nr_patches; p++) {
    fields_c_t *fields = psc_mfields_get_patch_c(mfields, p);
    char name[10]; sprintf(name, "p%d", p);

    hid_t groupp = H5Gcreate(group, name, H5P_DEFAULT, H5P_DEFAULT,
			     H5P_DEFAULT); H5_CHK(groupp);
    ierr = H5LTset_attribute_int(groupp, ".", "ib", fields->ib, 3); CE;
    ierr = H5LTset_attribute_int(groupp, ".", "im", fields->im, 3); CE;
    ierr = H5LTset_attribute_int(groupp, ".", "nr_comp", &fields->nr_comp, 1); CE;
    // write components separately instead?
    hsize_t hdims[4] = { fields->nr_comp, fields->im[2], fields->im[1], fields->im[0] };
    ierr = H5LTmake_dataset_double(groupp, "fields_c", 4, hdims, fields->flds); CE;
    ierr = H5Gclose(groupp); CE;
  }

  ierr = H5Gclose(group); CE;
}

static void
_psc_mfields_c_read(mfields_c_t *mfields, struct mrc_io *io)
{
  int ierr;
  const char *path = psc_mfields_name(mfields);
  mfields->domain = (struct mrc_domain *)
    mrc_io_read_obj_ref(io, path, "domain", &mrc_class_mrc_domain);
  mrc_io_read_attr_int(io, path, "nr_patches", &mfields->nr_patches);
  psc_mfields_setup(mfields);

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, path, H5P_DEFAULT); H5_CHK(group);
  for (int m = 0; m < mfields->nr_fields; m++) {
    char namec[10]; sprintf(namec, "m%d", m);
    
    hsize_t dims;
    H5T_class_t class;
    size_t sz;
    ierr = H5LTget_attribute_info(group, ".", namec, &dims, &class, &sz); CE;
    char *s = malloc(sz);
    ierr = H5LTget_attribute_string(group, ".", namec, s); CE;
    if (strcmp(s, "(null)") != 0) {
      psc_mfields_set_comp_name(mfields, m, s);
    }
  }

  for (int p = 0; p < mfields->nr_patches; p++) {
    fields_c_t *fields = psc_mfields_get_patch_c(mfields, p);
    char name[10]; sprintf(name, "p%d", p);

    hid_t groupp = H5Gopen(group, name, H5P_DEFAULT); H5_CHK(groupp);
    int ib[3], im[3], nr_comp;
    ierr = H5LTget_attribute_int(groupp, ".", "ib", ib); CE;
    ierr = H5LTget_attribute_int(groupp, ".", "im", im); CE;
    ierr = H5LTget_attribute_int(groupp, ".", "nr_comp", &nr_comp); CE;
    for (int d = 0; d < 3; d++) {
      assert(ib[d] == fields->ib[d]);
      assert(im[d] == fields->im[d]);
    }
    assert(nr_comp == fields->nr_comp);
    ierr = H5LTread_dataset_double(groupp, "fields_c", fields->flds); CE;
    ierr = H5Gclose(groupp); CE;
  }

  ierr = H5Gclose(group); CE;
}

#endif

static void
psc_mfields_c_copy_to_fortran(mfields_c_t *flds_c, mfields_fortran_t *flds_fortran, int mb, int me)
{
  psc_mfields_fortran_copy_from_c(flds_fortran, flds_c, mb, me);
}

static void
psc_mfields_c_copy_from_fortran(mfields_c_t *flds_c, mfields_fortran_t *flds_fortran, int mb, int me)
{
  psc_mfields_fortran_copy_to_c(flds_fortran, flds_c, mb, me);
}

#ifdef USE_CUDA
static void
psc_mfields_c_copy_to_cuda(mfields_c_t *flds_c, mfields_cuda_t *flds_cuda, int mb, int me)
{
  psc_mfields_cuda_copy_from_c(flds_cuda, flds_c, mb, me);
}

static void
psc_mfields_c_copy_from_cuda(mfields_c_t *flds_c, mfields_cuda_t *flds_cuda, int mb, int me)
{
  psc_mfields_cuda_copy_to_c(flds_cuda, flds_c, mb, me);
}
#endif

// ======================================================================
// psc_mfields: subclass "c"
  
struct psc_mfields_ops psc_mfields_c_ops = {
  .name                  = "c",
  .setup                 = _psc_mfields_c_setup,
  .destroy               = _psc_mfields_c_destroy,
#ifdef HAVE_LIBHDF5_HL
  .write                 = _psc_mfields_c_write,
  .read                  = _psc_mfields_c_read,
#endif
  .zero_comp             = _psc_mfields_c_zero_comp,
  .set_comp              = _psc_mfields_c_set_comp,
  .scale                 = _psc_mfields_c_scale,
  .copy_comp             = _psc_mfields_c_copy_comp,
  .axpy                  = _psc_mfields_c_axpy,
  .axpy_comp             = _psc_mfields_c_axpy_comp,
  .copy_to_fortran       = psc_mfields_c_copy_to_fortran,
  .copy_from_fortran     = psc_mfields_c_copy_from_fortran,
#ifdef USE_CUDA
  .copy_to_cuda          = psc_mfields_c_copy_to_cuda,
  .copy_from_cuda        = psc_mfields_c_copy_from_cuda,
#endif
};

