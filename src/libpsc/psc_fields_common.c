
// ----------------------------------------------------------------------
// psc_fields_setup

static void
PFX(setup)(struct psc_fields *pf)
{
  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    size *= pf->im[d];
  }

#if PSC_FIELDS_AS_FORTRAN
  fields_real_t **flds = calloc(pf->nr_comp, sizeof(*flds));
  flds[0] = calloc(size * pf->nr_comp, sizeof(flds[0]));
  for (int i = 1; i < pf->nr_comp; i++) {
    flds[i] = flds[0] + i * size;
  }
  pf->data = flds;
#elif PSC_FIELDS_AS_C && defined(USE_CBE)
  // The Cell processor translation can use the C fields with one modification:
  // the data needs to be 128 byte aligned (to speed off-loading to spes). This
  // change is roughly put in below.
  void *m;
  int ierr = posix_memalign(&m, 128, nr_comp * size * sizeof(fields_real_t));
  assert(ierr == 0);
  pf->flds = m; 
#else
  pf->data = calloc(pf->nr_comp * size, sizeof(fields_real_t));
#endif
}

// ----------------------------------------------------------------------
// psc_fields_destroy

static void
PFX(destroy)(struct psc_fields *pf)
{
#if PSC_FIELDS_AS_FORTRAN
  fields_real_t **flds = pf->data;
  free(flds[0]);

  for (int i = 0; i < pf->nr_comp; i++) {
    flds[i] = NULL;
  }
  free(flds);
#else
  free(pf->data);
#endif
}

// ----------------------------------------------------------------------
// psc_fields_zero_comp

static void
PFX(zero_comp)(struct psc_fields *pf, int m)
{
  memset(&F3(pf, m, pf->ib[0], pf->ib[1], pf->ib[2]), 0,
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_real_t));
}

// ----------------------------------------------------------------------
// psc_fields_set_comp

static void
PFX(set_comp)(struct psc_fields *pf, int m, double _val)
{
  fields_real_t val = _val;

  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3(pf, m, jx, jy, jz) = val;
      }
    }
  }
}

// ----------------------------------------------------------------------
// psc_fields_scale_comp

static void
PFX(scale_comp)(struct psc_fields *pf, int m, double _val)
{
  fields_real_t val = _val;

  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3(pf, m, jx, jy, jz) *= val;
      }
    }
  }
}

// ----------------------------------------------------------------------
// psc_fields_copy_comp

static void
PFX(copy_comp)(struct psc_fields *pto, int m_to, struct psc_fields *pfrom, int m_from)
{
  for (int jz = pto->ib[2]; jz < pto->ib[2] + pto->im[2]; jz++) {
    for (int jy = pto->ib[1]; jy < pto->ib[1] + pto->im[1]; jy++) {
      for (int jx = pto->ib[0]; jx < pto->ib[0] + pto->im[0]; jx++) {
	F3(pto, m_to, jx, jy, jz) = F3(pfrom, m_from, jx, jy, jz);
      }
    }
  }
}

// ----------------------------------------------------------------------
// psc_fields_axpy_comp

static void
PFX(axpy_comp)(struct psc_fields *y, int ym, double _a, struct psc_fields *x, int xm)
{
  fields_real_t a = _a;

  for (int jz = y->ib[2]; jz < y->ib[2] + y->im[2]; jz++) {
    for (int jy = y->ib[1]; jy < y->ib[1] + y->im[1]; jy++) {
      for (int jx = y->ib[0]; jx < y->ib[0] + y->im[0]; jx++) {
	F3(y, ym, jx, jy, jz) += a * F3(x, xm, jx, jy, jz);
      }
    }
  }
}

#if defined(HAVE_LIBHDF5_HL) && (PSC_FIELDS_AS_SINGLE || PSC_FIELDS_AS_C)

#include <mrc_io.h>

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------
// psc_fields_write

static void
PFX(write)(struct psc_fields *flds, struct mrc_io *io)
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
#if PSC_FIELDS_AS_SINGLE
  ierr = H5LTmake_dataset_float(group, "fields_single", 4, hdims, flds->data); CE;
#elif PSC_FIELDS_AS_C
  ierr = H5LTmake_dataset_double(group, "fields_c", 4, hdims, flds->data); CE;
#endif
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_fields_read

static void
PFX(read)(struct psc_fields *flds, struct mrc_io *io)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, flds), H5P_DEFAULT); H5_CHK(group);
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
  psc_fields_setup(flds);
#if PSC_FIELDS_AS_SINGLE
  ierr = H5LTread_dataset_float(group, "fields_single", flds->data); CE;
#elif PSC_FIELDS_AS_C
  ierr = H5LTread_dataset_double(group, "fields_c", flds->data); CE;
#endif
  ierr = H5Gclose(group); CE;
}

#endif // HAVE_LIBHDF5_HL

// ----------------------------------------------------------------------
// psc_fields: subclass ops
  
struct psc_fields_ops PFX(ops) = {
  .name                  = FIELDS_TYPE,
  .methods               = PFX(methods),
  .setup                 = PFX(setup),
  .destroy               = PFX(destroy),
#if defined(HAVE_LIBHDF5_HL) && (PSC_FIELDS_AS_SINGLE || PSC_FIELDS_AS_C)
  .write                 = PFX(write),
  .read                  = PFX(read),
#endif
  .zero_comp             = PFX(zero_comp),
};

// ----------------------------------------------------------------------
// psc_mfields_zero_comp

static void
MPFX(zero_comp)(struct psc_mfields *mflds, int m)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    PFX(zero_comp)(psc_mfields_get_patch(mflds, p), m);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_set_comp

static void
MPFX(set_comp)(struct psc_mfields *mflds, int m, double alpha)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    PFX(set_comp)(psc_mfields_get_patch(mflds, p), m, alpha);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_scale_comp

static void
MPFX(scale_comp)(struct psc_mfields *mflds, int m, double alpha)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    PFX(scale_comp)(psc_mfields_get_patch(mflds, p), m, alpha);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_copy_comp

static void
MPFX(copy_comp)(struct psc_mfields *to, int mto,
		struct psc_mfields *fr, int mfr)
{
  for (int p = 0; p < to->nr_patches; p++) {
    PFX(copy_comp)(psc_mfields_get_patch(to, p), mto,
		   psc_mfields_get_patch(fr, p), mfr);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_axpy_comp

static void
MPFX(axpy_comp)(struct psc_mfields *y, int my, double alpha,
		struct psc_mfields *x, int mx)
{
  for (int p = 0; p < y->nr_patches; p++) {
    PFX(axpy_comp)(psc_mfields_get_patch(y, p), my, alpha,
		   psc_mfields_get_patch(x, p), mx);
  }
}

// ----------------------------------------------------------------------
// psc_mfields: subclass ops
  
struct psc_mfields_ops MPFX(ops) = {
  .name                  = FIELDS_TYPE,
  .zero_comp             = MPFX(zero_comp),
  .set_comp              = MPFX(set_comp),
  .scale_comp            = MPFX(scale_comp),
  .copy_comp             = MPFX(copy_comp),
  .axpy_comp             = MPFX(axpy_comp),
};

