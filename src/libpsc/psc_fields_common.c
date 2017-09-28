
// ----------------------------------------------------------------------
// fields_t_zero_comp

static inline void
fields_t_zero_comp(fields_t flds, int m)
{
  for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
    for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
      for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	_F3(flds, m, jx,jy,jz) = 0;
      }
    }
  }
}

// ----------------------------------------------------------------------
// fields_t_set_comp

static inline void
fields_t_set_comp(fields_t flds, int m, fields_real_t val)
{
  for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
    for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
      for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	_F3(flds, m, jx,jy,jz) = val;
      }
    }
  }
}

// ----------------------------------------------------------------------
// fields_t_scale_comp

static inline void
fields_t_scale_comp(fields_t flds, int m, fields_real_t val)
{
  for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
    for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
      for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	_F3(flds, m, jx, jy, jz) *= val;
      }
    }
  }
}

// ----------------------------------------------------------------------
// fields_t_copy_comp

static inline void
fields_t_copy_comp(fields_t to, int m_to, fields_t from, int m_from)
{
  for (int jz = to.ib[2]; jz < to.ib[2] + to.im[2]; jz++) {
    for (int jy = to.ib[1]; jy < to.ib[1] + to.im[1]; jy++) {
      for (int jx = to.ib[0]; jx < to.ib[0] + to.im[0]; jx++) {
	_F3(to, m_to, jx,jy,jz) = _F3(from, m_from, jx,jy,jz);
      }
    }
  }
}

// ----------------------------------------------------------------------
// fields_t_axpy_comp

static inline void
fields_t_axpy_comp(fields_t y, int m_y, fields_real_t a, fields_t x, int m_x)
{
  for (int jz = y.ib[2]; jz < y.ib[2] + y.im[2]; jz++) {
    for (int jy = y.ib[1]; jy < y.ib[1] + y.im[1]; jy++) {
      for (int jx = y.ib[0]; jx < y.ib[0] + y.im[0]; jx++) {
	_F3(y, m_y, jx,jy,jz) += a * _F3(x, m_x, jx,jy,jz);
      }
    }
  }
}


// ======================================================================
// psc_mfields

struct MPFX(sub) {
};

// ----------------------------------------------------------------------
// psc_mfields_setup

static void
MPFX(setup)(struct psc_mfields *mflds)
{
  psc_mfields_setup_super(mflds);

  mflds->data = calloc(mflds->nr_patches, sizeof(*mflds->data));
  for (int p = 0; p < mflds->nr_patches; p++) {
    unsigned int size = 1;
    for (int d = 0; d < 3; d++) {
      size *= mflds->im[d];
    }
    
#if PSC_FIELDS_AS_FORTRAN
    fields_real_t **flds = calloc(mflds->nr_fields, sizeof(*flds));
    flds[0] = calloc(size * mflds->nr_fields, sizeof(flds[0]));
    for (int i = 1; i < mflds->nr_fields; i++) {
      flds[i] = flds[0] + i * size;
    }
    mflds->data[p] = flds;
#else
    mflds->data[p] = calloc(mflds->nr_fields * size, sizeof(fields_real_t));
#endif
  }
}

// ----------------------------------------------------------------------
// psc_mfields_destroy

static void
MPFX(destroy)(struct psc_mfields *mflds)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
#if PSC_FIELDS_AS_FORTRAN
    fields_real_t **flds = mflds->data[p];
    free(flds[0]);
    
    for (int i = 0; i < mflds->nr_fields; i++) {
      flds[i] = NULL;
    }
    free(flds);
#else
    free(mflds->data[p]);
#endif
  }
  free(mflds->data);
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
// psc_mfields_write

static void
MPFX(write)(struct psc_mfields *mflds, struct mrc_io *io)
{
  herr_t ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group0 = H5Gopen(h5_file, mrc_io_obj_path(io, mflds), H5P_DEFAULT); H5_CHK(group0);

  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t flds = fields_t_mflds(mflds, p);
    char name[20]; sprintf(name, "flds%d", p);
    hid_t group = H5Gcreate(group0, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
    ierr = H5LTset_attribute_int(group, ".", "ib", flds.ib, 3); CE;
    ierr = H5LTset_attribute_int(group, ".", "im", flds.im, 3); CE;
    ierr = H5LTset_attribute_int(group, ".", "nr_comp", &flds.nr_comp, 1); CE;
    // write components separately instead?
    hsize_t hdims[4] = { flds.nr_comp, flds.im[2], flds.im[1], flds.im[0] };
#if PSC_FIELDS_AS_SINGLE
    ierr = H5LTmake_dataset_float(group, "fields_single", 4, hdims, flds.data); CE;
#elif PSC_FIELDS_AS_C
    ierr = H5LTmake_dataset_double(group, "fields_c", 4, hdims, flds.data); CE;
#endif
    ierr = H5Gclose(group); CE;
  }

  ierr = H5Gclose(group0); CE;
}

// ----------------------------------------------------------------------
// psc_mfields_read

static void
MPFX(read)(struct psc_mfields *mflds, struct mrc_io *io)
{
  psc_mfields_read_super(mflds, io);

  psc_mfields_setup(mflds);

  herr_t ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group0 = H5Gopen(h5_file, mrc_io_obj_path(io, mflds), H5P_DEFAULT); H5_CHK(group0);

  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t flds = fields_t_mflds(mflds, p);
    char name[20]; sprintf(name, "flds%d", p);
    hid_t group = H5Gopen(group0, name, H5P_DEFAULT); H5_CHK(group);
    int ib[3], im[3], nr_comp;
    ierr = H5LTget_attribute_int(group, ".", "ib", ib); CE;
    ierr = H5LTget_attribute_int(group, ".", "im", im); CE;
    ierr = H5LTget_attribute_int(group, ".", "nr_comp", &nr_comp); CE;
    for (int d = 0; d < 3; d++) {
      assert(ib[d] == flds.ib[d]);
      assert(im[d] == flds.im[d]);
    }
    assert(nr_comp == flds.nr_comp);
#if PSC_FIELDS_AS_SINGLE
    ierr = H5LTread_dataset_float(group, "fields_single", flds.data); CE;
#elif PSC_FIELDS_AS_C
    ierr = H5LTread_dataset_double(group, "fields_c", flds.data); CE;
#endif
    ierr = H5Gclose(group); CE;
  }

  ierr = H5Gclose(group0); CE;
}

#endif

// ----------------------------------------------------------------------
// psc_mfields_zero_comp

static void
MPFX(zero_comp)(struct psc_mfields *mflds, int m)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t_zero_comp(fields_t_mflds(mflds, p), m);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_set_comp

static void
MPFX(set_comp)(struct psc_mfields *mflds, int m, double alpha)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t_set_comp(fields_t_mflds(mflds, p), m, alpha);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_scale_comp

static void
MPFX(scale_comp)(struct psc_mfields *mflds, int m, double alpha)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t_scale_comp(fields_t_mflds(mflds, p), m, alpha);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_copy_comp

static void
MPFX(copy_comp)(struct psc_mfields *to, int mto,
		struct psc_mfields *fr, int mfr)
{
  for (int p = 0; p < to->nr_patches; p++) {
    fields_t_copy_comp(fields_t_mflds(to, p), mto,
		       fields_t_mflds(fr, p), mfr);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_axpy_comp

static void
MPFX(axpy_comp)(struct psc_mfields *y, int my, double alpha,
		struct psc_mfields *x, int mx)
{
  for (int p = 0; p < y->nr_patches; p++) {
    fields_t_axpy_comp(fields_t_mflds(y, p), my, alpha,
		       fields_t_mflds(x, p), mx);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_get_field_t

fields_t
MPFX(get_field_t)(struct psc_mfields *mflds, int p)
{
  fields_t flds;

  flds.data = (fields_real_t *) mflds->data[p];
  for (int d = 0; d < 3; d++) {
    flds.ib[d] = mflds->ib[d];
    flds.im[d] = mflds->im[d];
  }
  flds.nr_comp = mflds->nr_fields;
  flds.first_comp = mflds->first_comp;

  return flds;
}

// ----------------------------------------------------------------------
// psc_mfields: subclass ops
  
struct psc_mfields_ops MPFX(ops) = {
  .name                  = FIELDS_TYPE,
  .size                  = sizeof(struct MPFX(sub)),
  .methods               = MPFX(methods),
  .setup                 = MPFX(setup),
  .destroy               = MPFX(destroy),
#if defined(HAVE_LIBHDF5_HL) && (PSC_FIELDS_AS_SINGLE || PSC_FIELDS_AS_C)
  .write                 = MPFX(write),
  .read                  = MPFX(read),
#endif
  .zero_comp             = MPFX(zero_comp),
  .set_comp              = MPFX(set_comp),
  .scale_comp            = MPFX(scale_comp),
  .copy_comp             = MPFX(copy_comp),
  .axpy_comp             = MPFX(axpy_comp),
};

