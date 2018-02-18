
#include "fields.hxx"

#include <limits>
#include <algorithm>
#include <math.h>

using Fields = Fields3d<fields_t>;

// ----------------------------------------------------------------------
// fields_t_zero_comp

static inline void
fields_t_zero_comp(fields_t flds, int m)
{
  Fields F(flds);
  for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
    for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
      for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	F(m, jx,jy,jz) = 0;
      }
    }
  }
}

// ----------------------------------------------------------------------
// fields_t_set_comp

static inline void
fields_t_set_comp(fields_t flds, int m, fields_t::real_t val)
{
  Fields F(flds);
  for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
    for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
      for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	F(m, jx,jy,jz) = val;
      }
    }
  }
}

// ----------------------------------------------------------------------
// fields_t_scale_comp

static inline void
fields_t_scale_comp(fields_t flds, int m, fields_t::real_t val)
{
  Fields F(flds);
  for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
    for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
      for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	F(m, jx, jy, jz) *= val;
      }
    }
  }
}

// ----------------------------------------------------------------------
// fields_t_copy_comp

static inline void
fields_t_copy_comp(fields_t to, int m_to, fields_t from, int m_from)
{
  Fields F_to(to), F_from(from);
  for (int jz = to.ib[2]; jz < to.ib[2] + to.im[2]; jz++) {
    for (int jy = to.ib[1]; jy < to.ib[1] + to.im[1]; jy++) {
      for (int jx = to.ib[0]; jx < to.ib[0] + to.im[0]; jx++) {
	F_to(m_to, jx,jy,jz) = F_from(m_from, jx,jy,jz);
      }
    }
  }
}

// ----------------------------------------------------------------------
// fields_t_axpy_comp

static inline void
fields_t_axpy_comp(fields_t y, int m_y, fields_t::real_t a, fields_t x, int m_x)
{
  Fields X(x), Y(y);
  for (int jz = y.ib[2]; jz < y.ib[2] + y.im[2]; jz++) {
    for (int jy = y.ib[1]; jy < y.ib[1] + y.im[1]; jy++) {
      for (int jx = y.ib[0]; jx < y.ib[0] + y.im[0]; jx++) {
	Y(m_y, jx,jy,jz) += a * X(m_x, jx,jy,jz);
      }
    }
  }
}


// ----------------------------------------------------------------------
// fields_t_max_comp

static inline double
fields_t_max_comp(fields_t f, int m)
{
  Fields F(f);
  fields_t::real_t rv = std::numeric_limits<fields_t::real_t>::min();
  for (int jz = f.ib[2]; jz < f.ib[2] + f.im[2]; jz++) {
    for (int jy = f.ib[1]; jy < f.ib[1] + f.im[1]; jy++) {
      for (int jx = f.ib[0]; jx < f.ib[0] + f.im[0]; jx++) {
	rv = std::max(rv, F(m, jx,jy,jz));
      }
    }
  }

  return rv;
}


// ======================================================================
// psc_mfields

// ----------------------------------------------------------------------
// psc_mfields_setup

static void
MPFX(setup)(struct psc_mfields *_mflds)
{
  mfields_t mflds(_mflds);

  psc_mfields_setup_super(_mflds);

  new(mflds.sub()) MPFX(sub)(ppsc->grid, mflds.n_fields(), _mflds->ibn);
}

// ----------------------------------------------------------------------
// psc_mfields_destroy

static void
MPFX(destroy)(struct psc_mfields *_mflds)
{
  mfields_t mflds(_mflds);

  mflds->~MPFX(sub)();
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
  mfields_t mf(mflds);
  herr_t ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group0 = H5Gopen(h5_file, mrc_io_obj_path(io, mflds), H5P_DEFAULT); H5_CHK(group0);

  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t flds = mf[p];
    char name[20]; sprintf(name, "flds%d", p);
    hid_t group = H5Gcreate(group0, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
    ierr = H5LTset_attribute_int(group, ".", "ib", flds.ib, 3); CE;
    ierr = H5LTset_attribute_int(group, ".", "im", flds.im, 3); CE;
    ierr = H5LTset_attribute_int(group, ".", "nr_comp", &flds.nr_comp, 1); CE;
    // write components separately instead?
    hsize_t hdims[4];
    hdims[0] = flds.nr_comp;
    hdims[1] = flds.im[2];
    hdims[2] = flds.im[1];
    hdims[3] = flds.im[0];
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

  mfields_t mf(mflds);

  herr_t ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group0 = H5Gopen(h5_file, mrc_io_obj_path(io, mflds), H5P_DEFAULT); H5_CHK(group0);

  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t flds = mf[p];
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
  mfields_t mf(mflds);
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t_zero_comp(mf[p], m);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_set_comp

static void
MPFX(set_comp)(struct psc_mfields *mflds, int m, double alpha)
{
  mfields_t mf(mflds);
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t_set_comp(mf[p], m, alpha);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_scale_comp

static void
MPFX(scale_comp)(struct psc_mfields *mflds, int m, double alpha)
{
  mfields_t mf(mflds);
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t_scale_comp(mf[p], m, alpha);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_copy_comp

static void
MPFX(copy_comp)(struct psc_mfields *to, int mto,
		struct psc_mfields *fr, int mfr)
{
  mfields_t mf_to(to), mf_from(fr);
  for (int p = 0; p < to->nr_patches; p++) {
    fields_t_copy_comp(mf_to[p], mto, mf_from[p], mfr);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_axpy_comp

static void
MPFX(axpy_comp)(struct psc_mfields *y, int my, double alpha,
		struct psc_mfields *x, int mx)
{
  mfields_t mf_y(y), mf_x(x);
  for (int p = 0; p < y->nr_patches; p++) {
    fields_t_axpy_comp(mf_y[p], my, alpha, mf_x[p], mx);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_max_comp

static double
MPFX(max_comp)(struct psc_mfields *mflds, int m)
{
  mfields_t mf(mflds);
  double rv = -1e37; // FIXME
  for (int p = 0; p < mflds->nr_patches; p++) {
    rv = fmax(rv, fields_t_max_comp(mf[p], m));
  }
  return rv;
}

// ----------------------------------------------------------------------
// psc_mfields_get_field_t

fields_t
MPFX(get_field_t)(struct psc_mfields *mflds, int p)
{
  //assert((struct psc_mfields_ops *) mflds->obj.ops == &MPFX(ops));
  MPFX(sub) *sub = mrc_to_subobj(mflds, MPFX(sub));
  fields_t flds;

#if PSC_FIELDS_AS_FORTRAN
  flds.data = sub->data[p][0];
#else
  flds.data = sub->data[p];
#endif
  for (int d = 0; d < 3; d++) {
    flds.ib[d] = sub->ib[d];
    flds.im[d] = sub->im[d];
  }
  flds.nr_comp = mflds->nr_fields;
  flds.first_comp = mflds->first_comp;

  return flds;
}

// ----------------------------------------------------------------------
// psc_mfields: subclass ops
  
struct MPFX(psc_mfields_ops) : psc_mfields_ops {
  MPFX(psc_mfields_ops)() {
    name                  = FIELDS_TYPE;
    size                  = sizeof(MPFX(sub));
    methods               = MPFX(methods);
    setup                 = MPFX(setup);
    destroy               = MPFX(destroy);
#if defined(HAVE_LIBHDF5_HL) && (PSC_FIELDS_AS_SINGLE || PSC_FIELDS_AS_C)
    write                 = MPFX(write);
    read                  = MPFX(read);
#endif
    zero_comp             = MPFX(zero_comp);
    set_comp              = MPFX(set_comp);
    scale_comp            = MPFX(scale_comp);
    copy_comp             = MPFX(copy_comp);
    axpy_comp             = MPFX(axpy_comp);
    max_comp              = MPFX(max_comp);
  }
} MPFX(ops);

