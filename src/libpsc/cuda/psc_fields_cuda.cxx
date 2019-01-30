
#include "psc.h"
#include "cuda_iface.h"
#include "cuda_iface_bnd.h"
#include "psc_fields_cuda.h"
#include "psc_fields_c.h"
#include "psc_fields_single.h"
#include "fields.hxx"

#include <mrc_params.h>

using FieldsH = Fields3d<fields_single_t>; // host
using FieldsS = Fields3d<fields_single_t> ;// host
using FieldsC = Fields3d<fields_c_t>; // host

// OPT, CUDA fields have too many ghostpoints, and 7 points in the invar direction!

// ======================================================================
// convert from/to "c"

static void psc_mfields_cuda_copy_from_c(MfieldsBase& mflds_cuda, MfieldsBase& mflds_c, int mb, int me)
{
  auto& mf_cuda = dynamic_cast<MfieldsCuda&>(mflds_cuda);
  auto& mf_c = dynamic_cast<MfieldsC&>(mflds_c);
  fields_single_t flds = mf_cuda.get_host_fields();
  FieldsH F(flds);

  for (int p = 0; p < mf_cuda.n_patches(); p++) {
    FieldsC F_c(mf_c[p]);
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib_[2]; jz < flds.ib_[2] + flds.im_[2]; jz++) {
	for (int jy = flds.ib_[1]; jy < flds.ib_[1] + flds.im_[1]; jy++) {
	  for (int jx = flds.ib_[0]; jx < flds.ib_[0] + flds.im_[0]; jx++) {
	    F(m, jx,jy,jz) = F_c( m, jx,jy,jz);
	  }
	}
      }
    }

    mf_cuda.copy_to_device(p, flds, mb, me);
  }
  
  flds.dtor();
}

static void psc_mfields_cuda_copy_to_c(MfieldsBase& mflds_cuda, MfieldsBase& mflds_c, int mb, int me)
{
  auto& mf_cuda = dynamic_cast<MfieldsCuda&>(mflds_cuda);
  auto& mf_c = dynamic_cast<MfieldsC&>(mflds_c);
  fields_single_t flds = mf_cuda.get_host_fields();
  FieldsH F(flds);

  for (int p = 0; p < mf_cuda.n_patches(); p++) {
    FieldsC F_c(mf_c[p]);
    mf_cuda.copy_from_device(p, flds, mb, me);
  
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib_[2]; jz < flds.ib_[2] + flds.im_[2]; jz++) {
	for (int jy = flds.ib_[1]; jy < flds.ib_[1] + flds.im_[1]; jy++) {
	  for (int jx = flds.ib_[0]; jx < flds.ib_[0] + flds.im_[0]; jx++) {
	    F_c(m, jx,jy,jz) = F(m, jx,jy,jz);
	  }
	}
      }
    }
  }

  flds.dtor();
}

static void psc_mfields_state_cuda_copy_from_c(MfieldsStateBase& mflds_cuda, MfieldsStateBase& mflds_c, int mb, int me)
{
  auto& mf_cuda = dynamic_cast<MfieldsStateCuda&>(mflds_cuda);
  auto& mf_c = dynamic_cast<MfieldsStateDouble&>(mflds_c);
  fields_single_t flds = mf_cuda.get_host_fields();
  FieldsH F(flds);

  for (int p = 0; p < mf_cuda.n_patches(); p++) {
    FieldsC F_c(mf_c[p]);
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib_[2]; jz < flds.ib_[2] + flds.im_[2]; jz++) {
	for (int jy = flds.ib_[1]; jy < flds.ib_[1] + flds.im_[1]; jy++) {
	  for (int jx = flds.ib_[0]; jx < flds.ib_[0] + flds.im_[0]; jx++) {
	    F(m, jx,jy,jz) = F_c( m, jx,jy,jz);
	  }
	}
      }
    }

    mf_cuda.copy_to_device(p, flds, mb, me);
  }
  
  flds.dtor();
}

static void psc_mfields_state_cuda_copy_to_c(MfieldsStateBase& mflds_cuda, MfieldsStateBase& mflds_c, int mb, int me)
{
  auto& mf_cuda = dynamic_cast<MfieldsStateCuda&>(mflds_cuda);
  auto& mf_c = dynamic_cast<MfieldsStateDouble&>(mflds_c);
  fields_single_t flds = mf_cuda.get_host_fields();
  FieldsH F(flds);

  for (int p = 0; p < mf_cuda.n_patches(); p++) {
    FieldsC F_c(mf_c[p]);
    mf_cuda.copy_from_device(p, flds, mb, me);
  
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib_[2]; jz < flds.ib_[2] + flds.im_[2]; jz++) {
	for (int jy = flds.ib_[1]; jy < flds.ib_[1] + flds.im_[1]; jy++) {
	  for (int jx = flds.ib_[0]; jx < flds.ib_[0] + flds.im_[0]; jx++) {
	    F_c(m, jx,jy,jz) = F(m, jx,jy,jz);
	  }
	}
      }
    }
  }

  flds.dtor();
}

// ======================================================================
// convert from/to "single"

static void psc_mfields_cuda_copy_from_single(MfieldsBase& mflds_cuda, MfieldsBase& mflds_single, int mb, int me)
{
  auto& mf_cuda = dynamic_cast<MfieldsCuda&>(mflds_cuda);
  auto& mf_single = dynamic_cast<MfieldsSingle&>(mflds_single);
  fields_single_t flds = mf_cuda.get_host_fields();
  FieldsH F(flds);
  
  for (int p = 0; p < mf_cuda.n_patches(); p++) {
    FieldsS F_s(mf_single[p]);

    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib_[2]; jz < flds.ib_[2] + flds.im_[2]; jz++) {
	for (int jy = flds.ib_[1]; jy < flds.ib_[1] + flds.im_[1]; jy++) {
	  for (int jx = flds.ib_[0]; jx < flds.ib_[0] + flds.im_[0]; jx++) {
	    F(m, jx,jy,jz) = F_s(m, jx,jy,jz);
	  }
	}
      }
    }

    mf_cuda.copy_to_device(p, flds, mb, me);
  }
  
  flds.dtor();
}

static void psc_mfields_state_cuda_copy_from_single(MfieldsStateBase& mflds_cuda, MfieldsStateBase& mflds_single, int mb, int me)
{
  auto& mf_cuda = dynamic_cast<MfieldsStateCuda&>(mflds_cuda);
  auto& mf_single = dynamic_cast<MfieldsStateSingle&>(mflds_single);
  fields_single_t flds = mf_cuda.get_host_fields();
  FieldsH F(flds);
  
  for (int p = 0; p < mf_cuda.n_patches(); p++) {
    FieldsS F_s(mf_single[p]);

    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib_[2]; jz < flds.ib_[2] + flds.im_[2]; jz++) {
	for (int jy = flds.ib_[1]; jy < flds.ib_[1] + flds.im_[1]; jy++) {
	  for (int jx = flds.ib_[0]; jx < flds.ib_[0] + flds.im_[0]; jx++) {
	    F(m, jx,jy,jz) = F_s(m, jx,jy,jz);
	  }
	}
      }
    }

    mf_cuda.copy_to_device(p, flds, mb, me);
  }
  
  flds.dtor();
}

static void psc_mfields_cuda_copy_to_single(MfieldsBase& mflds_cuda, MfieldsBase& mflds_single, int mb, int me)
{
  auto& mf_cuda = dynamic_cast<MfieldsCuda&>(mflds_cuda);
  auto& mf_single = dynamic_cast<MfieldsSingle&>(mflds_single);
  fields_single_t flds = mf_cuda.get_host_fields();
  FieldsH F(flds);

  for (int p = 0; p < mf_cuda.n_patches(); p++) {
    FieldsS F_s(mf_single[p]);
    mf_cuda.copy_from_device(p, flds, mb, me);
  
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib_[2]; jz < flds.ib_[2] + flds.im_[2]; jz++) {
	for (int jy = flds.ib_[1]; jy < flds.ib_[1] + flds.im_[1]; jy++) {
	  for (int jx = flds.ib_[0]; jx < flds.ib_[0] + flds.im_[0]; jx++) {
	    F_s(m, jx,jy,jz) = F(m, jx,jy,jz);
	  }
	}
      }
    }
  }

  flds.dtor();
}

static void psc_mfields_state_cuda_copy_to_single(MfieldsStateBase& mflds_cuda, MfieldsStateBase& mflds_single, int mb, int me)
{
  auto& mf_cuda = dynamic_cast<MfieldsStateCuda&>(mflds_cuda);
  auto& mf_single = dynamic_cast<MfieldsStateSingle&>(mflds_single);
  fields_single_t flds = mf_cuda.get_host_fields();
  FieldsH F(flds);

  for (int p = 0; p < mf_cuda.n_patches(); p++) {
    FieldsS F_s(mf_single[p]);
    mf_cuda.copy_from_device(p, flds, mb, me);
  
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib_[2]; jz < flds.ib_[2] + flds.im_[2]; jz++) {
	for (int jy = flds.ib_[1]; jy < flds.ib_[1] + flds.im_[1]; jy++) {
	  for (int jx = flds.ib_[0]; jx < flds.ib_[0] + flds.im_[0]; jx++) {
	    F_s(m, jx,jy,jz) = F(m, jx,jy,jz);
	  }
	}
      }
    }
  }

  flds.dtor();
}

// ======================================================================

#ifdef HAVE_LIBHDF5_HL

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
psc_mfields_cuda_write(MfieldsCuda& mflds, struct mrc_io *io)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group0 = H5Gopen(h5_file, mrc_io_obj_path(io, _mflds), H5P_DEFAULT); H5_CHK(group0);

  fields_single_t flds = mflds->get_host_fields();

  for (int p = 0; p < mflds.n_patches(); p++) {
    mflds->copy_from_device(p, flds, 0, flds.nr_comp);
    char name[20]; sprintf(name, "flds%d", p);
    hid_t group = H5Gcreate(group0, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
    
    ierr = H5LTset_attribute_int(group, ".", "ib", flds.ib, 3); CE;
    ierr = H5LTset_attribute_int(group, ".", "im", flds.im, 3); CE;
    ierr = H5LTset_attribute_int(group, ".", "nr_comp", &flds.nr_comp, 1); CE;
    // write components separately instead?
    hsize_t hdims[4];
    hdims[0] = flds.nr_comp;
    hdims[1] = flds.im_[2]; hdims[2] = flds.im_[1]; hdims[3] = flds.im_[0];
    ierr = H5LTmake_dataset_float(group, "fields_cuda", 4, hdims, flds.data); CE;
    ierr = H5Gclose(group); CE;
  }
  flds.dtor();

  ierr = H5Gclose(group0); CE;
}

// ----------------------------------------------------------------------
// psc_mfields_cuda_read

static void
psc_mfields_cuda_read(MfieldsCuda& mflds, struct mrc_io *io)
{
  psc_mfields_read_super(_mflds, io);
  
  psc_mfields_cuda_setup(_mflds);

  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group0 = H5Gopen(h5_file, mrc_io_obj_path(io, _mflds), H5P_DEFAULT); H5_CHK(group0);

  fields_single_t flds = mflds->get_host_fields();
  for (int p = 0; p < mflds.n_patches(); p++) {
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

    ierr = H5LTread_dataset_float(group, "fields_cuda", flds.data); CE;
    mflds->copy_to_device(p, flds, 0, flds.nr_comp);
    ierr = H5Gclose(group); CE;
  }
  flds.dtor();
  ierr = H5Gclose(group0); CE;
}

#endif

// ======================================================================
  
const MfieldsBase::Convert MfieldsCuda::convert_to_ = {
  { std::type_index(typeid(MfieldsC))     , psc_mfields_cuda_copy_to_c },
  { std::type_index(typeid(MfieldsSingle)), psc_mfields_cuda_copy_to_single },
};

const MfieldsBase::Convert MfieldsCuda::convert_from_ = {
  { std::type_index(typeid(MfieldsC))     , psc_mfields_cuda_copy_from_c },
  { std::type_index(typeid(MfieldsSingle)), psc_mfields_cuda_copy_from_single },
};

// ======================================================================
  
const MfieldsStateBase::Convert MfieldsStateCuda::convert_to_ = {
  { std::type_index(typeid(MfieldsStateDouble)), psc_mfields_state_cuda_copy_to_c },
  { std::type_index(typeid(MfieldsStateSingle)), psc_mfields_state_cuda_copy_to_single },
};

const MfieldsStateBase::Convert MfieldsStateCuda::convert_from_ = {
  { std::type_index(typeid(MfieldsStateDouble)), psc_mfields_state_cuda_copy_from_c },
  { std::type_index(typeid(MfieldsStateSingle)), psc_mfields_state_cuda_copy_from_single },
};

