
#include <mrc_io_private.h>
#include <mrc_params.h>
#include <mrc_list.h>
#include "mrc_io_xdmf_lib.h"

#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>


#ifdef H5_HAVE_PARALLEL
#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------

struct ph5 {
  hid_t h5_file;
  int sw;
  // parallel only
  bool use_independent_io;
  char *separator; // seperator between basename & numbers
};

#define VAR(x) (void *)offsetof(struct ph5, x)
static struct param ph5_descr[] = {
  { "sw"                     , VAR(sw)                      , PARAM_INT(0)           },
  { "independent"            , VAR(use_independent_io)      , PARAM_BOOL(true)       },
  { "separator"              , VAR(separator)               , PARAM_STRING(".")      },
  {},
};
#undef VAR

#define to_ph5(io) ((struct ph5 *)((io)->obj.subctx))

// ======================================================================
// hdf5_parallel


static void
phdf5_write_attr(struct mrc_io *io, const char *path, int type,
		const char *name, union param_u *pv)
{
  // the HDF5 docs claim attribute writing calls can be done collectively.
  // Not sure I believe them, but let's see...
  struct ph5 *ph5 = to_ph5(io);
  
  hid_t group;
  if (H5Lexists(ph5->h5_file, path, H5P_DEFAULT) > 0) {
    group = H5Gopen(ph5->h5_file, path, H5P_DEFAULT); H5_CHK(group);
  } else {
    group = H5Gcreate(ph5->h5_file, path, H5P_DEFAULT, H5P_DEFAULT,
		      H5P_DEFAULT); H5_CHK(group);
  }

  int ierr;

  switch (type) {
  case PT_SELECT:
  case PT_INT:
  case MRC_VAR_INT:
    ierr = H5LTset_attribute_int(group, ".", name, &pv->u_int, 1); CE;
    break;
  case PT_BOOL: 
  case MRC_VAR_BOOL: {
    int val = pv->u_bool;
    ierr = H5LTset_attribute_int(group, ".", name, &val, 1); CE;
    break;
  }
  case PT_FLOAT:
  case MRC_VAR_FLOAT:
    ierr = H5LTset_attribute_float(group, ".", name, &pv->u_float, 1); CE;
    break;
  case PT_DOUBLE:
  case MRC_VAR_DOUBLE:
    ierr = H5LTset_attribute_double(group, ".", name, &pv->u_double, 1); CE;
    break;
  case PT_STRING:
    if (pv->u_string) {
      ierr = H5LTset_attribute_string(group, ".", name, pv->u_string); CE;
    } else {
      ierr = H5LTset_attribute_string(group, ".", name, "(NULL)"); CE;
    }
    break;

  case PT_INT3:
    ierr = H5LTset_attribute_int(group, ".", name, pv->u_int3, 3); CE;
    break;

  case PT_FLOAT3:
    ierr = H5LTset_attribute_float(group, ".", name, pv->u_float3, 3); CE;
    break;

  case PT_DOUBLE3:
    ierr = H5LTset_attribute_double(group, ".", name, pv->u_double3, 3); CE;
    break;
  case PT_INT_ARRAY:
    ierr = H5LTset_attribute_int(group, ".", name, pv->u_int_array.vals,
				 pv->u_int_array.nr_vals); CE;
    break;
  case PT_PTR:
    break;
  default:
    mpi_printf(mrc_io_comm(io), "mrc_io_hdf5_parallel: not writing attr '%s' (type %d)\n",
	    name, type);
    assert(0);
  }
  ierr = H5Gclose(group); CE;

}

static void
phdf5_read_attr(struct mrc_io *io, const char *path, int type,
		const char *name, union param_u *pv)
{
  // the HDF5 docs claim attribute writing calls can be done collectively.
  // Not sure I believe them, but let's see...
  struct ph5 *ph5 = to_ph5(io);
  int ierr;
  hid_t group = H5Gopen(ph5->h5_file, path, H5P_DEFAULT); H5_CHK(group);
  switch (type) {
  case PT_SELECT:
  case PT_INT:
  case MRC_VAR_INT:
    ierr = H5LTget_attribute_int(group, ".", name, &pv->u_int); CE;
    break;
  case PT_BOOL:
  case MRC_VAR_BOOL: ;
    int val;
    ierr = H5LTget_attribute_int(group, ".", name, &val); CE;
    pv->u_bool = val;
    break;
  case PT_FLOAT:
    ierr = H5LTget_attribute_float(group, ".", name, &pv->u_float); CE;
    break;
  case PT_DOUBLE:
    ierr = H5LTget_attribute_double(group, ".", name, &pv->u_double); CE;
    break;
  case PT_STRING: ;
    hsize_t dims;
    H5T_class_t class;
    size_t sz;
    ierr = H5LTget_attribute_info(group, ".", name, &dims, &class, &sz); CE;
    pv->u_string = malloc(sz);
    ierr = H5LTget_attribute_string(group, ".", name, (char *)pv->u_string); CE;
    if (strcmp(pv->u_string, "(NULL)") == 0) {
      free((char *) pv->u_string);
      pv->u_string = NULL;
    }
    break;
  case PT_INT3:
    ierr = H5LTget_attribute_int(group, ".", name, pv->u_int3); CE;
    break;
  case PT_FLOAT3:
    ierr = H5LTget_attribute_float(group, ".", name, pv->u_float3); CE;
    break;
  case PT_DOUBLE3:
    ierr = H5LTget_attribute_double(group, ".", name, pv->u_double3); CE;
    break;
  case PT_INT_ARRAY: {
    int attr = H5Aopen(group, name, H5P_DEFAULT); H5_CHK(attr);
    H5A_info_t ainfo;
    ierr = H5Aget_info(attr, &ainfo); CE;
    ierr = H5Aclose(attr); CE;
    pv->u_int_array.nr_vals = ainfo.data_size / sizeof(int);
    pv->u_int_array.vals = calloc(pv->u_int_array.nr_vals, sizeof(int));
    ierr = H5LTget_attribute_int(group, ".", name, pv->u_int_array.vals); CE;
    break;
  }
  case PT_PTR:
    break;
  default:
    mpi_printf(mrc_io_comm(io), "mrc_io_hdf5_parallel: not reading attr '%s' (type %d)\n", name, type);
    assert(0);
    break;
  }
  ierr = H5Gclose(group); CE;
}


#if 0 
static void
xdmf_write_m1(struct mrc_io *io, const char *path, struct mrc_fld *f1)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;

  if (mrc_fld_nr_patches(f1) != 1) {
    return;
  }

  hid_t group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT);
  const int i1 = 1, i0 = 0;
  H5LTset_attribute_int(group0, ".", "nr_patches", &i1, 1);

  struct xdmf_spatial *xs = xdmf_spatial_find(&file->xdmf_spatial_list,
					      mrc_domain_name(f1->_domain));
  if (0 && !xs) {
    xs = xdmf_spatial_create_f1(&file->xdmf_spatial_list,
				mrc_domain_name(f1->_domain), f1->_domain);
    xdmf_spatial_write_crds(xs, file, f1->_domain);
  }

  for (int m = 0; m < mrc_fld_nr_comps(f1); m++) {
    const char *fld_name = mrc_fld_comp_name(f1, m);
    if (!fld_name) {
      char tmp_fld_name[10];
      fld_name = tmp_fld_name;
      // FIXME: warn
      MHERE;
      sprintf(tmp_fld_name, "m%d", m);
    }
    //    xdmf_spatial_save_fld_info(xs, strdup(fld_name), strdup(path), false);

    hid_t group_fld = H5Gcreate(group0, fld_name, H5P_DEFAULT,
				H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_fld);
    char s_patch[10];
    sprintf(s_patch, "p%d", 0);

    const int *dims = mrc_fld_dims(f1);
    hsize_t mdims[1] = { dims[0] };
    hsize_t fdims[1] = { dims[0] };
    hsize_t off[1] = { 0 };
    
    hid_t group = H5Gcreate(group_fld, s_patch, H5P_DEFAULT,
			    H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
    H5LTset_attribute_int(group, ".", "global_patch", &i0, 1);
    hid_t filespace = H5Screate_simple(1, fdims, NULL);
    hid_t memspace = H5Screate_simple(1, mdims, NULL);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, off, NULL, fdims, NULL);
    
    hid_t dset = H5Dcreate(group, "1d", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
			   H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT,
	     &MRC_F1(f1, m, 0));
    H5Dclose(dset);
    H5Gclose(group);

    H5Gclose(group_fld);
  }
  H5Gclose(group0);
}

#endif


static void
hdf5_parallel_open(struct mrc_io *io, const char *mode)
{
  struct ph5 *ph5 = to_ph5(io);
  // FIXME: There may be some future speed improvement on lustre if we use
  // a real mpi_info here and set romio_[cb/ds]_write bits

  char filename[strlen(io->par.outdir) + strlen(io->par.basename) + 20];
  // using '.' as a separator causes all sorts of problems with io/analysis in
  // mrc-v3, but I don't want to break compatibility with the other io types
  sprintf(filename, "%s/%s%s%06d.h5", io->par.outdir, io->par.basename,
	  ph5->separator, io->step);

  hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist, io->obj.comm, MPI_INFO_NULL);
  
  if (strcmp(mode, "w") == 0) {
    ph5->h5_file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
  } else if (strcmp(mode, "r") == 0) {
    ph5->h5_file = H5Fopen(filename, H5F_ACC_RDONLY, plist);
  } else {
    assert(0);
  }

  H5Pclose(plist);

}

// ----------------------------------------------------------------------
// xdmf_parallel_close

static void
hdf5_parallel_close(struct mrc_io *io)
{
  struct ph5 *ph5 = to_ph5(io);

  H5Fclose(ph5->h5_file);
}


// ----------------------------------------------------------------------
// xdmf_parallel_write_fld

static void
hdf5_parallel_write_fld(struct mrc_io *io, const char *path, struct mrc_fld *fld)
{

  // FIXME: This should be able to work without a domain.
  // maybe sub mpi size in for nr_patches later on? 
  assert(fld->_domain);
  
  struct ph5 *ph5 = to_ph5(io);

  hid_t group0;
  if (H5Lexists(ph5->h5_file, path, H5P_DEFAULT) > 0) {
    group0 = H5Gopen(ph5->h5_file, path, H5P_DEFAULT);
  } else {
    group0 = H5Gcreate(ph5->h5_file, path, H5P_DEFAULT,
		       H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group0);
  }
  
 
  int nr_patches;
  mrc_domain_get_nr_global_patches(fld->_domain, &nr_patches);
  H5LTset_attribute_int(group0, ".", "nr_patches", &nr_patches, 1);

  // FIXME: This dumps everything together (best for aos). There's no single
  // write cached way to do separate the fields in either aos or soa, though
  // aos has worse data fragmentation.
  
  // line up patches in file
  int nr_spatial_dims = fld->_nr_spatial_dims;
  int nr_file_dims = nr_spatial_dims + 2;
  const int *fld_dims = mrc_fld_dims(fld);
  int nr_local_patches = mrc_fld_nr_patches(fld);

  hsize_t fdims[nr_file_dims],
    mdims[nr_file_dims]; // size of the full memory block;

  fdims[0] = nr_patches;
  mdims[0] = nr_local_patches;
  for (int d=1; d<nr_file_dims; d++) {
    fdims[d] = fld_dims[nr_file_dims - 1 - d];
    mdims[d] = fld->_ghost_dims[nr_file_dims - 1 - d];
  }

  // FIXME: crds write ghost points, but this is a really hacky way to do this
  // I actually would like it if we added an interface to choose to write ghosts
  // for mrc_flds, since then I could just read trafo bits back in too.
  if (nr_spatial_dims == 1) {
    fdims[2] += 2*fld->_nr_ghosts;
  }

  hid_t filespace = H5Screate_simple(nr_file_dims, fdims, NULL);
  hid_t memspace = H5Screate_simple(nr_file_dims, mdims, NULL);

  hid_t datatype;
  switch (mrc_fld_data_type(fld)) {
  case MRC_NT_FLOAT:
    datatype = H5T_NATIVE_FLOAT;
    break;
  case MRC_NT_DOUBLE:
    datatype = H5T_NATIVE_DOUBLE;
    break;
  case MRC_NT_INT:
    datatype = H5T_NATIVE_INT;
    break;
  default:
    assert(0);
  }

  hid_t dset = H5Dcreate(group0, "3d", datatype, filespace, H5P_DEFAULT,
			 H5P_DEFAULT, H5P_DEFAULT);
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  if (ph5->use_independent_io) {
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
  } else {
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
  }

  hsize_t mcount[nr_file_dims],  // number of elements we'll transfer in each slab
    moff[nr_file_dims], // offset of that slab into the full memory block
    foff[nr_file_dims]; // offset of the corresponding slab in the file
  struct mrc_patch_info info; 
  mrc_domain_get_local_patch_info(fld->_domain, 0, &info);
  mcount[0] = nr_local_patches;
  moff[0] = 0;
  foff[0] = info.global_patch;
  for (int d=1; d < nr_file_dims; d++) {
    mcount[d] = fdims[d];
    moff[d] = fld->_sw.vals[nr_file_dims - 1 - d];
    foff[d] = 0;
  }

  // FIXME: another hack for writing crds with ghosts
  if (nr_spatial_dims == 1) {
    moff[2] -= fld->_nr_ghosts;
  }

  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, foff, NULL, mcount, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, moff, NULL, mcount, NULL);

  
  H5Dwrite(dset, datatype, memspace, filespace, dxpl, fld->_nd->arr);
  
  H5Sclose(memspace);
  H5Dclose(dset);
  H5Sclose(filespace);
  H5Pclose(dxpl);
  
  
  H5Gclose(group0);
}

static void
hdf5_parallel_read_fld(struct mrc_io *io, const char *path, struct mrc_fld *fld)
{

  // FIXME: This should be able to work without a domain.
  // maybe sub mpi size in for nr_patches later on? 
  assert(fld->_domain);
  
  struct ph5 *ph5 = to_ph5(io);

  hid_t group0;
  group0 = H5Gopen(ph5->h5_file, path, H5P_DEFAULT);
  
 
  int nr_patches;
  mrc_domain_get_nr_global_patches(fld->_domain, &nr_patches);
  int file_patches;
  H5LTget_attribute_int(group0, ".", "nr_patches", &file_patches);
  assert(nr_patches == file_patches);

  
  // line up patches in file
  int nr_spatial_dims = fld->_nr_spatial_dims;
  int nr_file_dims = nr_spatial_dims + 2;
  const int *fld_dims = mrc_fld_dims(fld);
  int nr_local_patches = mrc_fld_nr_patches(fld);

  hsize_t fdims[nr_file_dims],
    mdims[nr_file_dims]; // size of the full memory block;

  fdims[0] = nr_patches;
  mdims[0] = nr_local_patches;
  for (int d=1; d<nr_file_dims; d++) {
    fdims[d] = fld_dims[nr_file_dims - 1 - d];
    mdims[d] = fld->_ghost_dims[nr_file_dims - 1 - d];
  }
  // FIXME: crds write ghost points, but this is a really hacky way to do this
  // I actually would like it if we added an interface to choose to write/read ghosts
  // for mrc_flds, since then I could just read trafo bits back in too.
  if (nr_spatial_dims == 1) {
    fdims[2] += 2*fld->_nr_ghosts;
  }


  hid_t memspace = H5Screate_simple(nr_file_dims, mdims, NULL);

  hid_t datatype;
  switch (mrc_fld_data_type(fld)) {
  case MRC_NT_FLOAT:
    datatype = H5T_NATIVE_FLOAT;
    break;
  case MRC_NT_DOUBLE:
    datatype = H5T_NATIVE_DOUBLE;
    break;
  case MRC_NT_INT:
    datatype = H5T_NATIVE_INT;
    break;
  default:
    assert(0);
  }

  hid_t dset = H5Dopen(group0, "3d", H5P_DEFAULT);
  hid_t filespace = H5Dget_space(dset);

  // Temporary checks:
  int file_rank = H5Sget_simple_extent_ndims(filespace);
  hsize_t read_dims[file_rank];
  H5Sget_simple_extent_dims(filespace, read_dims, NULL);
  assert(file_rank == nr_file_dims);
  for (int d = 0; d < file_rank; d++) {
    assert(read_dims[d] == fdims[d]);
  }
  // end checks

  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  if (ph5->use_independent_io) {
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
  } else {
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
  }

  hsize_t mcount[nr_file_dims],  // number of elements we'll transfer in each slab
    moff[nr_file_dims], // offset of that slab into the full memory block
    foff[nr_file_dims]; // offset of the corresponding slab in the file
  struct mrc_patch_info info; 
  mrc_domain_get_local_patch_info(fld->_domain, 0, &info);
  mcount[0] = nr_local_patches;
  moff[0] = 0;
  foff[0] = info.global_patch;
  for (int d=1; d < nr_file_dims; d++) {
    mcount[d] = fdims[d];
    moff[d] = fld->_sw.vals[nr_file_dims - 1 - d];
    foff[d] = 0;
  }

  // FIXME: another hack for writing crds with ghosts
  if (nr_spatial_dims == 1) {
    moff[2] -= fld->_nr_ghosts;
  }

  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, foff, NULL, mcount, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, moff, NULL, mcount, NULL);

  H5Dread(dset, datatype, memspace, filespace, dxpl, fld->_nd->arr);
  
  H5Sclose(memspace);
  H5Dclose(dset);
  H5Sclose(filespace);
  H5Pclose(dxpl);
  
  
  H5Gclose(group0);
}


// ----------------------------------------------------------------------
// mrc_io_ops_hdf5_parallel

struct mrc_io_ops mrc_io_hdf5_parallel_ops = {
  .name          = "hdf5_parallel",
  .size          = sizeof(struct ph5),
  .param_descr   = ph5_descr,
  .parallel      = true,
  .open          = hdf5_parallel_open,
  .close         = hdf5_parallel_close,
  .write_fld     = hdf5_parallel_write_fld,
  .read_fld      = hdf5_parallel_read_fld,
  .write_attr    = phdf5_write_attr,
  .read_attr     = phdf5_read_attr,
};

#endif
