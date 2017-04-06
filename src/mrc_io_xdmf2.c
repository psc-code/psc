
#include <mrc_io_private.h>
#include <mrc_params.h>
#include <mrc_list.h>
#include "mrc_io_xdmf_lib.h"

#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------

struct xdmf_file {
  hid_t h5_file;
  list_t xdmf_spatial_list;
};

struct xdmf {
  struct xdmf_file file;
  struct xdmf_temporal *xdmf_temporal;
  int sw;
  // parallel only
  bool use_independent_io;
};

#define VAR(x) (void *)offsetof(struct xdmf, x)
static struct param xdmf2_descr[] = {
  { "sw"                     , VAR(sw)                      , PARAM_INT(0)           },
  {},
};
#undef VAR

#define VAR(x) (void *)offsetof(struct xdmf, x)
static struct param xdmf_parallel_descr[] = {
  { "use_independent_io"     , VAR(use_independent_io)      , PARAM_BOOL(false)      },
  {},
};
#undef VAR

#define to_xdmf(io) ((struct xdmf *)((io)->obj.subctx))

// ======================================================================
// xdmf

static void
xdmf_setup(struct mrc_io *io)
{
  struct xdmf *xdmf = to_xdmf(io);

  char filename[strlen(io->par.outdir) + strlen(io->par.basename) + 7];
  sprintf(filename, "%s/%s.xdmf", io->par.outdir, io->par.basename);
  xdmf->xdmf_temporal = xdmf_temporal_create(filename);

  mrc_io_setup_super(io);
}

static void
xdmf_destroy(struct mrc_io *io)
{
  struct xdmf *xdmf = to_xdmf(io);
  if (xdmf->xdmf_temporal) {
    xdmf_temporal_destroy(xdmf->xdmf_temporal);
  }
}

// ----------------------------------------------------------------------
// xdmf_open

static void
xdmf_open(struct mrc_io *io, const char *mode)
{
  struct xdmf *xdmf = to_xdmf(io);

  struct xdmf_file *file = &xdmf->file;
  char filename[strlen(io->par.outdir) + strlen(io->par.basename) + 20];
  sprintf(filename, "%s/%s.%06d_p%06d.h5", io->par.outdir, io->par.basename,
	  io->step, io->rank);


  if (strcmp(mode, "w") == 0) {
    file->h5_file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  } else if (strcmp(mode, "r") == 0) {
    file->h5_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  } else {
    assert(0);
  }

  xdmf_spatial_open(&file->xdmf_spatial_list);
}

static void
xdmf_close(struct mrc_io *io)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;

  H5Fclose(file->h5_file);
  xdmf_spatial_close(&file->xdmf_spatial_list, io, xdmf->xdmf_temporal);

  memset(file, 0, sizeof(*file));
}

static void
xdmf_write_attr(struct mrc_io *io, const char *path, int type,
		const char *name, union param_u *pv)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  
  hid_t group;
  if (H5Lexists(file->h5_file, path, H5P_DEFAULT) > 0) {
    group = H5Gopen(file->h5_file, path, H5P_DEFAULT);
  } else {
    group = H5Gcreate(file->h5_file, path, H5P_DEFAULT,
          H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
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
  case MRC_VAR_DOUBLE3:
    ierr = H5LTset_attribute_double(group, ".", name, pv->u_double3, 3); CE;
    break;
  case PT_INT_ARRAY: {
    hsize_t dims = pv->u_int_array.nr_vals;
    hid_t dataspace_id = H5Screate_simple(1, &dims, NULL); H5_CHK(dataspace_id);
    hid_t attr_id = H5Acreate(group, name, H5T_NATIVE_INT, dataspace_id,
			      H5P_DEFAULT, H5P_DEFAULT); H5_CHK(attr_id);
    if (dims > 0) {
      ierr = H5Awrite(attr_id, H5T_NATIVE_INT, pv->u_int_array.vals); CE;
    }
    ierr = H5Sclose(dataspace_id); CE;
    ierr = H5Aclose(attr_id); CE;
    break;
  }
  case PT_PTR:
    break;
  default:
    mpi_printf(mrc_io_comm(io), "mrc_io_xdmf2: not writing attr '%s' (type %d)\n",
      name, type);
    assert(0);
  }
  ierr = H5Gclose(group); CE;
}

static void
xdmf_read_attr(struct mrc_io *io, const char *path, int type,
    const char *name, union param_u *pv)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  
  hid_t group = H5Gopen(file->h5_file, path, H5P_DEFAULT); H5_CHK(group);

  int ierr;
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
    mpi_printf(mrc_io_comm(io), "mrc_io_xdmf2: not reading attr '%s' (type %d)\n", name, type);
    assert(0);
    break;
  }
  ierr = H5Gclose(group); CE;
}


// static void
// xdmf_spatial_write_mcrds_multi(struct mrc_io *io, struct xdmf_file *file,
// 			       struct mrc_domain *domain, int sw)
// {
//   struct mrc_crds *crds = mrc_domain_get_crds(domain);

//   for (int d = 0; d < 3; d++) {
//     struct mrc_fld *mcrd = crds->crd[d];

//     hid_t group_crd1 = H5Gopen(file->h5_file, mrc_io_obj_path(io, mcrd),
// 			       H5P_DEFAULT); H5_CHK(group_crd1);

//     mrc_m1_foreach_patch(mcrd, p) {
//       int im = mrc_fld_dims(mcrd)[0];
//       // get node-centered coordinates
//       float *crd_nc = calloc(im + 2*sw + 1, sizeof(*crd_nc));
//       if (mrc_fld_ghost_offs(mcrd)[0] < -sw) {
// 	for (int i = -sw; i <= im + sw; i++) {
// 	  crd_nc[i + sw] = .5 * (MRC_M1(mcrd,0, i-1, p) + MRC_M1(mcrd,0, i, p));
// 	}
//       } else {
// 	for (int i = 1-sw; i < im+sw; i++) {
// 	  crd_nc[i + sw] = .5 * (MRC_M1(mcrd,0, i-1, p) + MRC_M1(mcrd,0, i, p));
// 	}
// 	// extrapolate
// 	crd_nc[0      ] = MRC_M1(mcrd,0, -sw    , p) - .5 * (MRC_M1(mcrd,0, -sw+1  , p) - MRC_M1(mcrd,0, -sw    , p));
// 	crd_nc[im+2*sw] = MRC_M1(mcrd,0, im+sw-1, p) + .5 * (MRC_M1(mcrd,0, im+sw-1, p) - MRC_M1(mcrd,0, im+sw-2, p));
//       }
//       hsize_t im1 = im + 2*sw + 1;
//       char s_patch[10];
//       sprintf(s_patch, "p%d", p);
//       hid_t group_crdp = H5Gcreate(group_crd1, s_patch, H5P_DEFAULT,
// 				   H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_crdp);
//       H5LTmake_dataset_float(group_crdp, "1d", 1, &im1, crd_nc);
//       H5Gclose(group_crdp);

//       free(crd_nc);
//     }

//     H5Gclose(group_crd1);
//   }
// }

// static void
// xdmf_spatial_write_mcrds(struct mrc_io *io, struct xdmf_spatial *xs, struct xdmf_file *file,
// 			 struct mrc_domain *domain, int sw)
// {
//   if (xs->crds_done)
//     return;

//   xs->crds_done = true;

//   xdmf_spatial_write_mcrds_multi(io, file, domain, sw);
// }

// static void
// xdmf_spatial_write_crds_nonuni(struct xdmf_file *file,
// 			       struct mrc_domain *domain)
// {
//   struct mrc_crds *crds = mrc_domain_get_crds(domain);

//   for (int d = 0; d < 1; d++) {
//     struct mrc_fld *crd = crds->crd[d];

//     hid_t group_crd1 = H5Gcreate(file->h5_file, mrc_fld_name(crd), H5P_DEFAULT,
// 				 H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_crd1);

//     const int *dims = mrc_fld_dims(crd);
//     hsize_t hdims[1] = { dims[0] };
//     char s_patch[10];
//     sprintf(s_patch, "p%d", 0);
//     hid_t group_crdp = H5Gcreate(group_crd1, s_patch, H5P_DEFAULT,
// 				 H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_crdp);
//     H5LTmake_dataset_float(group_crdp, "1d", 1, hdims, &MRC_F1(crd,0, 0));
//     H5Gclose(group_crdp);

//     H5Gclose(group_crd1);
//   }
// }

// static void
// xdmf_spatial_write_crds(struct xdmf_spatial *xs, struct xdmf_file *file,
// 			struct mrc_domain *domain)
// {
//   if (xs->crds_done)
//     return;

//   xs->crds_done = true;

//   struct mrc_crds *crds = mrc_domain_get_crds(domain);
//   if (strcmp(mrc_crds_type(crds), "uniform") == 0) {
//     xdmf_spatial_write_crds_nonuni(file, domain); // FIXME
//   } else {
//     xdmf_spatial_write_crds_nonuni(file, domain);
//   }
// }

static void
xdmf_write_m3(struct mrc_io *io, const char *path, struct mrc_fld *m3_any)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  int ierr;

  // If we have an aos field, we need to get it as soa
  struct mrc_fld *m3 = m3_any;
  if (m3_any->_aos) {
    switch (mrc_fld_data_type(m3_any)) {
      case MRC_NT_FLOAT: m3 = mrc_fld_get_as(m3_any, "float"); break;
      case MRC_NT_DOUBLE: m3 = mrc_fld_get_as(m3_any, "double"); break;
      case MRC_NT_INT: m3 = mrc_fld_get_as(m3_any, "int"); break;
      default: assert(0);
    }
  }

  int gdims[3];
  mrc_domain_get_global_dims(m3->_domain, gdims);
  int sw[3];
  for (int d = 0; d < 3; d++) {
    sw[d] = gdims[d] > 1 ? xdmf->sw : 0;
  }
  hid_t group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT);
  int nr_patches = mrc_fld_nr_patches(m3);
  H5LTset_attribute_int(group0, ".", "nr_patches", &nr_patches, 1);

  struct xdmf_spatial *xs = xdmf_spatial_find(&file->xdmf_spatial_list,
                mrc_domain_name(m3->_domain));
  if (!xs) {
    xs = xdmf_spatial_create_m3(&file->xdmf_spatial_list,
        mrc_domain_name(m3->_domain), m3->_domain, io, xdmf->sw);
    //xdmf_spatial_write_mcrds(io, xs, file, m3->_domain, xdmf->sw);
  }

  for (int m = 0; m < mrc_fld_nr_comps(m3); m++) {
    char default_name[100];
    const char *compname;

    // If the comps aren't named just name them by their component number
    if ( !(compname = mrc_fld_comp_name(m3, m)) ) {
      sprintf(default_name, "_UNSET_%d", m);
      compname = (const char *) default_name;
    }

    xdmf_spatial_save_fld_info(xs, strdup(compname), strdup(path), false, mrc_fld_data_type(m3));

    hid_t group_fld = H5Gcreate(group0, compname, H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_fld);

    ierr = H5LTset_attribute_int(group_fld, ".", "m", &m, 1); CE;

    mrc_fld_foreach_patch(m3, p) {
      struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);

      char s_patch[10];
      sprintf(s_patch, "p%d", p);

      hsize_t mdims[3] = { m3->_ghost_dims[2], m3->_ghost_dims[1], m3->_ghost_dims[0] };
      hsize_t fdims[3] = { m3->_ghost_dims[2] + 2 * m3->_ghost_offs[2] + 2*sw[2],
         m3->_ghost_dims[1] + 2 * m3->_ghost_offs[1] + 2*sw[1],
         m3->_ghost_dims[0] + 2 * m3->_ghost_offs[0] + 2*sw[0]};
      hsize_t off[3] = { -m3->_ghost_offs[2]-sw[2], -m3->_ghost_offs[1]-sw[1], -m3->_ghost_offs[0]-sw[0] };
      /* mprintf("fdims %lld:%lld:%lld\n", fdims[0], fdims[1], fdims[2]); */
      /* mprintf("mdims %lld:%lld:%lld\n", mdims[0], mdims[1], mdims[2]); */
      /* mprintf("off   %lld:%lld:%lld\n", off[0], off[1], off[2]); */
      

      hid_t group = H5Gcreate(group_fld, s_patch, H5P_DEFAULT,
            H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(m3->_domain, p, &info);
      H5LTset_attribute_int(group, ".", "global_patch", &info.global_patch, 1);
      hid_t filespace = H5Screate_simple(3, fdims, NULL);
      hid_t memspace = H5Screate_simple(3, mdims, NULL);
      H5Sselect_hyperslab(memspace, H5S_SELECT_SET, off, NULL, fdims, NULL);

      hid_t dtype;
      void *arr;
      switch (mrc_fld_data_type(m3)) {
        case MRC_NT_FLOAT: {
          arr = (void *) &MRC_S5((m3p)->_fld, m3->_ghost_offs[0], m3->_ghost_offs[1], m3->_ghost_offs[2], m, (m3p)->_p);
          dtype = H5T_NATIVE_FLOAT; 
          break;
        }
        case MRC_NT_DOUBLE: {
          arr = (void *) &MRC_D5((m3p)->_fld, m3->_ghost_offs[0], m3->_ghost_offs[1], m3->_ghost_offs[2], m, (m3p)->_p);
          dtype = H5T_NATIVE_DOUBLE; 
          break;
        }
        case MRC_NT_INT: {
          arr = (void *) &MRC_I5((m3p)->_fld, m3->_ghost_offs[0], m3->_ghost_offs[1], m3->_ghost_offs[2], m, (m3p)->_p);
          dtype = H5T_NATIVE_INT; 
          break;
        }
        default: assert(0);
      }

      hid_t dset = H5Dcreate(group, "3d", dtype, filespace, H5P_DEFAULT,
           H5P_DEFAULT, H5P_DEFAULT);
      ierr = H5Dwrite(dset, dtype, memspace, filespace, H5P_DEFAULT, arr); CE;
      ierr = H5Dclose(dset); CE;
      ierr = H5Gclose(group); CE;
      mrc_fld_patch_put(m3);
    }
    H5Gclose(group_fld);
  }

  H5Gclose(group0);

  if (m3_any->_aos) {
    mrc_fld_put_as(m3, m3_any);
  }

}

struct read_m3_cb_data {
  struct mrc_io *io;
  struct mrc_fld *m3;
  int sw[3];
};

static herr_t
read_m3_cb(hid_t g_id, const char *name, const H5L_info_t *info, void *op_data)
{

  int ierr;
  struct read_m3_cb_data *data = op_data;
  struct mrc_fld *m3 = data->m3;
  int *sw = data->sw;

  hid_t group_fld = H5Gopen(g_id, name, H5P_DEFAULT); H5_CHK(group_fld);

  int m;
  ierr = H5LTget_attribute_int(group_fld, ".", "m", &m); CE;

  mrc_fld_foreach_patch(m3, p) {

    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);

    char s_patch[10];
    sprintf(s_patch, "p%d", p);

    hid_t group = H5Gopen(group_fld, s_patch, H5P_DEFAULT); H5_CHK(group);

    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(m3->_domain, p, &info);

    int file_gpatch;
    H5LTget_attribute_int(group, ".", "global_patch", &file_gpatch);
    assert(info.global_patch == file_gpatch);


    hsize_t mdims[3] = { m3->_ghost_dims[2], m3->_ghost_dims[1], m3->_ghost_dims[0] };
    hsize_t fdims[3] = { m3->_ghost_dims[2] + 2 * m3->_ghost_offs[2] + 2*sw[2],
       m3->_ghost_dims[1] + 2 * m3->_ghost_offs[1] + 2*sw[1],
       m3->_ghost_dims[0] + 2 * m3->_ghost_offs[0] + 2*sw[0]};
    hsize_t off[3] = { -m3->_ghost_offs[2]-sw[2], -m3->_ghost_offs[1]-sw[1], -m3->_ghost_offs[0]-sw[0] };

    hid_t memspace = H5Screate_simple(3, mdims, NULL);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, off, NULL, fdims, NULL);


    hid_t dtype;
    void *arr;
    switch (mrc_fld_data_type(m3)) {
      case MRC_NT_FLOAT: {
        arr = (void *) &MRC_S5((m3p)->_fld, m3->_ghost_offs[0], m3->_ghost_offs[1], m3->_ghost_offs[2], m, (m3p)->_p);
        dtype = H5T_NATIVE_FLOAT; 
        break;
      }
      case MRC_NT_DOUBLE: {
        arr = (void *) &MRC_D5((m3p)->_fld, m3->_ghost_offs[0], m3->_ghost_offs[1], m3->_ghost_offs[2], m, (m3p)->_p);
        dtype = H5T_NATIVE_DOUBLE; 
        break;
      }
      case MRC_NT_INT: {
        arr = (void *) &MRC_I5((m3p)->_fld, m3->_ghost_offs[0], m3->_ghost_offs[1], m3->_ghost_offs[2], m, (m3p)->_p);
        dtype = H5T_NATIVE_INT; 
        break;
      }
      default: assert(0);
    }

    hid_t dset = H5Dopen(group, "3d", H5P_DEFAULT);
    hid_t filespace = H5Dget_space(dset);
    ierr = H5Dread(dset, dtype, memspace, filespace, H5P_DEFAULT, arr); CE;
    ierr = H5Dclose(dset); CE;
    ierr = H5Gclose(group); CE;
    mrc_fld_patch_put(m3);
  }
  H5Gclose(group_fld);
  return 0;
}

static void
xdmf_read_m3(struct mrc_io *io, const char *path, struct mrc_fld *m3_any)
{

  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;

  // If we have an aos field, we need to get it as soa
  struct mrc_fld *m3 = m3_any;
  if (m3_any->_aos) {
    switch (mrc_fld_data_type(m3_any)) {
      case MRC_NT_FLOAT: m3 = mrc_fld_get_as(m3_any, "float"); break;
      case MRC_NT_DOUBLE: m3 = mrc_fld_get_as(m3_any, "double"); break;
      case MRC_NT_INT: m3 = mrc_fld_get_as(m3_any, "int"); break;
      default: assert(0);
    }
  }

  int gdims[3];
  mrc_domain_get_global_dims(m3->_domain, gdims);
  int sw[3];
  for (int d = 0; d < 3; d++) {
    sw[d] = gdims[d] > 1 ? xdmf->sw : 0;
  }
  hid_t group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT);
  int nr_patches = mrc_fld_nr_patches(m3);
  
  int file_patches;  
  H5LTget_attribute_int(group0, ".", "nr_patches", &file_patches);
  assert(nr_patches == file_patches);

  struct read_m3_cb_data cb_data = {
    .io        = io,
    .m3        = m3,
    .sw        = {sw[0], sw[1], sw[2]},
  };
  
  hsize_t idx = 0;
  H5Literate_by_name(group0, ".", H5_INDEX_NAME, H5_ITER_INC, &idx,
         read_m3_cb, &cb_data, H5P_DEFAULT);


  H5Gclose(group0);

  if (m3_any->_aos) {
    mrc_fld_put_as(m3, m3_any);
  }

}

static void
xdmf_write_m1(struct mrc_io *io, const char *path, struct mrc_fld *m1_any)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  int ierr;

  // If we have an aos field, we need to get it as soa
  struct mrc_fld *m1 = m1_any;
  if (m1_any->_aos) {
    switch (mrc_fld_data_type(m1_any)) {
      case MRC_NT_FLOAT: m1 = mrc_fld_get_as(m1_any, "float"); break;
      case MRC_NT_DOUBLE: m1 = mrc_fld_get_as(m1_any, "double"); break;
      case MRC_NT_INT: m1 = mrc_fld_get_as(m1_any, "int"); break;
      default: assert(0);
    }
  }

  int gdims[3];
  mrc_domain_get_global_dims(m1->_domain, gdims);
  int sw[1];
  for (int d = 0; d < 1; d++) {
    sw[d] = gdims[d] > 1 ? xdmf->sw : 0;
  }
  hid_t group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT);
  int nr_patches = mrc_fld_nr_patches(m1);
  H5LTset_attribute_int(group0, ".", "nr_patches", &nr_patches, 1);

  // struct xdmf_spatial *xs = xdmf_spatial_find(&file->xdmf_spatial_list,
  //               mrc_domain_name(m1->_domain));
  // if (!xs) {
  //   xs = xdmf_spatial_create_m1(&file->xdmf_spatial_list,
  //       mrc_domain_name(m1->_domain), m1->_domain, io, xdmf->sw);
  // }

  for (int m = 0; m < mrc_fld_nr_comps(m1); m++) {
    //xdmf_spatial_save_fld_info(xs, strdup(mrc_fld_comp_name(m1, m)), strdup(path), false);

    char default_name[100];
    const char *compname;

    // If the comps aren't named just name them by their component number
    if ( !(compname = mrc_fld_comp_name(m1, m)) ) {
      sprintf(default_name, "_UNSET_%d", m);
      compname = (const char *) default_name;
    }

    hid_t group_fld = H5Gcreate(group0, compname, H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_fld);

    ierr = H5LTset_attribute_int(group_fld, ".", "m", &m, 1); CE;

    mrc_fld_foreach_patch(m1, p) {
      struct mrc_fld_patch *m1p = mrc_fld_patch_get(m1, p);

      char s_patch[10];
      sprintf(s_patch, "p%d", p);

      hsize_t mdims[1] = { m1->_ghost_dims[0] };
      hsize_t fdims[1] = { m1->_ghost_dims[0] + 2 * m1->_ghost_offs[0] + 2*sw[0]};
      hsize_t off[1] = { -m1->_ghost_offs[0]-sw[0] };
      /* mprintf("fdims %lld:%lld:%lld\n", fdims[0], fdims[1], fdims[2]); */
      /* mprintf("mdims %lld:%lld:%lld\n", mdims[0], mdims[1], mdims[2]); */
      /* mprintf("off   %lld:%lld:%lld\n", off[0], off[1], off[2]); */
      

      hid_t group = H5Gcreate(group_fld, s_patch, H5P_DEFAULT,
            H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(m1->_domain, p, &info);
      H5LTset_attribute_int(group, ".", "global_patch", &info.global_patch, 1);
      hid_t filespace = H5Screate_simple(1, fdims, NULL);
      hid_t memspace = H5Screate_simple(1, mdims, NULL);
      H5Sselect_hyperslab(memspace, H5S_SELECT_SET, off, NULL, fdims, NULL);

      hid_t dtype;
      void *arr;
      switch (mrc_fld_data_type(m1)) {
        case MRC_NT_FLOAT: {
          arr = (void *) &MRC_S3((m1p)->_fld, m1->_ghost_offs[0], m, (m1p)->_p);
          dtype = H5T_NATIVE_FLOAT; 
          break;
        }
        case MRC_NT_DOUBLE: {
          arr = (void *) &MRC_D3((m1p)->_fld, m1->_ghost_offs[0], m, (m1p)->_p);
          dtype = H5T_NATIVE_DOUBLE; 
          break;
        }
        case MRC_NT_INT: {
          arr = (void *) &MRC_I3((m1p)->_fld, m1->_ghost_offs[0], m, (m1p)->_p);
          dtype = H5T_NATIVE_INT; 
          break;
        }
        default: assert(0);
      }

      hid_t dset = H5Dcreate(group, "1d", dtype, filespace, H5P_DEFAULT,
           H5P_DEFAULT, H5P_DEFAULT);
      ierr = H5Dwrite(dset, dtype, memspace, filespace, H5P_DEFAULT, arr); CE;
      ierr = H5Dclose(dset); CE;
      ierr = H5Gclose(group); CE;
      mrc_fld_patch_put(m1);
    }
    H5Gclose(group_fld);
  }

  H5Gclose(group0);
  if (m1_any->_aos) {
    mrc_fld_put_as(m1, m1_any);
  }

}

struct read_m1_cb_data {
  struct mrc_io *io;
  struct mrc_fld *m1;
  int sw[1];
};

static herr_t
read_m1_cb(hid_t g_id, const char *name, const H5L_info_t *info, void *op_data)
{
  int ierr;
  struct read_m1_cb_data *data = op_data;
  struct mrc_fld *m1 = data->m1;
  int *sw = data->sw;

  hid_t group_fld = H5Gopen(g_id, name, H5P_DEFAULT); H5_CHK(group_fld);

  int m;
  ierr = H5LTget_attribute_int(group_fld, ".", "m", &m); CE;

  mrc_fld_foreach_patch(m1, p) {
    struct mrc_fld_patch *m1p = mrc_fld_patch_get(m1, p);

    char s_patch[10];
    sprintf(s_patch, "p%d", p);

    hid_t group = H5Gopen(group_fld, s_patch, H5P_DEFAULT); H5_CHK(group);

    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(m1->_domain, p, &info);

    int file_gpatch;
    H5LTget_attribute_int(group, ".", "global_patch", &file_gpatch);
    assert(info.global_patch == file_gpatch);

    hsize_t mdims[1] = { m1->_ghost_dims[0] };
    hsize_t fdims[1] = { m1->_ghost_dims[0] + 2 * m1->_ghost_offs[0] + 2*sw[0]};
    hsize_t off[1] = { -m1->_ghost_offs[0]-sw[0] };


    hid_t memspace = H5Screate_simple(1, mdims, NULL);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, off, NULL, fdims, NULL);


    hid_t dtype;
    void *arr;
    switch (mrc_fld_data_type(m1)) {
      case MRC_NT_FLOAT: {
        arr = (void *) &MRC_S3((m1p)->_fld, m1->_ghost_offs[0], m, (m1p)->_p);
        dtype = H5T_NATIVE_FLOAT; 
        break;
      }
      case MRC_NT_DOUBLE: {
        arr = (void *) &MRC_D3((m1p)->_fld, m1->_ghost_offs[0], m, (m1p)->_p);
        dtype = H5T_NATIVE_DOUBLE; 
        break;
      }
      case MRC_NT_INT: {
        arr = (void *) &MRC_I3((m1p)->_fld, m1->_ghost_offs[0], m, (m1p)->_p);
        dtype = H5T_NATIVE_INT; 
        break;
      }
      default: assert(0);
    }

    hid_t dset = H5Dopen(group, "1d", H5P_DEFAULT);
    hid_t filespace = H5Dget_space(dset);
    ierr = H5Dread(dset, dtype, memspace, filespace, H5P_DEFAULT, arr); CE;
    ierr = H5Dclose(dset); CE;
    ierr = H5Gclose(group); CE;
    mrc_fld_patch_put(m1);
  }
  H5Gclose(group_fld);
  return 0;
}

static void
xdmf_read_m1(struct mrc_io *io, const char *path, struct mrc_fld *m1_any)
{

  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;

  // If we have an aos field, we need to get it as soa
  struct mrc_fld *m1 = m1_any;
  if (m1_any->_aos) {
    switch (mrc_fld_data_type(m1_any)) {
      case MRC_NT_FLOAT: m1 = mrc_fld_get_as(m1_any, "float"); break;
      case MRC_NT_DOUBLE: m1 = mrc_fld_get_as(m1_any, "double"); break;
      case MRC_NT_INT: m1 = mrc_fld_get_as(m1_any, "int"); break;
      default: assert(0);
    }
  }

  int gdims[3];
  mrc_domain_get_global_dims(m1->_domain, gdims);
  int sw[1];
  for (int d = 0; d < 1; d++) {
    sw[d] = gdims[d] > 1 ? xdmf->sw : 0;
  }
  hid_t group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT);
  int nr_patches = mrc_fld_nr_patches(m1);
  
  int file_patches;  
  H5LTget_attribute_int(group0, ".", "nr_patches", &file_patches);
  assert(nr_patches == file_patches);

  struct read_m1_cb_data cb_data = {
    .io        = io,
    .m1        = m1,
    .sw        = {sw[0]},
  };
  
  hsize_t idx = 0;
  H5Literate_by_name(group0, ".", H5_INDEX_NAME, H5_ITER_INC, &idx,
         read_m1_cb, &cb_data, H5P_DEFAULT);


  H5Gclose(group0);

  if (m1_any->_aos) {
    mrc_fld_put_as(m1, m1_any);
  }

}


// ----------------------------------------------------------------------
// mrc_io_ops_xdmf

struct mrc_io_ops mrc_io_xdmf2_ops = {
  .name          = "xdmf2",
  .size          = sizeof(struct xdmf),
  .param_descr   = xdmf2_descr,
  .parallel      = true,
  .setup         = xdmf_setup,
  .destroy       = xdmf_destroy,
  .open          = xdmf_open,
  .close         = xdmf_close,
  .write_attr    = xdmf_write_attr,
  .read_attr     = xdmf_read_attr,
  .write_m3      = xdmf_write_m3,
  .read_m3       = xdmf_read_m3,
  .write_m1      = xdmf_write_m1,
  .read_m1       = xdmf_read_m1,
};


// ======================================================================

#ifdef H5_HAVE_PARALLEL

// ----------------------------------------------------------------------
// xdmf_parallel_open

static void
xdmf_parallel_open(struct mrc_io *io, const char *mode)
{
  struct xdmf *xdmf = to_xdmf(io);
  assert(strcmp(mode, "w") == 0);

  struct xdmf_file *file = &xdmf->file;
  char filename[strlen(io->par.outdir) + strlen(io->par.basename) + 20];
  sprintf(filename, "%s/%s.%06d_p%06d.h5", io->par.outdir, io->par.basename,
	  io->step, 0);

  hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
#ifdef H5_HAVE_PARALLEL
  H5Pset_fapl_mpio(plist, io->obj.comm, MPI_INFO_NULL);
#endif
  file->h5_file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
  H5Pclose(plist);

  xdmf_spatial_open(&file->xdmf_spatial_list);
}

// ----------------------------------------------------------------------
// xdmf_parallel_close

static void
xdmf_parallel_close(struct mrc_io *io)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;

  H5Fclose(file->h5_file);
  xdmf_spatial_close(&file->xdmf_spatial_list, io, xdmf->xdmf_temporal);

  memset(file, 0, sizeof(*file));
}

// ----------------------------------------------------------------------
// xdmf_spatial_write_mcrds_multi_parallel

static void
xdmf_spatial_write_mcrds_multi_parallel(struct xdmf_file *file,
					struct mrc_domain *domain)
{
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);

  for (int d = 0; d < 3; d++) {
    struct mrc_fld *mcrd = crds->crd[d];

    hid_t group_crd1 = H5Gcreate(file->h5_file, mrc_fld_name(mcrd), H5P_DEFAULT,
				 H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_crd1);

    hid_t group_crdp = H5Gcreate(group_crd1, "p0", H5P_DEFAULT,
				 H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_crdp);
    
    hsize_t fdims[1] = { gdims[d] + 1 };
    hid_t filespace = H5Screate_simple(1, fdims, NULL);
    hid_t dset = H5Dcreate(group_crdp, "1d", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
			   H5P_DEFAULT, H5P_DEFAULT);

    mrc_m1_foreach_patch(mcrd, p) {
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(domain, p, &info);
      bool skip_write = false;
      for (int dd = 0; dd < 3; dd++) {
	if (d != dd && info.off[dd] != 0) {
	  skip_write = true;
	  break;
	}
      }
      // indep I/O only!
      // FIXME, do collective ?
      if (skip_write) {
	continue;
      }

      // get node-centered coordinates
      int im = info.ldims[d];
      float *crd_nc = calloc(im + 1, sizeof(*crd_nc));
      if (mrc_fld_ghost_offs(mcrd)[0] < 0) {
	for (int i = 0; i <= im; i++) {
	  crd_nc[i] = .5 * (MRC_M1(mcrd,0, i-1, p) + MRC_M1(mcrd,0, i, p));
	}
      } else {
	for (int i = 1; i < im; i++) {
	  crd_nc[i] = .5 * (MRC_M1(mcrd,0, i-1, p) + MRC_M1(mcrd,0, i, p));
	}
	// extrapolate
	crd_nc[0]  = MRC_M1(mcrd,0, 0   , p) - .5 * (MRC_M1(mcrd,0, 1   , p) - MRC_M1(mcrd,0, 0   , p));
	crd_nc[im] = MRC_M1(mcrd,0, im-1, p) + .5 * (MRC_M1(mcrd,0, im-1, p) - MRC_M1(mcrd,0, im-2, p));
      }

      hsize_t mdims[1] = { info.ldims[d] + (info.off[d] == 0 ? 1 : 0) };
      hsize_t off[1] = { info.off[d] + (info.off[d] == 0 ? 0 : 1) };
      hid_t memspace = H5Screate_simple(1, mdims, NULL);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, off, NULL, mdims, NULL);

      H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT,
	       &crd_nc[(info.off[d] == 0) ? 0 : 1]);

      H5Sclose(memspace);
    }
    
    H5Dclose(dset);
    H5Sclose(filespace);

    H5Gclose(group_crdp);
    H5Gclose(group_crd1);
  }
}

// ----------------------------------------------------------------------
// xdmf_spatial_write_crds_multi_parallel

static void
xdmf_spatial_write_crds_multi_parallel(struct mrc_io *io, struct mrc_domain *domain)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  int gdims[3], np[3], nr_global_patches, nr_patches;
  mrc_domain_get_global_dims(domain, gdims);
  mrc_domain_get_nr_procs(domain, np);
  mrc_domain_get_nr_global_patches(domain, &nr_global_patches);
  mrc_domain_get_patches(domain, &nr_patches);

  for (int d = 0; d < 3; d++) {
    struct mrc_fld *crd = crds->crd[d];

    MPI_Request *send_reqs = calloc(nr_patches, sizeof(*send_reqs));
    int nr_send_reqs = 0;
    assert(nr_patches == 1); // otherwise need to redo tmp_nc
    float *tmp_nc = NULL;
    for (int p = 0; p < nr_patches; p++) {
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(domain, p, &info);
      bool skip_write = false;
      for (int dd = 0; dd < 3; dd++) {
	if (d != dd && info.off[dd] != 0) {
	  skip_write = true;
	  break;
	}
      }
      if (!skip_write) {
	tmp_nc = calloc(info.ldims[d] + 1, sizeof(*tmp_nc));

	// get node-centered coordinates
	if (crd->_sw.vals[0] > 0) {
	  for (int i = 0; i <= info.ldims[d]; i++) {
	    tmp_nc[i] = .5 * (MRC_F1(crd,0, i-1) + MRC_F1(crd,0, i));
	  }
	} else {
	  int ld = info.ldims[d];
	  for (int i = 1; i < ld; i++) {
	    tmp_nc[i] = .5 * (MRC_F1(crd,0, i-1) + MRC_F1(crd,0, i));
	  }
	  // extrapolate
	  tmp_nc[0]  = MRC_F1(crd,0, 0) - .5 * (MRC_F1(crd,0, 1) - MRC_F1(crd,0, 0));
	  tmp_nc[ld] = MRC_F1(crd,0, ld-1) + .5 * (MRC_F1(crd,0, ld-1) - MRC_F1(crd,0, ld-2));
	}

 	/* mprintf("Isend off %d %d %d gp %d\n", info.off[0], info.off[1], info.off[2], */
	/* 	info.global_patch); */
	MPI_Isend(tmp_nc + (info.off[d] == 0 ? 0 : 1),
		  info.ldims[d] + (info.off[d] == 0 ? 1 : 0), MPI_FLOAT,
		  0, info.global_patch,
		  mrc_io_comm(io), &send_reqs[nr_send_reqs++]);
      }
    }

    int im = gdims[d];
    float *crd_nc = NULL;

    if (io->rank == 0) { // only on first writer
      crd_nc = calloc(im + 1, sizeof(*crd_nc));
      MPI_Request *recv_reqs = calloc(np[d], sizeof(*recv_reqs));
      int nr_recv_reqs = 0;
      for (int gp = 0; gp < nr_global_patches; gp++) {
	struct mrc_patch_info info;
	mrc_domain_get_global_patch_info(domain, gp, &info);
	bool skip_write = false;
	for (int dd = 0; dd < 3; dd++) {
	  if (d != dd && info.off[dd] != 0) {
	    skip_write = true;
	    break;
	  }
	}
	if (skip_write) {
	  continue;
	}
	/* mprintf("Irecv off %d %d %d gp %d\n", info.off[0], info.off[1], info.off[2], gp); */
	MPI_Irecv(&crd_nc[info.off[d] + (info.off[d] == 0 ? 0 : 1)],
		  info.ldims[d] + (info.off[d] == 0 ? 1 : 0), MPI_FLOAT, info.rank,
		  gp, mrc_io_comm(io), &recv_reqs[nr_recv_reqs++]);
      }
      assert(nr_recv_reqs == np[d]);

      MPI_Waitall(nr_recv_reqs, recv_reqs, MPI_STATUSES_IGNORE);
      free(recv_reqs);
    }

    MPI_Waitall(nr_send_reqs, send_reqs, MPI_STATUSES_IGNORE);
    free(tmp_nc);
    free(send_reqs);

    // this group has been generated by generic mrc_obj code
    hid_t group_crd1 = H5Gopen(file->h5_file, mrc_fld_name(crd),
			       H5P_DEFAULT); H5_CHK(group_crd1);

    hid_t group_crdp = H5Gcreate(group_crd1, "p0", H5P_DEFAULT,
				 H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_crd1);
    
    hsize_t fdims[1] = { gdims[d] + 1 };
    hid_t filespace = H5Screate_simple(1, fdims, NULL);
    hid_t dset = H5Dcreate(group_crdp, "1d", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
			   H5P_DEFAULT, H5P_DEFAULT);

    hid_t memspace;
    if (io->rank == 0) {
      memspace = H5Screate_simple(1, fdims, NULL);
    } else {
      memspace = H5Screate(H5S_NULL);
      H5Sselect_none(memspace);
      H5Sselect_none(filespace);
    }
    
    H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, crd_nc);
    
    H5Sclose(memspace);
    
    H5Dclose(dset);
    H5Sclose(filespace);

    H5Gclose(group_crdp);
    H5Gclose(group_crd1);

    free(crd_nc);
  }
}

// ----------------------------------------------------------------------
// xdmf_spatial_write_crds_uniform_parallel

static void
xdmf_spatial_write_crds_uniform_parallel(struct xdmf_file *file,
					 struct mrc_domain *domain)
{
  struct mrc_crds *crds = mrc_domain_get_crds(domain);

  double xl[3], xh[3];
  mrc_crds_get_param_double3(crds, "l", xl);
  mrc_crds_get_param_double3(crds, "h", xh);

  for (int d = 0; d < 3; d++) {
    struct mrc_fld *mcrd = crds->crd[d];

    hid_t group_crd1 = H5Gcreate(file->h5_file, mrc_fld_name(mcrd), H5P_DEFAULT,
				 H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_crd1);

    hid_t group_crdp = H5Gcreate(group_crd1, "p0", H5P_DEFAULT,
				 H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_crd1);
    H5LTset_attribute_double(group_crdp, ".", "xl", &xl[d], 1);
    H5LTset_attribute_double(group_crdp, ".", "xh", &xh[d], 1);
    H5Gclose(group_crdp);
    H5Gclose(group_crd1);
  }
}

// ----------------------------------------------------------------------
// xdmf_spatial_write_mcrds_parallel

static void
xdmf_spatial_write_mcrds_parallel(struct xdmf_spatial *xs,
				  struct mrc_io *io,
				  struct mrc_domain *domain)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  if (xs->crds_done)
    return;

  xs->crds_done = true;

  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  if (strcmp(mrc_crds_type(crds), "uniform") == 0) {
    xdmf_spatial_write_crds_uniform_parallel(file, domain);
  } else if (strcmp(mrc_crds_type(crds), "multi") == 0) {
    // FIXME, broken since not all are writers
    assert(0);
    xdmf_spatial_write_mcrds_multi_parallel(file, domain);
  } else if (strcmp(mrc_crds_type(crds), "rectilinear") == 0) {
    // FIXME, should do XDMF uniform, or rather just use m1, not f1
    xdmf_spatial_write_crds_multi_parallel(io, domain);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// xdmf_parallel_write_m3

static void
xdmf_parallel_write_m3(struct mrc_io *io, const char *path, struct mrc_fld *m3)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;

  hid_t group0;
  if (H5Lexists(file->h5_file, path, H5P_DEFAULT) > 0) {
    group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT);
    MHERE;
  } else {
    group0 = H5Gcreate(file->h5_file, path, H5P_DEFAULT,
		       H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group0);
  }
  int nr_1 = 1;
  H5LTset_attribute_int(group0, ".", "nr_patches", &nr_1, 1);

  struct xdmf_spatial *xs = xdmf_spatial_find(&file->xdmf_spatial_list,
					      mrc_domain_name(m3->_domain));
  int gdims[3];
  mrc_domain_get_global_dims(m3->_domain, gdims);

  if (!xs) {
    int off[3] = {};
    xs = xdmf_spatial_create_m3_parallel(&file->xdmf_spatial_list,
					 mrc_domain_name(m3->_domain),
					 m3->_domain, off, gdims, io);
    xdmf_spatial_write_mcrds_parallel(xs, io, m3->_domain);
  }

  int nr_patches;
  mrc_domain_get_patches(m3->_domain, &nr_patches);
  int nr_patches_max;
  // FIXME, mrc_domain may know / cache
  MPI_Allreduce(&nr_patches, &nr_patches_max, 1, MPI_INT, MPI_MAX,
		mrc_domain_comm(m3->_domain));

  for (int m = 0; m < mrc_fld_nr_comps(m3); m++) {
    xdmf_spatial_save_fld_info(xs, strdup(mrc_fld_comp_name(m3, m)), strdup(path), false, mrc_fld_data_type(m3));

    hid_t group_fld = H5Gcreate(group0, mrc_fld_comp_name(m3, m), H5P_DEFAULT,
				H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_fld);
    hid_t group = H5Gcreate(group_fld, "p0", H5P_DEFAULT,
			    H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
    int i0 = 0;
    H5LTset_attribute_int(group, ".", "global_patch", &i0, 1);

    hsize_t fdims[3] = { gdims[2], gdims[1], gdims[0] };
    hid_t filespace = H5Screate_simple(3, fdims, NULL);
    hid_t dset = H5Dcreate(group, "3d", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
			   H5P_DEFAULT, H5P_DEFAULT);
    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
#ifdef H5_HAVE_PARALLEL
    if (xdmf->use_independent_io) {
      H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
    } else {
      H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
    }
#endif
    for (int p = 0; p < nr_patches_max; p++) {
      if (p >= nr_patches) {
	if (xdmf->use_independent_io)
	  continue;

	// for collective I/O write nothing if no patch left,
	// but still call H5Dwrite()
	H5Sselect_none(filespace);
	hid_t memspace = H5Screate(H5S_NULL);
	H5Sselect_none(memspace);
	H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, dxpl, NULL);
	H5Sclose(memspace);
	continue;
      }
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(m3->_domain, p, &info);
      struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);

      hsize_t mdims[3] = { m3->_ghost_dims[2], m3->_ghost_dims[1], m3->_ghost_dims[0] };
      hsize_t mcount[3] = { info.ldims[2], info.ldims[1], info.ldims[0] };
      hsize_t moff[3] = { -m3->_ghost_offs[2], -m3->_ghost_offs[1], -m3->_ghost_offs[0] };
      hsize_t foff[3] = { info.off[2], info.off[1], info.off[0] };

      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, foff, NULL, mcount, NULL);
      hid_t memspace = H5Screate_simple(3, mdims, NULL);
      H5Sselect_hyperslab(memspace, H5S_SELECT_SET, moff, NULL, mcount, NULL);

      H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, dxpl,
	       &MRC_M3(m3p, m, m3->_ghost_offs[0], m3->_ghost_offs[1], m3->_ghost_offs[2]));

      H5Sclose(memspace);

      mrc_fld_patch_put(m3);
    }
    H5Dclose(dset);
    H5Sclose(filespace);
    H5Pclose(dxpl);

    H5Gclose(group);
    H5Gclose(group_fld);
  }
  H5Gclose(group0);
}


// ----------------------------------------------------------------------
// mrc_io_ops_xdmf_parallel

struct mrc_io_ops mrc_io_xdmf2_parallel_ops = {
  .name          = "xdmf2_parallel",
  .size          = sizeof(struct xdmf),
  .param_descr   = xdmf_parallel_descr,
  .parallel      = true,
  .setup         = xdmf_setup,
  .destroy       = xdmf_destroy,
  .open          = xdmf_parallel_open,
  .close         = xdmf_parallel_close,
#if 0
  .write_attr    = xdmf_write_attr,
#endif
  .write_m3      = xdmf_parallel_write_m3,
};

#endif
