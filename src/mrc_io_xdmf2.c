
#include <mrc_io_private.h>
#include <mrc_params.h>
#include <mrc_list.h>

#ifdef HAVE_HDF5_H

#include "mrc_io_xdmf_lib.h"

#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#define CE assert(ierr == 0)

// ----------------------------------------------------------------------

struct xdmf_file {
  hid_t h5_file;
  list_t xdmf_spatial_list;
};

struct xdmf {
  struct xdmf_file file;
  struct xdmf_temporal *xdmf_temporal;
};

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
  assert(strcmp(mode, "w") == 0);

  struct xdmf_file *file = &xdmf->file;
  char filename[strlen(io->par.outdir) + strlen(io->par.basename) + 20];
  sprintf(filename, "%s/%s.%06d_p%06d.h5", io->par.outdir, io->par.basename,
	  io->step, io->rank);
  file->h5_file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

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
    group = H5Gcreate(file->h5_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  switch (type) {
  case PT_SELECT:
  case PT_INT:
    H5LTset_attribute_int(group, ".", name, &pv->u_int, 1);
    break;
  case PT_BOOL: {
    int val = pv->u_bool;
    H5LTset_attribute_int(group, ".", name, &val, 1);
    break;
  }
  case PT_FLOAT:
    H5LTset_attribute_float(group, ".", name, &pv->u_float, 1);
    break;
  case PT_DOUBLE:
    H5LTset_attribute_double(group, ".", name, &pv->u_double, 1);
    break;
  case PT_STRING:
    H5LTset_attribute_string(group, ".", name, pv->u_string);
    break;
  case PT_INT3:
    H5LTset_attribute_int(group, ".", name, pv->u_int3, 3);
    break;
  case PT_FLOAT3:
    H5LTset_attribute_float(group, ".", name, pv->u_float3, 3);
    break;
  }
  H5Gclose(group);
}

static void
xdmf_spatial_write_mcrds(struct xdmf_spatial *xs, struct xdmf_file *file,
			 struct mrc_domain *domain)
{
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  if (xs->crds_done)
    return;

  xs->crds_done = true;

  for (int d = 0; d < 3; d++) {
    struct mrc_m1 *mcrd = crds->mcrd[d];

    hid_t group_crd1 = H5Gopen(file->h5_file, mrc_m1_name(mcrd), H5P_DEFAULT);

    mrc_m1_foreach_patch(mcrd, p) {
      struct mrc_m1_patch *mcrdp = mrc_m1_patch_get(mcrd, p);
      int im = mcrdp->im[0] + 2 * mcrdp->ib[0];
      // get node-centered coordinates
      float *crd_nc = calloc(im + 1, sizeof(*crd_nc));
      if (mcrdp->ib[0] < 0) {
	for (int i = 0; i <= im; i++) {
	  crd_nc[i] = .5 * (MRC_M1(mcrdp,0, i-1) + MRC_M1(mcrdp,0, i));
	}
      } else {
	for (int i = 1; i < im; i++) {
	  crd_nc[i] = .5 * (MRC_M1(mcrdp,0, i-1) + MRC_M1(mcrdp,0, i));
	}
	// extrapolate
	crd_nc[0]  = MRC_M1(mcrdp,0, 0) - .5 * (MRC_M1(mcrdp,0, 1) - MRC_M1(mcrdp,0, 0));
	crd_nc[im] = MRC_M1(mcrdp,0, im-1) + .5 * (MRC_M1(mcrdp,0, im-1) - MRC_M1(mcrdp,0, im-2));
      }
      hsize_t im1 = im + 1;
      char s_patch[10];
      sprintf(s_patch, "p%d", p);
      hid_t group_crdp = H5Gcreate(group_crd1, s_patch, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5LTmake_dataset_float(group_crdp, "1d", 1, &im1, crd_nc);
      H5Gclose(group_crdp);

      free(crd_nc);
      mrc_m1_patch_put(mcrd);
    }

    H5Gclose(group_crd1);
  }
}

static void
xdmf_write_m3(struct mrc_io *io, const char *path, struct mrc_m3 *m3)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  int ierr;

  hid_t group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT);
  H5LTset_attribute_int(group0, ".", "nr_patches", &m3->nr_patches, 1);

  struct xdmf_spatial *xs = xdmf_spatial_find(&file->xdmf_spatial_list,
					      mrc_domain_name(m3->domain));
  if (!xs) {
    xs = xdmf_spatial_create_m3(&file->xdmf_spatial_list,
				mrc_domain_name(m3->domain), m3->domain);
    xdmf_spatial_write_mcrds(xs, file, m3->domain);
  }

  for (int m = 0; m < m3->nr_comp; m++) {
    xdmf_spatial_save_fld_info(xs, strdup(m3->name[m]), strdup(path), false);

    hid_t group_fld = H5Gcreate(group0, m3->name[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    mrc_m3_foreach_patch(m3, p) {
      struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, p);

      char s_patch[10];
      sprintf(s_patch, "p%d", p);

      hsize_t mdims[3] = { m3p->im[2], m3p->im[1], m3p->im[0] };
      hsize_t fdims[3] = { m3p->im[2] + 2 * m3p->ib[2],
			   m3p->im[1] + 2 * m3p->ib[1],
			   m3p->im[0] + 2 * m3p->ib[0] };
      hsize_t off[3] = { -m3p->ib[2], -m3p->ib[1], -m3p->ib[0] };

      hid_t group = H5Gcreate(group_fld, s_patch, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(m3->domain, p, &info);
      H5LTset_attribute_int(group, ".", "global_patch", &info.global_patch, 1);
      hid_t filespace = H5Screate_simple(3, fdims, NULL);
      hid_t memspace = H5Screate_simple(3, mdims, NULL);
      H5Sselect_hyperslab(memspace, H5S_SELECT_SET, off, NULL, fdims, NULL);

      hid_t dset = H5Dcreate(group, "3d", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
			     H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT,
	       &MRC_M3(m3p, m, m3p->ib[0], m3p->ib[1], m3p->ib[2]));
      H5Dclose(dset);
      ierr = H5Gclose(group); CE;
      mrc_m3_patch_put(m3);
    }
    H5Gclose(group_fld);
  }

  H5Gclose(group0);
}


// ----------------------------------------------------------------------
// mrc_io_ops_xdmf

struct mrc_io_ops mrc_io_xdmf2_ops = {
  .name          = "xdmf2",
  .size          = sizeof(struct xdmf),
  .parallel      = true,
  .setup         = xdmf_setup,
  .destroy       = xdmf_destroy,
  .open          = xdmf_open,
  .close         = xdmf_close,
  .write_attr    = xdmf_write_attr,
  .write_m3      = xdmf_write_m3,
};


// ======================================================================


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
  H5Pset_fapl_mpio(plist, io->obj.comm, MPI_INFO_NULL);
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
// xdmf_spatial_write_mcrds_parallel

static void
xdmf_spatial_write_mcrds_parallel(struct xdmf_spatial *xs,
				  struct xdmf_file *file,
				  struct mrc_domain *domain)
{
  if (xs->crds_done)
    return;

  xs->crds_done = true;

  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);

  for (int d = 0; d < 3; d++) {
    struct mrc_m1 *mcrd = crds->mcrd[d];

    hid_t group_crd1 = H5Gcreate(file->h5_file, mrc_m1_name(mcrd), H5P_DEFAULT,
				 H5P_DEFAULT, H5P_DEFAULT);

    hid_t group_crdp = H5Gcreate(group_crd1, "p0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    hsize_t fdims[1] = { gdims[d] + 1 };
    hid_t filespace = H5Screate_simple(1, fdims, NULL);
    hid_t dset = H5Dcreate(group_crdp, "1d", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
			   H5P_DEFAULT, H5P_DEFAULT);

    mrc_m1_foreach_patch(mcrd, p) {
      struct mrc_m1_patch *mcrdp = mrc_m1_patch_get(mcrd, p);
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
      if (skip_write) {
	continue;
      }

      // get node-centered coordinates
      int im = info.ldims[d];
      float *crd_nc = calloc(im + 1, sizeof(*crd_nc));
      if (mcrdp->ib[0] < 0) {
	for (int i = 0; i <= im; i++) {
	  crd_nc[i] = .5 * (MRC_M1(mcrdp,0, i-1) + MRC_M1(mcrdp,0, i));
	}
      } else {
	for (int i = 1; i < im; i++) {
	  crd_nc[i] = .5 * (MRC_M1(mcrdp,0, i-1) + MRC_M1(mcrdp,0, i));
	}
	// extrapolate
	crd_nc[0]  = MRC_M1(mcrdp,0, 0) - .5 * (MRC_M1(mcrdp,0, 1) - MRC_M1(mcrdp,0, 0));
	crd_nc[im] = MRC_M1(mcrdp,0, im-1) + .5 * (MRC_M1(mcrdp,0, im-1) - MRC_M1(mcrdp,0, im-2));
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
// xdmf_parallel_write_m3

static void
xdmf_parallel_write_m3(struct mrc_io *io, const char *path, struct mrc_m3 *m3)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;

  hid_t group0;
  if (H5Lexists(file->h5_file, path, H5P_DEFAULT) > 0) {
    group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT);
  } else {
    group0 = H5Gcreate(file->h5_file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  int nr_1 = 1;
  H5LTset_attribute_int(group0, ".", "nr_patches", &nr_1, 1);

  struct xdmf_spatial *xs = xdmf_spatial_find(&file->xdmf_spatial_list,
					      mrc_domain_name(m3->domain));
  if (!xs) {
    xs = xdmf_spatial_create_m3_parallel(&file->xdmf_spatial_list,
					 mrc_domain_name(m3->domain),
					 m3->domain);
    xdmf_spatial_write_mcrds_parallel(xs, file, m3->domain);
  }

  int gdims[3];
  mrc_domain_get_global_dims(m3->domain, gdims);

  for (int m = 0; m < m3->nr_comp; m++) {
    xdmf_spatial_save_fld_info(xs, strdup(m3->name[m]), strdup(path), false);

    hid_t group_fld = H5Gcreate(group0, m3->name[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t group = H5Gcreate(group_fld, "p0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int i0 = 0;
    H5LTset_attribute_int(group, ".", "global_patch", &i0, 1);

    hsize_t fdims[3] = { gdims[0], gdims[1], gdims[2] };
    hid_t filespace = H5Screate_simple(3, fdims, NULL);
    hid_t dset = H5Dcreate(group, "3d", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
			   H5P_DEFAULT, H5P_DEFAULT);
    mrc_m3_foreach_patch(m3, p) {
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(m3->domain, p, &info);
      struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, p);

      hsize_t mdims[3] = { m3p->im[2], m3p->im[1], m3p->im[0] };
      hsize_t mcount[3] = { info.ldims[2], info.ldims[1], info.ldims[0] };
      hsize_t moff[3] = { -m3p->ib[2], -m3p->ib[1], -m3p->ib[0] };
      hsize_t foff[3] = { info.off[2], info.off[1], info.off[0] };

      hid_t filespace = H5Screate_simple(3, fdims, NULL);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, foff, NULL, mcount, NULL);
      hid_t memspace = H5Screate_simple(3, mdims, NULL);
      H5Sselect_hyperslab(memspace, H5S_SELECT_SET, moff, NULL, mcount, NULL);

      H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT,
	       &MRC_M3(m3p, m, m3p->ib[0], m3p->ib[1], m3p->ib[2]));

      mrc_m3_patch_put(m3);
    }
    H5Dclose(dset);
    H5Sclose(filespace);

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
