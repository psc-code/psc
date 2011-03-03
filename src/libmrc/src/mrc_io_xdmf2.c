
#include <mrc_io_private.h>
#include <mrc_params.h>
#include <mrc_list.h>

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

#define MAX_XDMF_FLD_INFO (30)

struct xdmf_fld_info {
  char *name;
  char *path;
  bool is_vec;
};

struct xdmf_spatial {
  char *name; //< from domain::name

  bool crds_done;

  int nr_global_patches;
  struct mrc_patch_info *patch_infos;

  int nr_fld_info;
  struct xdmf_fld_info fld_info[MAX_XDMF_FLD_INFO];

  list_t entry; //< on xdmf_file::xdmf_spatial_list
};

struct xdmf_temporal_step {
  list_t entry;
  char filename[0];
};

struct xdmf_temporal {
  char *filename;
  list_t timesteps;
};

#define to_xdmf(io) ((struct xdmf *)((io)->obj.subctx))

// ======================================================================
// xdmf_write_header

static void
xdmf_write_header(FILE *f)
{
  fprintf(f, "<?xml version='1.0' ?>\n");
  fprintf(f, "<Xdmf xmlns:xi='http://www.w3.org/2001/XInclude' Version='2.0'>\n");
}

static void
xdmf_write_topology_m3(FILE *f, int im[3], const char *filename, int p)
{
  // FIXME crd[012] hardcoded, should use mrc_m1_name()
  fprintf(f, "     <Topology TopologyType=\"3DRectMesh\" Dimensions=\"%d %d %d\"/>\n",
	  im[2] + 1, im[1] + 1, im[0] + 1);
  fprintf(f, "     <Geometry GeometryType=\"VXVYVZ\">\n");
  fprintf(f, "     <DataItem Name=\"VX\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[0] + 1);
  fprintf(f, "        %s:/crd0/p%d/1d\n", filename, p);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VY\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[1] + 1);
  fprintf(f, "        %s:/crd1/p%d/1d\n", filename, p);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VZ\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[2] + 1);
  fprintf(f, "        %s:/crd2/p%d/1d\n", filename, p);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     </Geometry>\n");
  fprintf(f, "\n");
}

static void
xdmf_write_fld_m3(FILE *f, struct xdmf_fld_info *fld_info, int im[3],
		  const char *filename, int p)
{
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Cell\">\n",
	  fld_info->name, fld_info->is_vec ? "Vector" : "Scalar");
  fprintf(f, "       <DataItem Dimensions=\"%d %d %d%s\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
	  im[2], im[1], im[0], fld_info->is_vec ? " 3" : "");
  fprintf(f, "        %s:/%s/%s/p%d/3d\n", filename, fld_info->path, fld_info->name, p);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
}

// ======================================================================
// xdmf_temporal

static struct xdmf_temporal *
xdmf_temporal_create(const char *filename)
{
  struct xdmf_temporal *xt = calloc(1, sizeof(*xt));
  INIT_LIST_HEAD(&xt->timesteps);
  xt->filename = strdup(filename);
  return xt;
}

static void
xdmf_temporal_destroy(struct xdmf_temporal *xt)
{
  while (!list_empty(&xt->timesteps)) {
    struct xdmf_temporal_step *xt_step =
      list_entry(xt->timesteps.next, struct xdmf_temporal_step, entry);
    list_del(&xt_step->entry);
    free(xt_step);
  }
  free(xt->filename);
  free(xt);
}

static void
xdmf_temporal_append(struct xdmf_temporal *xt, const char *fname_spatial)
{
  struct xdmf_temporal_step *xt_step =
    malloc(sizeof(*xt_step) + strlen(fname_spatial) + 1);
  strcpy(xt_step->filename, fname_spatial);
  list_add_tail(&xt_step->entry, &xt->timesteps);
}

static void
xdmf_temporal_write(struct xdmf_temporal *xt)
{
  // It'd be easier to create those files line by line as time goes by.
  // However, then we won't be able to get the timeseries into Paraview
  // until the solution is all finished.
  // So this version rewrites the xdmf file completely every timestep,
  // which needs however some way to figure out what times we've written
  // before.

  FILE *f = fopen(xt->filename, "w");

  xdmf_write_header(f);
  fprintf(f, "<Domain>\n");
  fprintf(f, "  <Grid GridType='Collection' CollectionType='Temporal'>\n");
  struct xdmf_temporal_step *xt_step;
  list_for_each_entry(xt_step, &xt->timesteps, entry) {
    fprintf(f, "  <xi:include href='%s' xpointer='xpointer(//Xdmf/Domain/Grid)'/>\n",
	    xt_step->filename);
  }
  fprintf(f, "  </Grid>\n");
  fprintf(f, "  </Domain>\n");
  fprintf(f, "</Xdmf>\n");
  fclose(f);
}

// ======================================================================
// xdmf_spatial -- corresponds to one mrc_domain

static void
xdmf_spatial_open(struct xdmf_file *file)
{
  INIT_LIST_HEAD(&file->xdmf_spatial_list);
}

static struct xdmf_spatial *
xdmf_spatial_find(struct xdmf_file *file, const char *name)
{
  struct xdmf_spatial *xs;
  list_for_each_entry(xs, &file->xdmf_spatial_list, entry) {
    if (strcmp(xs->name, name) == 0) {
      return xs;
    }
  }
  return NULL;
}

static struct xdmf_spatial *
xdmf_spatial_create_m3d(struct xdmf_file *file, const char *name, 
			struct mrc_domain *domain)
{
  // OPT, we could skip this on procs which aren't writing xdmf

  struct xdmf_spatial *xs = calloc(1, sizeof(*xs));
  xs->name = strdup(name);
  mrc_domain_get_nr_global_patches(domain, &xs->nr_global_patches);
  xs->patch_infos = calloc(xs->nr_global_patches, sizeof(*xs->patch_infos));

  for (int gp = 0; gp < xs->nr_global_patches; gp++) {
    mrc_domain_get_global_patch_info(domain, gp, &xs->patch_infos[gp]);
  }

  list_add_tail(&xs->entry, &file->xdmf_spatial_list);
  return xs;
}

static void
xdmf_spatial_save_fld_info(struct xdmf_spatial *xs, char *fld_name,
			   char *path, bool is_vec)
{
  assert(xs->nr_fld_info < MAX_XDMF_FLD_INFO);
  struct xdmf_fld_info *fld_info = &xs->fld_info[xs->nr_fld_info++];

  fld_info->name = fld_name;
  fld_info->path = path;
  fld_info->is_vec = is_vec;
}


static void
xdmf_spatial_write(struct xdmf_spatial *xs, const char *filename,
		   struct mrc_io *io)
{
  FILE *f = fopen(filename, "w");
  xdmf_write_header(f);
  fprintf(f, "<Domain>\n");
  fprintf(f, "<Grid GridType=\"Collection\" CollectionType=\"Spatial\">\n");
  fprintf(f, "   <Time Type=\"Single\" Value=\"%g\" />\n", io->time);
  for (int s = 0; s < xs->nr_global_patches; s++) {
    fprintf(f, "   <Grid Name=\"patch-%s-%d\" GridType=\"Uniform\">\n", xs->name, s);
    
    char fname[strlen(io->par.outdir) + strlen(io->par.basename) + 20];
    int rank = xs->patch_infos[s].rank;
    int patch = xs->patch_infos[s].patch;
    sprintf(fname, "%s.%06d_p%06d.h5", io->par.basename, io->step, rank);
    int *ldims = xs->patch_infos[s].ldims;
    xdmf_write_topology_m3(f, ldims, fname, patch);
    
    for (int m = 0; m < xs->nr_fld_info; m++) {
      xdmf_write_fld_m3(f, &xs->fld_info[m], ldims, fname, patch);
    }
    fprintf(f, "   </Grid>\n");
  }
  fprintf(f, "</Grid>\n");
  fprintf(f, "</Domain>\n");
  fprintf(f, "</Xdmf>\n");
  fclose(f);
}

static void
xdmf_spatial_destroy(struct xdmf_spatial *xs)
{
  list_del(&xs->entry);
  free(xs->patch_infos);
  free(xs);
}

static void
xdmf_spatial_close(struct xdmf_file *file, struct mrc_io *io,
		   struct xdmf_temporal *xt)
{
  while (!list_empty(&file->xdmf_spatial_list)) {
    struct xdmf_spatial *xs = list_entry(file->xdmf_spatial_list.next, typeof(*xs), entry);
    
    char fname_spatial[strlen(io->par.outdir) + strlen(io->par.basename) + 20];
    sprintf(fname_spatial, "%s/%s.%06d.xdmf", io->par.outdir, io->par.basename, io->step);
    xdmf_spatial_write(xs, fname_spatial, io);
    xdmf_spatial_destroy(xs);

    if (io->rank == 0) {
      xdmf_temporal_append(xt, fname_spatial);
      xdmf_temporal_write(xt);
    }
  }
}

// ======================================================================
// xdmf

static void
xdmf_setup(struct mrc_obj *obj)
{
  struct mrc_io *io = to_mrc_io(obj);
  struct xdmf *xdmf = to_xdmf(io);

  char filename[strlen(io->par.outdir) + strlen(io->par.basename) + 7];
  sprintf(filename, "%s/%s.xdmf", io->par.outdir, io->par.basename);
  xdmf->xdmf_temporal = xdmf_temporal_create(filename);
}

static void
xdmf_destroy(struct mrc_obj *obj)
{
  struct mrc_io *io = to_mrc_io(obj);
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

  xdmf_spatial_open(file);
}

static void
xdmf_close(struct mrc_io *io)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;

  H5Fclose(file->h5_file);
  xdmf_spatial_close(file, io, xdmf->xdmf_temporal);

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
			 struct mrc_domain *domain, int sw)
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
      int im = mcrdp->im[0];
      // get node-centered coordinates
      float *crd_nc = calloc(im + 1, sizeof(*crd_nc));
      if (sw > 0) {
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
      char name[20];
      sprintf(name, "p%d", p);
      hid_t group_crdp = H5Gcreate(group_crd1, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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

  struct xdmf_spatial *xs = xdmf_spatial_find(file, mrc_domain_name(m3->domain));
  if (!xs) {
    xs = xdmf_spatial_create_m3d(file, mrc_domain_name(m3->domain), m3->domain);
    xdmf_spatial_write_mcrds(xs, file, m3->domain, m3->sw);
  }

  for (int m = 0; m < m3->nr_comp; m++) {
    hid_t group_fld = H5Gcreate(group0, m3->name[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    mrc_m3_foreach_patch(m3, p) {
      struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, p);

      char s_patch[10];
      sprintf(s_patch, "p%d", p);

      xdmf_spatial_save_fld_info(xs, strdup(m3->name[m]), strdup(path), false);
      hsize_t hdims[3] = { m3p->im[2], m3p->im[1], m3p->im[0] };
      hid_t group = H5Gcreate(group_fld, s_patch, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      ierr = H5LTmake_dataset_float(group, "3d", 3, hdims, 
				    &MRC_M3(m3p, m, m3p->ib[0], m3p->ib[1], m3p->ib[2])); CE;
      ierr = H5Gclose(group); CE;
      mrc_m3_patch_put(m3);
    }

    H5Gclose(group_fld);
  }

  H5Gclose(group0);
}


// ----------------------------------------------------------------------
// mrc_io_ops_xdmf

static struct mrc_io_ops mrc_io_ops_xdmf2 = {
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

void
libmrc_io_register_xdmf2()
{
  libmrc_io_register(&mrc_io_ops_xdmf2);
  //  libmrc_io_register(&mrc_io_ops_xdmf_to_one);
#ifdef H5_HAVE_PARALLEL
  //  libmrc_io_register(&mrc_io_ops_xdmf_parallel);
#endif
}

