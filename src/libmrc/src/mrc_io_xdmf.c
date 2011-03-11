
#include <mrc_io_private.h>
#include <mrc_params.h>
#include <mrc_list.h>

// TODO:

// make vectors for slices
// H5LT
// simplify? support multiple domains?
// combine into fewer files..
// don't write XDMF multi-block if it's a single block
// use hyperslabs, eventually
// cell vs node

//FIXME
#ifdef HAVE_HDF5_H

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#define CE assert(ierr == 0)

#define MAX_FLD_INFO (30)

struct fld_info {
  char *name;
  char *path;
  bool is_vec;
  int dim;
};

struct xdmf_temporal_step {
  list_t entry;
  char filename[0];
};

struct xdmf_temporal {
  char *filename;
  list_t timesteps;
};

struct xdmf_subdomain {
  int im[3];
};

struct xdmf_spatial {
  char *sfx;
  float sheet;
  int p;
  int nr_subdomains;
  struct xdmf_subdomain *subdomains;
  int nr_fld_info;
  struct fld_info fld_info[MAX_FLD_INFO];
  void (*write_topology)(FILE *f, int im[3], const char *filename, float sheet, int p);
  void (*write_fld)(FILE *f, struct fld_info *fld_info, int im[3], const char *filename);
  list_t entry;
};

struct diag_hdf5_params {
  bool use_independent_io;
};

#define VAR(x) (void *)offsetof(struct diag_hdf5_params, x)

static struct param diag_hdf5_params_descr[] = {
  { "use_independent_io"     , VAR(use_independent_io)      , PARAM_BOOL(false)      },
  {},
};

#undef VAR

struct diag_hdf5 {
  struct diag_hdf5_params par;
  char *mode;
  hid_t file;
  bool crd_written[3]; // this should be per mrc_domain, but this will do for crds
  bool crds_done;
  int gdims[3], nr_procs[3];
  int off[3], ldims[3];
  struct mrc_f3 *vfld;
  struct xdmf_temporal *xdmf_temporal; // lives from create() to destroy()
  list_t xdmf_spatial_list;  // lives from open to close
  hid_t group_crd;
};

// ======================================================================
// filename helpers

static void
make_path(char *path, const char *sfx)
{
  sprintf(path, "%s", sfx);
}

static void
make_sfx(char *sfx, int outtype, float sheet)
{
  int out_sheet = (sheet>=0.0) ? (int)(sheet*10.0001): -(int)(-sheet*10.0001);
  if (outtype == DIAG_TYPE_2D_X) {
    sprintf(sfx, "px_%d", out_sheet);
  } else if (outtype == DIAG_TYPE_2D_Y) {
    sprintf(sfx, "py_%d", out_sheet);
  } else if (outtype == DIAG_TYPE_2D_Z) {
    sprintf(sfx, "pz_%d", out_sheet);
  } else if (outtype == DIAG_TYPE_2D_IONO) {
    sprintf(sfx, "iof");
  } else {
    fprintf(stderr, "unknown outtype %d!\n", outtype);
    assert(0);
  }
}

// ======================================================================
// generic hdf5

#define diag_hdf5(io) ((struct diag_hdf5 *)io->obj.subctx)

static void
hdf5_open(struct mrc_io *io) 
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);
  char filename[strlen(io->par.outdir) + strlen(io->par.basename) + 20];
  sprintf(filename, "%s/%s.%06d_p%06d.h5", io->par.outdir, io->par.basename,
	  io->step, io->rank);
  hdf5->file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

static void
hdf5_close(struct mrc_io *io) 
{
  herr_t ierr;
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

  hid_t group = H5Gcreate(hdf5->file, "mrc_io", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_int(group, ".", "step", &io->step, 1);
  H5LTset_attribute_float(group, ".", "time", &io->time, 1);
  H5Gclose(group);
  
  group = H5Gcreate(hdf5->file, "mrc_info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_int(group, ".", "global_dims", hdf5->gdims, 3);
  
  H5LTset_attribute_int(group, ".", "proc_rank", &io->rank, 1);
  // put in the number of procs and global indices
  // Note that the nproc variables are the number of procs writing output,
  // not the number of procs for the simulation.  
  // Simlarly for the global_?_idx+{min,max}! 
  
  int g_high[3];
  for (int d = 0; d < 3; d++) {
    g_high[d] = hdf5->off[d] + hdf5->ldims[d];
  }
  H5LTset_attribute_int(group, ".", "nproc", hdf5->nr_procs, 3);
  H5LTset_attribute_int(group, ".", "global_idx_low", hdf5->off, 3);
  H5LTset_attribute_int(group, ".", "global_idx_high", g_high, 3);

  ierr = H5Gclose(group); CE;

  ierr = H5Fclose(hdf5->file); CE;
}

static void
hdf5_write_field2d_serial(struct mrc_io *io, float scale, struct mrc_f2 *fld,
			  const char *path)
{
  herr_t ierr;
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

  hsize_t fdims[2] = { fld->im[1], fld->im[0] };
  //  printf("[%d] diagsrv: write '%s'\n", info->rank, fld_name);

  hid_t group0 = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  hid_t group = H5Gcreate(group0, fld->name[0], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // FIXME H5lt
  hid_t filespace = H5Screate_simple(2, fdims, NULL);
  hid_t dataset = H5Dcreate(group, "2d", H5T_NATIVE_FLOAT, filespace,
			    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  ierr = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, fld->arr); CE;

  ierr = H5Dclose(dataset); CE;
  ierr = H5Sclose(filespace); CE;
  ierr = H5Gclose(group); CE;
  ierr = H5Gclose(group0); CE;
}

// ======================================================================
// XDMF
//
// works in two modes:
// -- every procs writes its own HDF5, and one proc creates an XDMF
//    description file describing the global topology
// -- use diagsrv, so only one proc writes, but still create XDMF description

static void
xdmf_write_header(FILE *f)
{
  fprintf(f, "<?xml version='1.0' ?>\n");
  fprintf(f, "<Xdmf xmlns:xi='http://www.w3.org/2001/XInclude' Version='2.0'>\n");
}

static void
xdmf_write_topology_3d(FILE *f, int im[3], const char *filename, float sheet, int p)
{
  fprintf(f, "     <Topology TopologyType=\"3DRectMesh\" Dimensions=\"%d %d %d\"/>\n",
	  im[2] + 1, im[1] + 1, im[0] + 1);
  fprintf(f, "     <Geometry GeometryType=\"VXVYVZ\">\n");
  fprintf(f, "     <DataItem Name=\"VX\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[0] + 1);
  if (p >= 0) {
    fprintf(f, "        ./%s:/crd/x-%d\n", filename, p);
  } else {
    fprintf(f, "        ./%s:/crd/x\n", filename);
  }
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VY\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[1] + 1);
  if (p >= 0) {
    fprintf(f, "        ./%s:/crd/y-%d\n", filename, p);
  } else {
    fprintf(f, "        ./%s:/crd/y\n", filename);
  }
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VZ\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[2] + 1);
  if (p >= 0) {
    fprintf(f, "        ./%s:/crd/z-%d\n", filename, p);
  } else {
    fprintf(f, "        ./%s:/crd/z\n", filename);
  }
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     </Geometry>\n");
  fprintf(f, "\n");
}

static void
xdmf_write_topology_2d_x(FILE *f, int im[2], const char *filename, float sheet, int p)
{
  fprintf(f, "     <Topology TopologyType=\"3DRectMesh\" Dimensions=\"%d %d %d\"/>\n",
	  im[1] + 1, im[0] + 1, 2);
  fprintf(f, "     <Geometry GeometryType=\"VXVYVZ\">\n");
  fprintf(f, "     <DataItem Name=\"VX\" DataType=\"Float\" Dimensions=\"%d\">\n", 2);
  fprintf(f, "        %g %g\n", sheet - .01, sheet + .01);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VY\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[0] + 1);
  fprintf(f, "        ./%s:/crd/y\n", filename);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VZ\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[1] + 1);
  fprintf(f, "        ./%s:/crd/z\n", filename);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     </Geometry>\n");
  fprintf(f, "\n");
}

static void
xdmf_write_topology_2d_y(FILE *f, int im[2], const char *filename, float sheet, int p)
{
  fprintf(f, "     <Topology TopologyType=\"3DRectMesh\" Dimensions=\"%d %d %d\"/>\n",
	  im[1] + 1, 2, im[0] + 1);
  fprintf(f, "     <Geometry GeometryType=\"VXVYVZ\">\n");
  fprintf(f, "     <DataItem Name=\"VX\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[0] + 1);
  fprintf(f, "        ./%s:/crd/x\n", filename);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VY\" DataType=\"Float\" Dimensions=\"%d\">\n", 2);
  fprintf(f, "        %g %g\n", sheet - .01, sheet + .01);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VZ\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[1] + 1);
  fprintf(f, "        ./%s:/crd/z\n", filename);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     </Geometry>\n");
  fprintf(f, "\n");
}

static void
xdmf_write_topology_2d_z(FILE *f, int im[2], const char *filename, float sheet, int p)
{
  fprintf(f, "     <Topology TopologyType=\"3DRectMesh\" Dimensions=\"%d %d %d\"/>\n",
	  2, im[1] + 1, im[0] + 1);
  fprintf(f, "     <Geometry GeometryType=\"VXVYVZ\">\n");
  fprintf(f, "     <DataItem Name=\"VX\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[0] + 1);
  fprintf(f, "        ./%s:/crd/x\n", filename);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VY\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[1] + 1);
  fprintf(f, "        ./%s:/crd/y\n", filename);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VZ\" DataType=\"Float\" Dimensions=\"%d\">\n", 2);
  fprintf(f, "        %g %g\n", sheet - .01, sheet + .01);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     </Geometry>\n");
  fprintf(f, "\n");
}

static void
xdmf_write_topology_iono(FILE *f, int im[2], const char *filename, float sheet, int p)
{
  fprintf(f, "     <Topology TopologyType=\"2DSMesh\" Dimensions=\"%d %d\"/>\n", im[1], im[0]);
  fprintf(f, "     <Geometry GeometryType=\"XYZ\">\n");
  fprintf(f, "       <DataItem Format=\"XML\" Dimensions=\"%d %d 3\">\n", im[1], im[0]);
  // FIXME, write to HDF5
  float dphi = 2.*M_PI / (im[0] - 1);
  float dtheta = M_PI / (im[1] - 1);
  for (int itheta = 0; itheta < im[1]; itheta++) {
    for (int iphi = 0; iphi < im[0]; iphi++) {
      double phi = iphi * dphi, theta = itheta * dtheta;
      fprintf(f, "        %g %g %g\n", cos(phi)*sin(theta), sin(phi)*sin(theta),cos(theta)
);
    }
  }
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Geometry>\n");
  fprintf(f, "\n");
}

static void
xdmf_write_fld_2d_iono(FILE *f, struct fld_info *fld_info, int im[2], const char *filename)
{
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",
	  fld_info->name);
  fprintf(f, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
	  im[1], im[0]);
  fprintf(f, "        %s:/%s/%s/2d\n", filename, fld_info->path, fld_info->name);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
}

static void
xdmf_write_fld_2d_x(FILE *f, struct fld_info *fld_info, int im[2], const char *filename)
{
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
	  fld_info->name);
  fprintf(f, "       <DataItem Dimensions=\"1 %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
	  im[1], im[0]);
  fprintf(f, "        %s:/%s/%s/2d\n", filename, fld_info->path, fld_info->name);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
}

static void
xdmf_write_fld_2d_y(FILE *f, struct fld_info *fld_info, int im[2], const char *filename)
{
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
	  fld_info->name);
  fprintf(f, "       <DataItem Dimensions=\"%d 1 %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
	  im[1], im[0]);
  fprintf(f, "        %s:/%s/%s/2d\n", filename, fld_info->path, fld_info->name);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
}

static void
xdmf_write_fld_2d_z(FILE *f, struct fld_info *fld_info, int im[2], const char *filename)
{
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
	  fld_info->name);
  fprintf(f, "       <DataItem Dimensions=\"%d %d 1\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
	  im[1], im[0]);
  fprintf(f, "        %s:/%s/%s/2d\n", filename, fld_info->path, fld_info->name);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
}

static void
xdmf_write_fld_3d(FILE *f, struct fld_info *fld_info, int im[3], const char *filename)
{
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Cell\">\n",
	  fld_info->name, fld_info->is_vec ? "Vector" : "Scalar");
  fprintf(f, "       <DataItem Dimensions=\"%d %d %d%s\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
	  im[2], im[1], im[0], fld_info->is_vec ? " 3" : "");
  fprintf(f, "        %s:/%s/%s/3d\n", filename, fld_info->path, fld_info->name);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
}

static struct xdmf_spatial *
xdmf_spatial_create(struct mrc_io *io, const char *sfx, float sheet, int p, int nr_subdomains,
		    struct xdmf_subdomain *subdomains, int size,
		    void (*write_topology)(FILE *f, int im[3], const char *filename, float sheet, int p),
		    void (*write_fld)(FILE *f, struct fld_info *fld_info, int im[3],
				      const char *filename))
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);
  struct xdmf_spatial *xs = malloc(sizeof(*xs));
  memset(xs, 0, sizeof(*xs));
  xs->sfx            = strdup(sfx);
  xs->sheet          = sheet;
  xs->nr_subdomains  = nr_subdomains;
  xs->subdomains     = subdomains;
  xs->write_topology = write_topology;
  xs->write_fld      = write_fld;
  
  list_add_tail(&xs->entry, &hdf5->xdmf_spatial_list);

  return xs;
}

static struct xdmf_spatial *
xdmf_spatial_find(struct mrc_io *io, const char *sfx, float sheet)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);
  struct xdmf_spatial *xs;
  list_for_each_entry(xs, &hdf5->xdmf_spatial_list, entry) {
    if (xs->sheet == sheet && strcmp(xs->sfx, sfx) == 0) {
      return xs;
    }
  }
  return NULL;
}

static void
xdmf_spatial_destroy(struct xdmf_spatial *xs)
{
  list_del(&xs->entry);

  for (int i = 0; i < xs->nr_fld_info; i++) {
    free(xs->fld_info[i].path);
    free(xs->fld_info[i].name);
  }
  free(xs->subdomains);
  free(xs->sfx);
  free(xs);
}

static void
xdmf_write_spatial_collection(struct mrc_io *io, struct xdmf_spatial *xs, const char *filename)
{
  if (io->rank > 0)
    return;

  mprintf("diag: writing '%s'\n", filename);

  FILE *f = fopen(filename, "w");
  xdmf_write_header(f);
  fprintf(f, "<Domain>\n");
  fprintf(f, "<Grid GridType=\"Collection\" CollectionType=\"Spatial\">\n");
  fprintf(f, "   <Time Type=\"Single\" Value=\"%g\" />\n", io->time);
  for (int s = 0; s < xs->nr_subdomains; s++) {
    fprintf(f, "   <Grid Name=\"mesh-%s-%d\" GridType=\"Uniform\">\n", xs->sfx, s);
    
    char fname[strlen(io->par.outdir) + strlen(io->par.basename) + 20];
    sprintf(fname, "%s/%s.%06d_p%06d.h5", io->par.outdir, io->par.basename, io->step, s);
    int *im = xs->subdomains[s].im;
    xs->write_topology(f, im, fname, xs->sheet, xs->p);
    
    for (int m = 0; m < xs->nr_fld_info; m++) {
      xs->write_fld(f, &xs->fld_info[m], im, fname);
    }
    fprintf(f, "   </Grid>\n");
  }
  fprintf(f, "</Grid>\n");
  fprintf(f, "</Domain>\n");
  fprintf(f, "</Xdmf>\n");
  fclose(f);
}

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

static void
xdmf_open(struct mrc_io *io) 
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);
  hdf5_open(io);

  hdf5->group_crd = H5Gcreate(hdf5->file, "crd", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for (int d = 0; d < 3; d++) {
    hdf5->crd_written[d] = false;
  }
  hdf5->crds_done = false;
  INIT_LIST_HEAD(&hdf5->xdmf_spatial_list);
}

static void
xdmf_close(struct mrc_io *io)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);
  while (!list_empty(&hdf5->xdmf_spatial_list)) {
    struct xdmf_spatial *xs = list_entry(hdf5->xdmf_spatial_list.next, typeof(*xs), entry);
    
    char fname_spatial[strlen(io->par.outdir) + strlen(io->par.basename) + 20];
    sprintf(fname_spatial, "%s/%s.%06d.xdmf", io->par.outdir, io->par.basename, io->step);
    xdmf_write_spatial_collection(io, xs, fname_spatial); // FIXME, rm dir part
    xdmf_spatial_destroy(xs);

    if (io->rank == 0) {
      struct xdmf_temporal *xt = hdf5->xdmf_temporal;
      xdmf_temporal_append(xt, fname_spatial);
      xdmf_temporal_write(xt);
    }
  }
}

static void
ds_xdmf_setup(struct mrc_io *io)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);
  char filename[strlen(io->par.outdir) + strlen(io->par.basename) + 7];
  sprintf(filename, "%s/%s.xdmf", io->par.outdir, io->par.basename);
  hdf5->xdmf_temporal = xdmf_temporal_create(filename);
}

static void
ds_xdmf_destroy(struct mrc_io *io)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

  if (hdf5->xdmf_temporal) {
    xdmf_temporal_destroy(hdf5->xdmf_temporal);
  }
}

static void
ds_xdmf_open(struct mrc_io *io, const char *mode) 
{
  assert(strcmp(mode, "w") == 0);

  xdmf_open(io);
}

static void
hdf5_write_crds(struct mrc_io *io, int im[3], struct mrc_domain *domain, int sw)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

  if (!hdf5->gdims[0]) {
    mrc_domain_get_global_dims(domain, hdf5->gdims);
    mrc_domain_get_nr_procs(domain, hdf5->nr_procs);
    int nr_patches;
    struct mrc_patch *patches = mrc_domain_get_patches(domain, &nr_patches);
    assert(nr_patches == 1);
    for (int d = 0; d < 3; d++) {
      hdf5->off[d] = patches[0].off[d];
      hdf5->ldims[d] = patches[0].ldims[d];
    }
  }

  const char *xyz[3] = { "x", "y", "z" };
  
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  for (int d = 0; d < 3; d++) {
    if (im[d] == 0 || hdf5->crd_written[d])
      continue;

    float *crd = &MRC_CRD(crds, d, sw);
    float *crd_nc = calloc(im[d] + 1, sizeof(*crd_nc));
    if (sw > 0) {
      for (int i = 0; i <= im[d]; i++) {
	crd_nc[i] = .5 * (crd[i-1] + crd[i]);
      }
    } else {
      for (int i = 1; i < im[d]; i++) {
	crd_nc[i] = .5 * (crd[i-1] + crd[i]);
      }
      // extrapolate
      crd_nc[0]  = crd[0] - .5 * (crd[1] - crd[0]);
      crd_nc[im[d]] = crd[im[d]-1] + .5 * (crd[im[d]-1] - crd[im[d]-2]);
    }
    hsize_t im1 = im[d] + 1;
    H5LTmake_dataset_float(hdf5->group_crd, xyz[d], 1, &im1, crd_nc);
    free(crd_nc);
    hdf5->crd_written[d] = true;
  }
}

static void
hdf5_write_mcrds(struct mrc_io *io, struct mrc_domain *domain, int sw)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

  if (!hdf5->gdims[0]) {
    mrc_domain_get_global_dims(domain, hdf5->gdims);
    int nr_patches;
    struct mrc_patch *patches = mrc_domain_get_patches(domain, &nr_patches);
    assert(nr_patches == 1);
    for (int d = 0; d < 3; d++) {
      hdf5->off[d] = patches[0].off[d];
      hdf5->ldims[d] = patches[0].ldims[d];
    }
  }

  
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  for (int d = 0; d < 3; d++) {
    if (hdf5->crd_written[d])
      continue;

    struct mrc_m1 *mcrd = crds->mcrd[d];
    mrc_m1_foreach_patch(mcrd, p) {
      struct mrc_m1_patch *mcrdp = mrc_m1_patch_get(mcrd, p);
      int im = mcrdp->im[0];
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
      sprintf(name, "%c-%d", 'x' + d, p);
      H5LTmake_dataset_float(hdf5->group_crd, name, 1, &im1, crd_nc);

      free(crd_nc);
      mrc_m1_patch_put(mcrd);
    }
    hdf5->crd_written[d] = true;
  }
}

static void
ds_xdmf_close(struct mrc_io *io)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

  H5Gclose(hdf5->group_crd);
  hdf5_close(io);
  xdmf_close(io);
}

static void
ds_xdmf_read_attr(struct mrc_io *io, const char *path, int type,
		  const char *name, union param_u *pv)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);
  
  hid_t group = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  switch (type) {
  case PT_SELECT:
  case PT_INT:
    H5LTget_attribute_int(group, ".", name, &pv->u_int);
    break;
  case PT_BOOL: {
    int val;
    H5LTget_attribute_int(group, ".", name, &val);
    pv->u_bool = val;
    break;
  }
  case PT_FLOAT:
    H5LTget_attribute_float(group, ".", name, &pv->u_float);
    break;
  case PT_DOUBLE:
    H5LTget_attribute_double(group, ".", name, &pv->u_double);
    break;
  case PT_STRING: ;
    hsize_t dims;
    H5T_class_t class;
    size_t sz;
    H5LTget_attribute_info(group, ".", name, &dims, &class, &sz);
    pv->u_string = malloc(sz);
    H5LTget_attribute_string(group, ".", name, (char *)pv->u_string);
    break;
  case PT_INT3:
    H5LTget_attribute_int(group, ".", name, pv->u_int3);
    break;
  case PT_FLOAT3:
    H5LTget_attribute_float(group, ".", name, pv->u_float3);
    break;
  }
  H5Gclose(group);
}

static void
ds_xdmf_write_attr(struct mrc_io *io, const char *path, int type,
		   const char *name, union param_u *pv)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);
  
  hid_t group;
  if (H5Lexists(hdf5->file, path, H5P_DEFAULT) > 0) {
    group = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  } else {
    group = H5Gcreate(hdf5->file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
copy_and_scale(struct mrc_f3 *vfld, int m, struct mrc_f3 *fld, int fld_m,
	       float scale, int bnd)
{
  for (int iz = 0; iz < vfld->im[2]; iz++) {
    for (int iy = 0; iy < vfld->im[1]; iy++) {
      for (int ix = 0; ix < vfld->im[0]; ix++) {
	// cannot use MRC_F3, because the layout is different (for vecs, fast component)!
	vfld->arr[(((iz * vfld->im[1]) + iy) * vfld->im[0] + ix) * vfld->nr_comp + m] =
	  scale * MRC_F3(fld, fld_m, ix+bnd,iy+bnd,iz+bnd);
      }
    }
  }
}

static void
copy_back(struct mrc_f3 *vfld, int m, struct mrc_f3 *fld, int fld_m, int bnd)
{
  for (int iz = 0; iz < vfld->im[2]; iz++) {
    for (int iy = 0; iy < vfld->im[1]; iy++) {
      for (int ix = 0; ix < vfld->im[0]; ix++) {
	// cannot use MRC_F3, because the layout is different (for vecs, fast component)!
	MRC_F3(fld, fld_m, ix+bnd,iy+bnd,iz+bnd) = 
	  vfld->arr[(((iz * vfld->im[1]) + iy) * vfld->im[0] + ix) * vfld->nr_comp + m];
      }
    }
  }
}

static void
save_fld_info(struct xdmf_spatial *xs, char *fld_name, char *path, bool is_vec)
{
  assert(xs->nr_fld_info < MAX_FLD_INFO);
  struct fld_info *fld_info = &xs->fld_info[xs->nr_fld_info++];

  fld_info->name = fld_name;
  fld_info->path = path;
  fld_info->is_vec = is_vec;
}

static struct xdmf_spatial *
xdmf_spatial_create_3d(struct mrc_io *io, int im[3], int p, int nr_subdomains)
{
  // OPT, we could skip this on procs which aren't writing xdmf

  struct xdmf_subdomain *subdomains = calloc(nr_subdomains, sizeof(*subdomains));
  for (int s = 0; s < nr_subdomains; s++) {
    for (int d = 0; d < 3; d++) {
      subdomains[s].im[d] = im[d];
    }
  }

  return xdmf_spatial_create(io, "3df", -1, p, nr_subdomains, subdomains, io->size,
			     xdmf_write_topology_3d, xdmf_write_fld_3d);
}

static struct xdmf_spatial *
xdmf_spatial_create_2d_x(struct mrc_io *io, int im[2], 
			 const char *sfx, int sheet, int nr_subdomains)
{
  // OPT, we could skip this on procs which aren't writing xdmf

  struct xdmf_subdomain *subdomains = calloc(nr_subdomains, sizeof(*subdomains));
  for (int s = 0; s < nr_subdomains; s++) {
    for (int d = 0; d < 2; d++) {
      subdomains[s].im[d] = im[d];
    }
  }

  return xdmf_spatial_create(io, sfx, sheet, -1, nr_subdomains, subdomains, io->size,
			     xdmf_write_topology_2d_x, xdmf_write_fld_2d_x);
}

static struct xdmf_spatial *
xdmf_spatial_create_2d_y(struct mrc_io *io, int im[2], 
			 const char *sfx, int sheet, int nr_subdomains)
{
  // OPT, we could skip this on procs which aren't writing xdmf

  struct xdmf_subdomain *subdomains = calloc(nr_subdomains, sizeof(*subdomains));
  for (int s = 0; s < nr_subdomains; s++) {
    for (int d = 0; d < 2; d++) {
      subdomains[s].im[d] = im[d];
    }
  }

  return xdmf_spatial_create(io, sfx, sheet, -1, nr_subdomains, subdomains, io->size,
			     xdmf_write_topology_2d_y, xdmf_write_fld_2d_y);
}

static struct xdmf_spatial *
xdmf_spatial_create_2d_z(struct mrc_io *io, int im[2], 
			 const char *sfx, int sheet, int nr_subdomains)
{
  // OPT, we could skip this on procs which aren't writing xdmf

  struct xdmf_subdomain *subdomains = calloc(nr_subdomains, sizeof(*subdomains));
  for (int s = 0; s < nr_subdomains; s++) {
    for (int d = 0; d < 2; d++) {
      subdomains[s].im[d] = im[d];
    }
  }

  return xdmf_spatial_create(io, sfx, sheet, -1, nr_subdomains, subdomains, io->size,
			     xdmf_write_topology_2d_z, xdmf_write_fld_2d_z);
}

static struct xdmf_spatial *
xdmf_spatial_create_iono(struct mrc_io *io, int im[2])
{
  struct xdmf_subdomain *subdomain = calloc(1, sizeof(*subdomain));
  for (int d = 0; d < 2; d++) {
    subdomain->im[d] = im[d];
  }

  return xdmf_spatial_create(io, "iof", -1, -1, 1, subdomain, 1,
			     xdmf_write_topology_iono, xdmf_write_fld_2d_iono);
}

static void
ds_xdmf_write_field(struct mrc_io *io, const char *path,
		    float scale, struct mrc_f3 *fld, int m)
{
  herr_t ierr;
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

  // on diagsrv, ghost points have already been dropped
  int sw = fld->sw;
  int im[3] = { fld->im[0] - 2 * sw, fld->im[1] - 2 * sw, fld->im[2] - 2 * sw };
  struct xdmf_spatial *xs;
  xs = xdmf_spatial_find(io, "3df", -1);
  if (!xs) {
    xs = xdmf_spatial_create_3d(io, im, -1, io->size);
    hdf5_write_crds(io, im, fld->domain, fld->sw);
  }

  hid_t group0 = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  char c = fld->name[m][strlen(fld->name[m]) - 1];
  if (c == 'x') {
    assert(!hdf5->vfld);
    hdf5->vfld = mrc_f3_alloc(mrc_f3_comm(fld), NULL, im);
    mrc_f3_set_nr_comps(hdf5->vfld, 3);
    mrc_f3_setup(hdf5->vfld);
    copy_and_scale(hdf5->vfld, 0, fld, m, scale, sw);
  } else if (c == 'y') {
    copy_and_scale(hdf5->vfld, 1, fld, m, scale, sw);
  } else if (c == 'z') {
    copy_and_scale(hdf5->vfld, 2, fld, m, scale, sw);
    char *vec_name = strdup(fld->name[m]);
    vec_name[strlen(fld->name[m])-1] = 0;
    save_fld_info(xs, vec_name, strdup(path), true);
    struct mrc_f3 *vfld = hdf5->vfld;
    hsize_t hdims[4] = { vfld->im[2], vfld->im[1], vfld->im[0], vfld->nr_comp };
    hid_t group = H5Gcreate(group0, vec_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    ierr = H5LTmake_dataset_float(group, "3d", 4, hdims, vfld->arr); CE;
    ierr = H5Gclose(group); CE;

    mrc_f3_destroy(vfld);
    hdf5->vfld = NULL;
  } else { // scalar
    save_fld_info(xs, strdup(fld->name[m]), strdup(path), false);
    if (io->rank < 0) { // FIXME not happening anymore, still optimize this case
      // on diag srv, ghost points are already gone and scaling is done
      assert(scale == 1.f);
      hsize_t hdims[3] = { fld->im[2], fld->im[1], fld->im[0] };
      hid_t group = H5Gcreate(group0, fld->name[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      ierr = H5LTmake_dataset_float(group, "3d", 3, hdims, &MRC_F3(fld, m, 0,0,0)); CE;
      ierr = H5Gclose(group); CE;
    } else {
      struct mrc_f3 *sfld = mrc_f3_alloc(mrc_f3_comm(fld), NULL, im);
      mrc_f3_setup(sfld);
      copy_and_scale(sfld, 0, fld, m, scale, sw);

      hsize_t hdims[3] = { sfld->im[2], sfld->im[1], sfld->im[0] };
      hid_t group = H5Gcreate(group0, fld->name[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      ierr = H5LTmake_dataset_float(group, "3d", 3, hdims, sfld->arr); CE;
      ierr = H5Gclose(group); CE;
      mrc_f3_destroy(sfld);
    }
  }
  H5Gclose(group0);
}

static void
ds_xdmf_write_field2d(struct mrc_io *io, float scale, struct mrc_f2 *fld,
		      int outtype, float sheet)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);
  int ierr;

  // diagsrv xdmf_serial / single proc only for now
  int bnd = 0;
  int im[2] = { fld->im[0] - bnd, fld->im[1] - bnd };

  struct xdmf_spatial *xs = NULL;
  char sfx[10]; make_sfx(sfx, outtype, sheet);
  char path[100]; make_path(path, sfx);
  xs = xdmf_spatial_find(io, sfx, sheet);
  if (!xs) {
    if (outtype == DIAG_TYPE_2D_X) {
      xs = xdmf_spatial_create_2d_x(io, im, sfx, sheet, io->size);
      hdf5_write_crds(io, (int[]) { 0, im[0], im[1] }, fld->domain, fld->sw);
    } else if (outtype == DIAG_TYPE_2D_Y) {
      xs = xdmf_spatial_create_2d_y(io, im, sfx, sheet, io->size);
      hdf5_write_crds(io, (int[]) { im[0], 0, im[1] }, fld->domain, fld->sw);
    } else if (outtype == DIAG_TYPE_2D_Z) {
      xs = xdmf_spatial_create_2d_z(io, im, sfx, sheet, io->size);
      hdf5_write_crds(io, (int[]) { im[0], im[1], 0 }, fld->domain, fld->sw);
    } else if (outtype == DIAG_TYPE_2D_IONO) {
      xs = xdmf_spatial_create_iono(io, im);
    }
    hid_t group = H5Gcreate(hdf5->file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    ierr = H5Gclose(group); CE;
  }

  save_fld_info(xs, strdup(fld->name[0]), strdup(path), false);
  hdf5_write_field2d_serial(io, scale, fld, path);
}

static void
ds_xdmf_write_m3(struct mrc_io *io, const char *path, struct mrc_m3 *m3)
{
  int ierr;
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

  hid_t group0 = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  H5LTset_attribute_int(group0, ".", "nr_patches", &m3->nr_patches, 1);

  mrc_m3_foreach_patch(m3, p) {
    struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, p);

    struct xdmf_spatial *xs = xdmf_spatial_find(io, "3df", p);
    if (!xs) {
      int np[3];
      mrc_domain_get_param_int3(m3->domain, "np", np);
      xs = xdmf_spatial_create_3d(io, m3p->im, p, np[0] * np[1] * np[2]);
      if (p == 0) {
	hdf5_write_mcrds(io, m3->domain, m3->sw);
      }
    }

    for (int m = 0; m < m3->nr_comp; m++) {
      char fld_name[strlen(m3->name[m]) + 5];

      sprintf(fld_name, "%s-%d", m3->name[m], p);
      mprintf("fld_name %s\n", fld_name);

      save_fld_info(xs, strdup(fld_name), strdup(path), false);
      hsize_t hdims[3] = { m3p->im[2], m3p->im[1], m3p->im[0] };
      hid_t group = H5Gcreate(group0, fld_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      ierr = H5LTmake_dataset_float(group, "3d", 3, hdims, m3p->arr); CE;
      ierr = H5Gclose(group); CE;
      mrc_m3_patch_put(m3);
    }
  }

  H5Gclose(group0);
}

struct mrc_io_ops mrc_io_xdmf_ops = {
  .name          = "xdmf",
  .size          = sizeof(struct diag_hdf5),
  .parallel      = true,
  .destroy       = ds_xdmf_destroy,
  .setup         = ds_xdmf_setup,
  .open          = ds_xdmf_open,
  .close         = ds_xdmf_close,
  .write_field   = ds_xdmf_write_field,
  .write_field2d = ds_xdmf_write_field2d,
  .write_attr    = ds_xdmf_write_attr,
  .write_m3      = ds_xdmf_write_m3,
};

struct mrc_io_ops mrc_io_xdmf_serial_ops = {
  .name          = "xdmf_serial",
  .size          = sizeof(struct diag_hdf5),
  .destroy       = ds_xdmf_destroy,
  .setup         = ds_xdmf_setup,
  .open          = ds_xdmf_open,
  .close         = ds_xdmf_close,
  .write_field   = ds_xdmf_write_field,
  .write_field2d = ds_xdmf_write_field2d,
  .write_attr    = ds_xdmf_write_attr,
};

// ======================================================================

enum {
  TAG_OFF_DIMS = 150000,
  TAG_DATA,
  TAG_CRDX,
  TAG_CRDY,
  TAG_CRDZ,
};

static void
ds_xdmf_to_one_open(struct mrc_io *io, const char *mode) 
{
  assert(strcmp(mode, "w") == 0);
#ifndef NDEBUG
  MPI_Barrier(io->obj.comm);
#endif
  
  struct diag_hdf5 *hdf5 = diag_hdf5(io);
  hdf5->crds_done = false;
  if (io->rank != 0)
    return;

  xdmf_open(io);
}

static void
ds_xdmf_to_one_close(struct mrc_io *io)
{
#ifndef NDEBUG
  MPI_Barrier(io->obj.comm);
#endif
  
  if (io->rank != 0)
    return;

  struct diag_hdf5 *hdf5 = diag_hdf5(io);

  H5Gclose(hdf5->group_crd);
  hdf5_close(io);
}

static void
ds_xdmf_to_one_write_attr(struct mrc_io *io, const char *path, int type,
			 const char *name, union param_u *pv)
{
  if (io->rank == 0) {
    ds_xdmf_write_attr(io, path, type, name, pv);
  }
}

static void
communicate_crds(struct mrc_io *io, struct mrc_f3 *gfld, struct mrc_f3 *lfld)
{
  int sw = gfld->sw;

  struct mrc_domain *gdomain = gfld->domain;
  struct mrc_crds *gcrds = mrc_domain_get_crds(gdomain);
  if (io->rank != 0) {
    int iw[6], *ib = iw, *im = iw + 3;
    int nr_patches;
    struct mrc_patch *patches = mrc_domain_get_patches(gdomain, &nr_patches);
    assert(nr_patches == 1);
    for (int d = 0; d < 3; d++) {
      ib[d] = patches[0].off[d];
      im[d] = patches[0].ldims[d];
    }
    MPI_Send(iw, 6, MPI_INT, 0, TAG_OFF_DIMS, io->obj.comm);
    for (int d = 0; d < 3; d++) {
      MPI_Send(&MRC_CRD(gcrds, d, sw), im[d], MPI_FLOAT, 0,
	       TAG_CRDX + d, io->obj.comm);
    }
  } else { // io->rank == 0
    struct mrc_domain *ldomain = lfld->domain;
    struct mrc_crds *lcrds = mrc_domain_get_crds(ldomain);
    for (int n = 0; n < io->size; n++) {
      float *recv_crds[3];
      int iw[6], *ib = iw, *im = iw + 3;
      if (n == 0) {
	int nr_patches;
	struct mrc_patch *patches = mrc_domain_get_patches(gdomain, &nr_patches);
	assert(nr_patches == 1);
	for (int d = 0; d < 3; d++) {
	  ib[d] = patches[0].off[d];
	  im[d] = patches[0].ldims[d];
	}
	for (int d = 0; d < 3; d++) {
	  recv_crds[d] = &MRC_CRD(gcrds, d, sw);
	}
      } else {
	MPI_Recv(iw, 6, MPI_INT, n, TAG_OFF_DIMS, io->obj.comm, MPI_STATUS_IGNORE);
	for (int d = 0; d < 3; d++) {
	  recv_crds[d] = calloc(im[d], sizeof(*recv_crds[d]));
	  MPI_Recv(recv_crds[d], im[d], MPI_FLOAT, n, TAG_CRDX + d, io->obj.comm,
		   MPI_STATUS_IGNORE);
	}
      }

      for (int d = 0; d < 3; d++) {
	for (int i = 0; i < im[d]; i++) {
	  MRC_CRD(lcrds, d, i + ib[d]) = recv_crds[d][i];
	}
      }
      
      if (n != 0) {
	for (int d = 0; d < 3; d++) {
	  free(recv_crds[d]);
	}
      }
    }
  }
}

static void
communicate_fld(struct mrc_io *io, struct mrc_f3 *gfld, int m, float scale,
		struct mrc_f3 *lfld)
{
  int sw = gfld->sw;

  struct mrc_f3 *send_fld = mrc_domain_f3_create(gfld->domain, SW_0);
  mrc_f3_setup(send_fld);
  copy_and_scale(send_fld, 0, gfld, m, scale, sw);

  if (io->rank != 0) {
    int iw[6], *ib = iw, *im = iw + 3;
    int nr_patches;
    struct mrc_patch *patches = mrc_domain_get_patches(send_fld->domain, &nr_patches);
    assert(nr_patches == 1);
    for (int d = 0; d < 3; d++) {
      ib[d] = patches[0].off[d];
      im[d] = patches[0].ldims[d];
    }
    MPI_Send(iw, 6, MPI_INT, 0, TAG_OFF_DIMS, io->obj.comm);
    MPI_Send(send_fld->arr, send_fld->len, MPI_FLOAT, 0, TAG_DATA, io->obj.comm);
  } else { // io->rank == 0
    for (int n = 0; n < io->size; n++) {
      struct mrc_f3 *recv_fld;
      if (n == 0) {
	recv_fld = mrc_f3_alloc(MPI_COMM_SELF, send_fld->ib, send_fld->im);
	mrc_f3_set_array(recv_fld, send_fld->arr);
	mrc_f3_setup(recv_fld);
      } else {
	int iw[6], *ib = iw, *im = iw + 3;
	MPI_Recv(iw, 6, MPI_INT, n, TAG_OFF_DIMS, io->obj.comm, MPI_STATUS_IGNORE);
	recv_fld = mrc_f3_alloc(MPI_COMM_SELF, ib, im);
	mrc_f3_setup(recv_fld);
	MPI_Recv(recv_fld->arr, recv_fld->len, MPI_FLOAT, n, TAG_DATA, io->obj.comm,
		 MPI_STATUS_IGNORE);
      }
      
      mrc_f3_foreach(recv_fld, ix,iy,iz, 0, 0) {
	MRC_F3(lfld,0, ix,iy,iz) = MRC_F3(recv_fld,0, ix,iy,iz);
      } mrc_f3_foreach_end;
      
      mrc_f3_destroy(recv_fld);
    }
  }

  mrc_f3_destroy(send_fld);
}

static void
ds_xdmf_to_one_write_field(struct mrc_io *io, const char *path,
			   float scale, struct mrc_f3 *fld, int m)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

#ifndef NDEBUG
  MPI_Barrier(io->obj.comm);
#endif
  int gdims[3];
  struct mrc_domain *ldomain = NULL;
  struct mrc_f3 *lfld = NULL;

  if (io->rank == 0) {
    mrc_domain_get_global_dims(fld->domain, gdims);
    ldomain = mrc_domain_create(MPI_COMM_SELF);
    mrc_domain_set_type(ldomain, "simple");
    mrc_domain_set_param_int3(ldomain, "m", gdims);
    struct mrc_crds *crds = mrc_domain_get_crds(ldomain);
    mrc_crds_set_type(crds, "rectilinear");
    mrc_domain_setup(ldomain);
    lfld = mrc_domain_f3_create(ldomain, SW_0);
    mrc_f3_setup(lfld);
  }

  communicate_fld(io, fld, m, scale, lfld);

  if (!hdf5->crds_done) {
    communicate_crds(io, fld, lfld);
    hdf5->crds_done = true;
  }

  if (io->rank != 0) {
    return;
  }

  struct xdmf_spatial *xs = xdmf_spatial_find(io, "3df", -1);
  if (!xs) {
    xs = xdmf_spatial_create_3d(io, gdims, -1, 1);
    hdf5_write_crds(io, gdims, ldomain, fld->sw);
  }
  save_fld_info(xs, strdup(fld->name[m]), strdup(path), false);

  hid_t group0 = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  hsize_t hdims[3] = { lfld->im[2], lfld->im[1], lfld->im[0] };
  hid_t group = H5Gcreate(group0, fld->name[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTmake_dataset_float(group, "3d", 3, hdims, lfld->arr);
  H5Gclose(group);
  H5Gclose(group0);

  mrc_f3_destroy(lfld);
}


struct mrc_io_ops mrc_io_xdmf_to_one_ops = {
  .name          = "xdmf_to_one",
  .size          = sizeof(struct diag_hdf5),
  .parallel      = true,
  .destroy       = ds_xdmf_destroy,
  .setup         = ds_xdmf_setup,
  .open          = ds_xdmf_to_one_open,
  .close         = ds_xdmf_to_one_close,
  .write_field   = ds_xdmf_to_one_write_field,
  .write_attr    = ds_xdmf_to_one_write_attr,
};

#ifdef H5_HAVE_PARALLEL

// ======================================================================
// xdmf_parallel

static void
ds_xdmf_parallel_open(struct mrc_io *io, const char *mode) 
{
#ifndef NDEBUG
  MPI_Barrier(io->obj.comm);
#endif
  
  struct diag_hdf5 *hdf5 = diag_hdf5(io);
  hdf5->mode = strdup(mode);
  hdf5->crds_done = false;

  char filename[strlen(io->par.outdir) + strlen(io->par.basename) + 20];
  sprintf(filename, "%s/%s.%06d_p%06d.h5", io->par.outdir, io->par.basename,
	  io->step, 0);

  hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist, io->obj.comm, MPI_INFO_NULL);
  if (strcmp(mode, "w") == 0) {
    hdf5->file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    hdf5->group_crd = H5Gcreate(hdf5->file, "crd", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  } else if (strcmp(mode, "r") == 0) {
    hdf5->file = H5Fopen(filename, H5F_ACC_RDONLY, plist);
    hdf5->group_crd = H5Gopen(hdf5->file, "crd", H5P_DEFAULT);
  } else {
    assert(0);
  }
  H5Pclose(plist);

  for (int d = 0; d < 3; d++) {
    hdf5->crd_written[d] = false;
  }

  INIT_LIST_HEAD(&hdf5->xdmf_spatial_list);
}

static void
ds_xdmf_parallel_close(struct mrc_io *io)
{
#ifndef NDEBUG
  MPI_Barrier(io->obj.comm);
#endif
  
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

  H5Gclose(hdf5->group_crd);
  if (strcmp(hdf5->mode, "w") == 0) {
    hdf5_close(io);
  } else {
    H5Fclose(hdf5->file);
  }

  if (io->rank != 0)
    return;

  xdmf_close(io);
  free(hdf5->mode);
}

static void
hdf5_write_crds_parallel(struct mrc_io *io, struct mrc_f3 *fld)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

  const char *xyz[3] = { "x", "y", "z" };

  int gdims[3], nr_procs[3], idx[3];
  mrc_domain_get_global_dims(fld->domain, gdims);
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(fld->domain, &nr_patches);
  assert(nr_patches == 1);
  int *ldims = patches[0].ldims, *off = patches[0].off;
  mrc_domain_get_local_idx(fld->domain, idx);
  mrc_domain_get_nr_procs(fld->domain, nr_procs);
  for (int d = 0; d < 3; d++) {
    bool skip_write = false;
    for (int dd = 0; dd < 3; dd++) {
      if (dd != d && idx[dd] != 0) {
	skip_write = true;
	break;
      }
    }
    hsize_t hgdims[1] = { gdims[d] + 1 };
    hsize_t hldims[1] = { ldims[d] };
    hsize_t hoff[1] = { off[d] };
    if (idx[d] == nr_procs[d] - 1) {
      hldims[0]++;
    }
    hid_t filespace = H5Screate_simple(1, hgdims, NULL);
    hid_t memspace = H5Screate_simple(1, hldims, NULL);
    hid_t dset = H5Dcreate(hdf5->group_crd, xyz[d], H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
			   H5P_DEFAULT, H5P_DEFAULT);
    if (!skip_write) {
      // FIXME make sure it is indep I/O
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, hoff, NULL, hldims, NULL);
      struct mrc_crds *crds = mrc_domain_get_crds(fld->domain);
      float *crd = &MRC_CRD(crds, d, fld->sw);
      float *crd_nc = calloc(ldims[d] + 1, sizeof(*crd_nc));

      if (fld->sw > 0) { // FIXME, shouldn't be done here, and should
	                 // be based on sw for the crds
	for (int i = 0; i <= ldims[d]; i++) {
	  crd_nc[i] = .5 * (crd[i-1] + crd[i]);
	}
      } else {
	if (ldims[d] > 1) {
	  for (int i = 1; i < ldims[d]; i++) {
	    crd_nc[i] = .5 * (crd[i-1] + crd[i]);
	  }
	  // extrapolate
	  crd_nc[0]  = crd[0] - .5 * (crd[1] - crd[0]);
	  crd_nc[ldims[d]] = crd[ldims[d]-1] + .5 * (crd[ldims[d]-1] - crd[ldims[d]-2]);
	} else { // this is an invariant direction...
	  crd_nc[0] = crd[0] - 5e-3;
	  crd_nc[1] = crd[0] + 5e-3;
	}
      }
      H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, crd_nc);
    }
    H5Dclose(dset);
    H5Sclose(filespace);
    H5Sclose(memspace);
  }
}

struct read_f1_cb_data {
  struct mrc_io *io;
  struct mrc_f1 *fld;
  struct mrc_f1 *lfld;
  hid_t filespace;
  hid_t memspace;
  hid_t dxpl;
};

static herr_t
read_f1_cb(hid_t g_id, const char *name, const H5L_info_t *info, void *op_data)
{
  struct read_f1_cb_data *data = op_data;

  hid_t group = H5Gopen(g_id, name, H5P_DEFAULT);
  int m;
  H5LTget_attribute_int(group, ".", "m", &m);

  hid_t dset = H5Dopen(group, "1d", H5P_DEFAULT);
  H5Dread(dset, H5T_NATIVE_FLOAT, data->memspace, data->filespace, data->dxpl, data->fld->arr);
  H5Dclose(dset);

  H5Gclose(group);

  return 0;
}

static void
ds_xdmf_parallel_read_f1(struct mrc_io *io, const char *path, struct mrc_f1 *f1)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

  assert(io->size == 1);
#ifndef NDEBUG
  MPI_Barrier(io->obj.comm);
#endif

  hid_t group0 = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  hsize_t hdims[1] = { f1->im[0] };

  hid_t filespace = H5Screate_simple(1, hdims, NULL);
  hid_t memspace = H5Screate_simple(1, hdims, NULL);
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  if (hdf5->par.use_independent_io) {
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
  } else {
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
  }

  struct read_f1_cb_data cb_data = {
    .io        = io,
    .fld       = f1,
    .filespace = filespace,
    .memspace  = memspace,
    .dxpl      = dxpl,
  };
  hsize_t idx = 0;
  H5Literate_by_name(group0, ".", H5_INDEX_NAME, H5_ITER_INC, &idx,
		     read_f1_cb, &cb_data, H5P_DEFAULT);

  H5Pclose(dxpl);
  H5Sclose(filespace);
  H5Sclose(memspace);

  H5Gclose(group0);
}

struct read_f3_cb_data {
  struct mrc_io *io;
  struct mrc_f3 *fld;
  struct mrc_f3 *lfld;
  hid_t filespace;
  hid_t memspace;
  hid_t dxpl;
};

static herr_t
read_f3_cb(hid_t g_id, const char *name, const H5L_info_t *info, void *op_data)
{
  struct read_f3_cb_data *data = op_data;

  hid_t group = H5Gopen(g_id, name, H5P_DEFAULT);
  int m;
  H5LTget_attribute_int(group, ".", "m", &m);

  hid_t dset = H5Dopen(group, "3d", H5P_DEFAULT);
  H5Dread(dset, H5T_NATIVE_FLOAT, data->memspace, data->filespace, data->dxpl, data->lfld->arr);
  H5Dclose(dset);

  // FIXME, could be done w/hyperslab, vectors...
  copy_back(data->lfld, 0, data->fld, m, data->fld->sw);

  H5Gclose(group);

  return 0;
}

static void
ds_xdmf_parallel_read_f3(struct mrc_io *io, const char *path, struct mrc_f3 *fld)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

#ifndef NDEBUG
  MPI_Barrier(io->obj.comm);
#endif

  int gdims[3];
  mrc_domain_get_global_dims(fld->domain, gdims);
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(fld->domain, &nr_patches);
  assert(nr_patches == 1);
  int *off = patches[0].off, *ldims = patches[0].ldims;

  struct mrc_f3 *lfld = mrc_f3_alloc(MPI_COMM_SELF, NULL, ldims);
  mrc_f3_setup(lfld);

  hid_t group0 = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  hsize_t hgdims[3] = { gdims[2], gdims[1], gdims[0] };
  hsize_t hldims[3] = { ldims[2], ldims[1], ldims[0] };
  hsize_t hoff[3]   = { off[2], off[1], off[0] };

  hid_t filespace = H5Screate_simple(3, hgdims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, hoff, NULL, hldims, NULL);
  hid_t memspace = H5Screate_simple(3, hldims, NULL);
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  if (hdf5->par.use_independent_io) {
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
  } else {
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
  }

  struct read_f3_cb_data cb_data = {
    .io        = io,
    .fld       = fld,
    .lfld      = lfld,
    .filespace = filespace,
    .memspace  = memspace,
    .dxpl      = dxpl,
  };
  hsize_t idx = 0;
  H5Literate_by_name(group0, ".", H5_INDEX_NAME, H5_ITER_INC, &idx,
		     read_f3_cb, &cb_data, H5P_DEFAULT);

  H5Pclose(dxpl);
  H5Sclose(filespace);
  H5Sclose(memspace);

  H5Gclose(group0);

  mrc_f3_destroy(lfld);
}

static void
ds_xdmf_parallel_write_f1(struct mrc_io *io, const char *path, struct mrc_f1 *f1)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

#ifndef NDEBUG
  MPI_Barrier(io->obj.comm);
#endif

  if (io->size > 1)
    return; // FIXME

  assert(io->size == 1);
  hid_t group0 = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  hsize_t hdims[1] = { f1->im[0] };
  for (int m = 0; m < f1->nr_comp; m++) {
    hid_t group = H5Gcreate(group0, f1->name[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5LTset_attribute_int(group, ".", "m", &m, 1);
    
    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
    if (hdf5->par.use_independent_io) {
      H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
    } else {
      H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
    }
    hid_t filespace = H5Screate_simple(1, hdims, NULL);
    hid_t memspace = H5Screate_simple(1, hdims, NULL);
    hid_t dset = H5Dcreate(group, "1d", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
			   H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, dxpl, f1->arr);
    H5Dclose(dset);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(dxpl);

    H5Gclose(group);
  }
  H5Gclose(group0);
}

static void
ds_xdmf_parallel_write_field(struct mrc_io *io, const char *path,
			     float scale, struct mrc_f3 *fld, int m)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(io);

#ifndef NDEBUG
  MPI_Barrier(io->obj.comm);
#endif

  int gdims[3];
  mrc_domain_get_global_dims(fld->domain, gdims);
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(fld->domain, &nr_patches);
  assert(nr_patches == 1);
  int *off = patches[0].off, *ldims = patches[0].ldims;

  // strip boundary, could be done through hyperslab, but
  // still have to scale, anyway
  struct mrc_f3 *lfld = mrc_f3_alloc(MPI_COMM_SELF, NULL, ldims);
  mrc_f3_setup(lfld);
  copy_and_scale(lfld, 0, fld, m, scale, fld->sw);

  if (!hdf5->crds_done) {
    hdf5_write_crds_parallel(io, fld);
    hdf5->crds_done = true;
  }

  if (io->rank == 0) {
    struct xdmf_spatial *xs = xdmf_spatial_find(io, "3df", -1);
    if (!xs) {
      xs = xdmf_spatial_create_3d(io, gdims, -1, 1);
    }
    save_fld_info(xs, strdup(fld->name[m]), strdup(path), false);
  }

  hid_t group0 = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  hsize_t hgdims[3] = { gdims[2], gdims[1], gdims[0] };
  hsize_t hldims[3] = { ldims[2], ldims[1], ldims[0] };
  hsize_t hoff[3]   = { off[2], off[1], off[0] };
  hid_t group = H5Gcreate(group0, fld->name[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_int(group, ".", "m", &m, 1);

  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  if (hdf5->par.use_independent_io) {
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
  } else {
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
  }
  hid_t filespace = H5Screate_simple(3, hgdims, NULL);
  hid_t memspace = H5Screate_simple(3, hldims, NULL);
  hid_t dset = H5Dcreate(group, "3d", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
			 H5P_DEFAULT, H5P_DEFAULT);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, hoff, NULL, hldims, NULL);
  H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, dxpl, lfld->arr);
  H5Dclose(dset);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(dxpl);

  H5Gclose(group);
  H5Gclose(group0);

  mrc_f3_destroy(lfld);
}

struct mrc_io_ops mrc_io_xdmf_parallel_ops = {
  .name          = "xdmf_parallel",
  .size          = sizeof(struct diag_hdf5),
  .param_descr   = diag_hdf5_params_descr,
  .param_offset  = offsetof(struct diag_hdf5, par),
  .parallel      = true,
  .destroy       = ds_xdmf_destroy,
  .setup         = ds_xdmf_setup,
  .open          = ds_xdmf_parallel_open,
  .close         = ds_xdmf_parallel_close,
  .read_f1       = ds_xdmf_parallel_read_f1,
  .read_f3       = ds_xdmf_parallel_read_f3,
  .write_f1      = ds_xdmf_parallel_write_f1,
  .write_field   = ds_xdmf_parallel_write_field,
  .read_attr     = ds_xdmf_read_attr,
  .write_attr    = ds_xdmf_write_attr,
};

#endif

// ======================================================================

#endif

