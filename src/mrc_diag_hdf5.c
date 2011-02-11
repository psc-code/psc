
#include "mrc_diag_private.h"
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
#define MAX_TIMESTEPS (100000)

struct fld_info {
  char *name;
  bool is_vec;
  int dim;
};

struct xdmf_temporal {
  char *sfx;
  int nr_timestep;
  int *timestep;
  list_t list;
};

struct xdmf_subdomain {
  int im[3];
  int fileno;
};

struct xdmf_spatial {
  char *sfx;
  float sheet;
  int nr_subdomains;
  int size;
  struct xdmf_subdomain *subdomains;
  int nr_fld_info;
  struct fld_info fld_info[MAX_FLD_INFO];
  void (*write_topology)(FILE *f, int im[3], const char *filename, float sheet);
  void (*write_fld)(FILE *f, struct fld_info *fld_info, int im[3], const char *filename,
		    const char *path);
  list_t list;
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
  hid_t file;
  int step;
  bool crd_written[3]; // this should be per mrc_domain, but this will do for crds
  bool crds_done;
  int gdims[3], nr_procs[3];
  int off[3], ldims[3];
  struct mrc_f3 vfld;
  list_t xdmf_temporal_list; // lives as long as diag_hdf5
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

#define diag_hdf5(format) ((struct diag_hdf5 *)format->obj.subctx)
#define diag_format(obj) (container_of(obj, struct diag_format, obj))

static void
hdf5_open(struct diag_format *format, float sheet, int outtype, int step) 
{
  struct diag_hdf5 *hdf5 = diag_hdf5(format);
  hdf5->step = step;

  char *filename = diagc_make_filename(format, ".h5", sheet, outtype, step);
  hdf5->file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  free(filename);
}

static void
hdf5_close(struct diag_format *format) 
{
  herr_t ierr;
  struct diag_hdf5 *hdf5 = diag_hdf5(format);
  struct diag_info *diag_info = &format->diag_info;

  hid_t group = H5Gcreate(hdf5->file, "openggcm", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_string(group, ".", "run", diag_info->run);
  H5LTset_attribute_int(group, ".", "step", &hdf5->step, 1);
  H5LTset_attribute_string(group, ".", "time_str", format->time_str);
  H5LTset_attribute_float(group, ".", "time", &format->time, 1);

  H5LTset_attribute_int(group, ".", "global_dims", hdf5->gdims, 3);
  
  int size;
  MPI_Comm_size(format->obj.comm, &size);
  int parallel_write = (size > 1);
  H5LTset_attribute_int(group, ".", "parallel_write", &parallel_write, 1);
  H5LTset_attribute_int(group, ".", "proc_rank", &format->rank, 1);
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
hdf5_write_field2d_serial(struct diag_format *format, float scale, struct mrc_f2 *fld,
			  const char *fld_name, const char *path)
{
  herr_t ierr;
  struct diag_hdf5 *hdf5 = diag_hdf5(format);

  hsize_t fdims[2] = { fld->im[1], fld->im[0] };
  //  printf("[%d] diagsrv: write '%s'\n", info->rank, fld_name);

  hid_t group0 = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  hid_t group = H5Gcreate(group0, fld_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
xdmf_write_topology_3d(FILE *f, int im[3], const char *filename, float sheet)
{
  fprintf(f, "     <Topology TopologyType=\"3DRectMesh\" Dimensions=\"%d %d %d\"/>\n",
	  im[2] + 1, im[1] + 1, im[0] + 1);
  fprintf(f, "     <Geometry GeometryType=\"VXVYVZ\">\n");
  fprintf(f, "     <DataItem Name=\"VX\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[0] + 1);
  fprintf(f, "        ./%s:/crd/x\n", filename);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VY\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[1] + 1);
  fprintf(f, "        ./%s:/crd/y\n", filename);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VZ\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[2] + 1);
  fprintf(f, "        ./%s:/crd/z\n", filename);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     </Geometry>\n");
  fprintf(f, "\n");
}

static void
xdmf_write_topology_2d_x(FILE *f, int im[2], const char *filename, float sheet)
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
xdmf_write_topology_2d_y(FILE *f, int im[2], const char *filename, float sheet)
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
xdmf_write_topology_2d_z(FILE *f, int im[2], const char *filename, float sheet)
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
xdmf_write_topology_iono(FILE *f, int im[2], const char *filename, float sheet)
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
xdmf_write_fld_2d_iono(FILE *f, struct fld_info *fld_info, int im[2], const char *filename,
		       const char *path)
{
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",
	  fld_info->name);
  fprintf(f, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
	  im[1], im[0]);
  fprintf(f, "        %s:/%s/%s/2d\n", filename, path, fld_info->name);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
}

static void
xdmf_write_fld_2d_x(FILE *f, struct fld_info *fld_info, int im[2], const char *filename,
		    const char *path)
{
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
	  fld_info->name);
  fprintf(f, "       <DataItem Dimensions=\"1 %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
	  im[1], im[0]);
  fprintf(f, "        %s:/%s/%s/2d\n", filename, path, fld_info->name);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
}

static void
xdmf_write_fld_2d_y(FILE *f, struct fld_info *fld_info, int im[2], const char *filename,
		    const char *path)
{
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
	  fld_info->name);
  fprintf(f, "       <DataItem Dimensions=\"%d 1 %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
	  im[1], im[0]);
  fprintf(f, "        %s:/%s/%s/2d\n", filename, path, fld_info->name);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
}

static void
xdmf_write_fld_2d_z(FILE *f, struct fld_info *fld_info, int im[2], const char *filename,
		    const char *path)
{
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
	  fld_info->name);
  fprintf(f, "       <DataItem Dimensions=\"%d %d 1\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
	  im[1], im[0]);
  fprintf(f, "        %s:/%s/%s/2d\n", filename, path, fld_info->name);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
}

static void
xdmf_write_fld_3d(FILE *f, struct fld_info *fld_info, int im[3], const char *filename,
		  const char *path)
{
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Cell\">\n",
	  fld_info->name, fld_info->is_vec ? "Vector" : "Scalar");
  fprintf(f, "       <DataItem Dimensions=\"%d %d %d%s\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
	  im[2], im[1], im[0], fld_info->is_vec ? " 3" : "");
  fprintf(f, "        ./%s:/%s/%s/3d\n", filename, path, fld_info->name);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
}

static struct xdmf_spatial *
xdmf_spatial_create(struct diag_format *format, const char *sfx, float sheet, int nr_subdomains,
		    struct xdmf_subdomain *subdomains, int size,
		    void (*write_topology)(FILE *f, int im[3], const char *filename, float sheet),
		    void (*write_fld)(FILE *f, struct fld_info *fld_info, int im[3],
				      const char *filename, const char *path))
{
  struct diag_hdf5 *hdf5 = diag_hdf5(format);
  struct xdmf_spatial *xs = malloc(sizeof(*xs));
  memset(xs, 0, sizeof(*xs));
  xs->sfx            = strdup(sfx);
  xs->sheet          = sheet;
  xs->nr_subdomains  = nr_subdomains;
  xs->subdomains     = subdomains;
  xs->size           = size;
  xs->write_topology = write_topology;
  xs->write_fld      = write_fld;
  
  list_add_tail(&xs->list, &hdf5->xdmf_spatial_list);

  return xs;
}

static struct xdmf_spatial *
xdmf_spatial_find(struct diag_format *format, const char *sfx, float sheet)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(format);
  struct xdmf_spatial *xs;
  list_for_each_entry(xs, &hdf5->xdmf_spatial_list, list) {
    if (xs->sheet == sheet && strcmp(xs->sfx, sfx) == 0) {
      return xs;
    }
  }
  return NULL;
}

static void
xdmf_spatial_destroy(struct xdmf_spatial *xs)
{
  list_del(&xs->list);

  for (int i = 0; i < xs->nr_fld_info; i++) {
    free(xs->fld_info[i].name);
  }
  free(xs->subdomains);
  free(xs->sfx);
  free(xs);
}

static void
xdmf_write_spatial_collection(struct diag_format *format, int step, struct xdmf_spatial *xs)
{
  if (format->rank > 0)
    return;

  struct diag_info *diag_info = &format->diag_info;

  char filename[strlen(diag_info->outdir) + strlen(diag_info->run) + 30];
  sprintf(filename, "%s/%s.%s.%06d.xdmf", diag_info->outdir, diag_info->run, xs->sfx, step);
  mprintf("diag: writing '%s'\n", filename);

  FILE *f = fopen(filename, "w");
  xdmf_write_header(f);
  fprintf(f, "<Domain>\n");
  fprintf(f, "<Grid GridType=\"Collection\" CollectionType=\"Spatial\">\n");
  fprintf(f, "   <Time Type=\"Single\" Value=\"%g\" />\n", format->time);
  for (int s = 0; s < xs->nr_subdomains; s++) {
    fprintf(f, "   <Grid Name=\"mesh-%s-%d-%d\" GridType=\"Uniform\">\n", xs->sfx, s, step);
    
    char filename[strlen(diag_info->run) + 30];
    if (xs->size == 1) {
      sprintf(filename, "%s.%s.%06d.h5", diag_info->run, xs->sfx, step);
    } else {
      sprintf(filename, "%s.%s.%06d_p%06d.h5", diag_info->run, xs->sfx, step, xs->subdomains[s].fileno);
    }
    
    int *im = xs->subdomains[s].im;
    xs->write_topology(f, im, filename, xs->sheet);
    
    char path[100];
    make_path(path, xs->sfx);
    for (int m = 0; m < xs->nr_fld_info; m++) {
      xs->write_fld(f, &xs->fld_info[m], im, filename, path);
    }
    fprintf(f, "   </Grid>\n");
  }
  fprintf(f, "</Grid>\n");
  fprintf(f, "</Domain>\n");
  fprintf(f, "</Xdmf>\n");
  fclose(f);
}

static void
xdmf_write_temporal_collection(struct diag_format *format, const char *sfx, int step)
{
  struct diag_info *diag_info = &format->diag_info;
  struct diag_hdf5 *hdf5 = diag_hdf5(format);

  if (format->rank > 0) {
    return;
  }

  struct xdmf_temporal *xt;
  list_for_each_entry(xt, &hdf5->xdmf_temporal_list, list) {
    if (strcmp(xt->sfx, sfx) == 0)
      goto found;
  }
  xt = malloc(sizeof(*xt));
  memset(xt, 0, sizeof(*xt));
  xt->timestep = calloc(sizeof(*xt->timestep), MAX_TIMESTEPS);
  xt->sfx = strdup(sfx);
  list_add_tail(&xt->list, &hdf5->xdmf_temporal_list);

 found:
  assert(xt->nr_timestep < MAX_TIMESTEPS);
  xt->timestep[xt->nr_timestep++] = step;

  // It'd be easier to create those files line by line as time goes by.
  // However, then we won't be able to get the timeseries into Paraview
  // until the solution is all finished.
  // So this version rewrites the xdmf file completely every timestep,
  // which needs however some way to figure out what times we've written
  // before.

  char fname[strlen(diag_info->outdir) + strlen(diag_info->run) + 20];
  sprintf(fname, "%s/%s.%s.xdmf", diag_info->outdir, diag_info->run, sfx);
  FILE *f = fopen(fname, "w");

  xdmf_write_header(f);
  fprintf(f, "<Domain>\n");
  fprintf(f, "  <Grid GridType='Collection' CollectionType='Temporal'>\n");
  for (int i = 0; i < xt->nr_timestep; i++) {
    fprintf(f, "  <xi:include href='%s.%s.%06d.xdmf' xpointer='xpointer(//Xdmf/Domain/Grid)'/>\n", diag_info->run, sfx, xt->timestep[i]);
  }
  fprintf(f, "  </Grid>\n");
  fprintf(f, "  </Domain>\n");
  fprintf(f, "</Xdmf>\n");
  fclose(f);
}

static void
ds_xdmf_create(struct mrc_obj *obj)
{
  struct diag_format *format = diag_format(obj);
  struct diag_hdf5 *hdf5 = diag_hdf5(format);
  INIT_LIST_HEAD(&hdf5->xdmf_temporal_list);
}

static void
ds_xdmf_destroy(struct mrc_obj *obj)
{
  struct diag_format *format = (struct diag_format *)obj;
  struct diag_hdf5 *hdf5 = diag_hdf5(format);

  while (!list_empty(&hdf5->xdmf_temporal_list)) {
    struct xdmf_temporal *xt = list_entry(hdf5->xdmf_temporal_list.next, typeof(*xt), list);
    list_del(&xt->list);
    free(xt->sfx);
    free(xt->timestep);
    free(xt);
  }
}

static void
ds_xdmf_open(struct diag_format *format, float sheet, int outtype, int step) 
{
  struct diag_hdf5 *hdf5 = diag_hdf5(format);
  hdf5_open(format, sheet, outtype, step);
  hdf5->group_crd = H5Gcreate(hdf5->file, "crd", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for (int d = 0; d < 3; d++) {
    hdf5->crd_written[d] = false;
  }
  hdf5->crds_done = false;

  INIT_LIST_HEAD(&hdf5->xdmf_spatial_list);
}

static void
hdf5_write_crds(struct diag_format *format, int im[3], struct mrc_domain *domain, int sw)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(format);

  if (!hdf5->gdims[0]) {
    mrc_domain_get_global_dims(domain, hdf5->gdims);
    mrc_domain_get_nr_procs(domain, hdf5->nr_procs);
    mrc_domain_get_local_offset_dims(domain, hdf5->off, hdf5->ldims);
  }

  const char *xyz[3] = { "x", "y", "z" };
  
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  for (int d = 0; d < 3; d++) {
    if (im[d] == 0 || hdf5->crd_written[d])
      continue;

    float *crd = crds->crd[d] + sw;
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
ds_xdmf_close(struct diag_format *format)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(format);

  H5Gclose(hdf5->group_crd);
  hdf5_close(format);

  while (!list_empty(&hdf5->xdmf_spatial_list)) {
    struct xdmf_spatial *xs = list_entry(hdf5->xdmf_spatial_list.next, typeof(*xs), list);

    xdmf_write_spatial_collection(format, hdf5->step, xs);
    xdmf_write_temporal_collection(format, xs->sfx, hdf5->step);
    xdmf_spatial_destroy(xs);
  }
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
save_fld_info(struct xdmf_spatial *xs, char *fld_name, bool is_vec)
{
  assert(xs->nr_fld_info < MAX_FLD_INFO);
  struct fld_info *fld_info = &xs->fld_info[xs->nr_fld_info++];

  fld_info->name = fld_name;
  fld_info->is_vec = is_vec;
}

static struct xdmf_spatial *
xdmf_spatial_create_3d(struct diag_format *format, int im[3], int nr_subdomains)
{
  // OPT, we could skip this on procs which aren't writing xdmf

  int size_format;
  MPI_Comm_size(format->obj.comm, &size_format);

  struct xdmf_subdomain *subdomains = calloc(nr_subdomains, sizeof(*subdomains));
  for (int s = 0; s < nr_subdomains; s++) {
    for (int d = 0; d < 3; d++) {
      subdomains[s].im[d] = im[d];
    }
    subdomains[s].fileno = s;
  }

  return xdmf_spatial_create(format, "3df", -1, nr_subdomains, subdomains, size_format,
			     xdmf_write_topology_3d, xdmf_write_fld_3d);
}

static struct xdmf_spatial *
xdmf_spatial_create_2d_x(struct diag_format *format, int im[2], 
			 const char *sfx, int sheet, int nr_subdomains)
{
  // OPT, we could skip this on procs which aren't writing xdmf

  int size_format;
  MPI_Comm_size(format->obj.comm, &size_format);

  struct xdmf_subdomain *subdomains = calloc(nr_subdomains, sizeof(*subdomains));
  for (int s = 0; s < nr_subdomains; s++) {
    for (int d = 0; d < 2; d++) {
      subdomains[s].im[d] = im[d];
    }
    subdomains[s].fileno = s;
  }

  return xdmf_spatial_create(format, sfx, sheet, nr_subdomains, subdomains, size_format,
			     xdmf_write_topology_2d_x, xdmf_write_fld_2d_x);
}

static struct xdmf_spatial *
xdmf_spatial_create_2d_y(struct diag_format *format, int im[2], 
			 const char *sfx, int sheet, int nr_subdomains)
{
  // OPT, we could skip this on procs which aren't writing xdmf

  int size_format;
  MPI_Comm_size(format->obj.comm, &size_format);

  struct xdmf_subdomain *subdomains = calloc(nr_subdomains, sizeof(*subdomains));
  for (int s = 0; s < nr_subdomains; s++) {
    for (int d = 0; d < 2; d++) {
      subdomains[s].im[d] = im[d];
    }
    subdomains[s].fileno = s;
  }

  return xdmf_spatial_create(format, sfx, sheet, nr_subdomains, subdomains, size_format,
			     xdmf_write_topology_2d_y, xdmf_write_fld_2d_y);
}

static struct xdmf_spatial *
xdmf_spatial_create_2d_z(struct diag_format *format, int im[2], 
			 const char *sfx, int sheet, int nr_subdomains)
{
  // OPT, we could skip this on procs which aren't writing xdmf

  int size_format;
  MPI_Comm_size(format->obj.comm, &size_format);

  struct xdmf_subdomain *subdomains = calloc(nr_subdomains, sizeof(*subdomains));
  for (int s = 0; s < nr_subdomains; s++) {
    for (int d = 0; d < 2; d++) {
      subdomains[s].im[d] = im[d];
    }
    subdomains[s].fileno = s;
  }

  return xdmf_spatial_create(format, sfx, sheet, nr_subdomains, subdomains, size_format,
			     xdmf_write_topology_2d_z, xdmf_write_fld_2d_z);
}

static struct xdmf_spatial *
xdmf_spatial_create_iono(struct diag_format *format, int im[2])
{
  struct xdmf_subdomain *subdomain = calloc(1, sizeof(*subdomain));
  for (int d = 0; d < 2; d++) {
    subdomain->im[d] = im[d];
  }
  subdomain->fileno = -1;

  return xdmf_spatial_create(format, "iof", -1, 1, subdomain, 1,
			     xdmf_write_topology_iono, xdmf_write_fld_2d_iono);
}

static void
ds_xdmf_write_field(struct diag_format *format, float scale, struct mrc_f3 *fld,
		    int m)
{
  herr_t ierr;
  struct diag_hdf5 *hdf5 = diag_hdf5(format);

  // on diagsrv, ghost points have already been dropped
  int sw = fld->sw;
  int im[3] = { fld->im[0] - 2 * sw, fld->im[1] - 2 * sw, fld->im[2] - 2 * sw };

  char path[100]; make_path(path, "3df");

  struct xdmf_spatial *xs;
  xs = xdmf_spatial_find(format, "3df", -1);
  if (!xs) {
    int size; MPI_Comm_size(format->obj.comm, &size);
    xs = xdmf_spatial_create_3d(format, im, size);
    hid_t group0 = H5Gcreate(hdf5->file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group0);
    hdf5_write_crds(format, im, fld->domain, fld->sw);
  }

  hid_t group0 = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  struct mrc_f3 *vfld = &hdf5->vfld;
  char c = fld->name[m][strlen(fld->name[m]) - 1];
  if (c == 'x') {
    assert(!vfld->arr);
    mrc_f3_alloc(vfld, NULL, im, 3);
    copy_and_scale(vfld, 0, fld, m, scale, sw);
  } else if (c == 'y') {
    copy_and_scale(vfld, 1, fld, m, scale, sw);
  } else if (c == 'z') {
    copy_and_scale(vfld, 2, fld, m, scale, sw);
    char *vec_name = strdup(fld->name[m]);
    vec_name[strlen(fld->name[m])-1] = 0;
    save_fld_info(xs, vec_name, true);
    hsize_t hdims[4] = { vfld->im[2], vfld->im[1], vfld->im[0], vfld->nr_comp };
    hid_t group = H5Gcreate(group0, vec_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    ierr = H5LTmake_dataset_float(group, "3d", 4, hdims, vfld->arr); CE;
    ierr = H5Gclose(group); CE;

    mrc_f3_free(vfld);
  } else { // scalar
    save_fld_info(xs, strdup(fld->name[m]), false);
    if (format->rank < 0) { // FIXME not happening anymore, still optimize this case
      // on diag srv, ghost points are already gone and scaling is done
      assert(scale == 1.f);
      hsize_t hdims[3] = { fld->im[2], fld->im[1], fld->im[0] };
      hid_t group = H5Gcreate(group0, fld->name[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      ierr = H5LTmake_dataset_float(group, "3d", 3, hdims, &MRC_F3(fld, m, 0,0,0)); CE;
      ierr = H5Gclose(group); CE;
    } else {
      struct mrc_f3 sfld;
      mrc_f3_alloc(&sfld, NULL, im, 1);
      copy_and_scale(&sfld, 0, fld, m, scale, sw);

      hsize_t hdims[3] = { sfld.im[2], sfld.im[1], sfld.im[0] };
      hid_t group = H5Gcreate(group0, fld->name[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      ierr = H5LTmake_dataset_float(group, "3d", 3, hdims, sfld.arr); CE;
      ierr = H5Gclose(group); CE;
      mrc_f3_free(&sfld);
    }
  }
  H5Gclose(group0);
}

static void
ds_xdmf_write_field2d(struct diag_format *format, float scale, struct mrc_f2 *fld,
		      const char *fld_name, int outtype, float sheet)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(format);
  int ierr;

  // diagsrv xdmf_serial / single proc only for now
  int bnd = 0;
  int im[2] = { fld->im[0] - bnd, fld->im[1] - bnd };

  struct xdmf_spatial *xs = NULL;
  char sfx[10]; make_sfx(sfx, outtype, sheet);
  char path[100]; make_path(path, sfx);
  xs = xdmf_spatial_find(format, sfx, sheet);
  if (!xs) {
    int size; MPI_Comm_size(format->obj.comm, &size);
    if (outtype == DIAG_TYPE_2D_X) {
      xs = xdmf_spatial_create_2d_x(format, im, sfx, sheet, size);
      hdf5_write_crds(format, (int[]) { 0, im[0], im[1] }, fld->domain, fld->sw);
    } else if (outtype == DIAG_TYPE_2D_Y) {
      xs = xdmf_spatial_create_2d_y(format, im, sfx, sheet, size);
      hdf5_write_crds(format, (int[]) { im[0], 0, im[1] }, fld->domain, fld->sw);
    } else if (outtype == DIAG_TYPE_2D_Z) {
      xs = xdmf_spatial_create_2d_z(format, im, sfx, sheet, size);
      hdf5_write_crds(format, (int[]) { im[0], im[1], 0 }, fld->domain, fld->sw);
    } else if (outtype == DIAG_TYPE_2D_IONO) {
      xs = xdmf_spatial_create_iono(format, im);
    }
    hid_t group = H5Gcreate(hdf5->file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    ierr = H5Gclose(group); CE;
  }

  save_fld_info(xs, strdup(fld_name), false);
  hdf5_write_field2d_serial(format, scale, fld, fld_name, path);
}

static struct diag_format_ops ds_xdmf_ops = {
  .name          = "xdmf",
  .size          = sizeof(struct diag_hdf5),
  .parallel      = true,
  .create        = ds_xdmf_create,
  .destroy       = ds_xdmf_destroy,
  .open          = ds_xdmf_open,
  .close         = ds_xdmf_close,
  .write_field   = ds_xdmf_write_field,
  .write_field2d = ds_xdmf_write_field2d,
};

static struct diag_format_ops ds_xdmf_serial_ops = {
  .name          = "xdmf_serial",
  .size          = sizeof(struct diag_hdf5),
  .create        = ds_xdmf_create,
  .destroy       = ds_xdmf_destroy,
  .open          = ds_xdmf_open,
  .close         = ds_xdmf_close,
  .write_field   = ds_xdmf_write_field,
  .write_field2d = ds_xdmf_write_field2d,
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
ds_xdmf_to_one_open(struct diag_format *format, float sheet, int outtype, int step) 
{
#ifndef NDEBUG
  MPI_Barrier(format->obj.comm);
#endif
  
  if (format->rank != 0)
    return;

  struct diag_hdf5 *hdf5 = diag_hdf5(format);
  hdf5_open(format, sheet, outtype, step);
  hdf5->group_crd = H5Gcreate(hdf5->file, "crd", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for (int d = 0; d < 3; d++) {
    hdf5->crd_written[d] = false;
  }

  INIT_LIST_HEAD(&hdf5->xdmf_spatial_list);
}

static void
ds_xdmf_to_one_close(struct diag_format *format)
{
#ifndef NDEBUG
  MPI_Barrier(format->obj.comm);
#endif
  
  if (format->rank != 0)
    return;

  struct diag_hdf5 *hdf5 = diag_hdf5(format);

  H5Gclose(hdf5->group_crd);
  hdf5_close(format);

  while (!list_empty(&hdf5->xdmf_spatial_list)) {
    struct xdmf_spatial *xs = list_entry(hdf5->xdmf_spatial_list.next, typeof(*xs), list);

    xdmf_write_spatial_collection(format, hdf5->step, xs);
    xdmf_write_temporal_collection(format, xs->sfx, hdf5->step);
    xdmf_spatial_destroy(xs);
  }
}

static void
communicate_crds(struct diag_format *format, struct mrc_f3 *gfld, struct mrc_f3 *lfld)
{
  int sw = gfld->sw;

  struct mrc_domain *gdomain = gfld->domain;
  struct mrc_crds *gcrds = mrc_domain_get_crds(gdomain);
  if (format->rank != 0) {
    int iw[6], *ib = iw, *im = iw + 3;
    mrc_domain_get_local_offset_dims(gdomain, ib, im);
    MPI_Send(iw, 6, MPI_INT, 0, TAG_OFF_DIMS, format->obj.comm);
    for (int d = 0; d < 3; d++) {
      MPI_Send(gcrds->crd[d] + sw, im[d], MPI_FLOAT, 0,
	       TAG_CRDX + d, format->obj.comm);
    }
  } else { // format->rank == 0
    struct mrc_domain *ldomain = lfld->domain;
    struct mrc_crds *lcrds = mrc_domain_get_crds(ldomain);
    int size;
    MPI_Comm_size(format->obj.comm, &size);
    for (int n = 0; n < size; n++) {
      float *recv_crds[3];
      int iw[6], *ib = iw, *im = iw + 3;
      if (n == 0) {
	mrc_domain_get_local_offset_dims(gdomain, ib, im);
	for (int d = 0; d < 3; d++) {
	  recv_crds[d] = gcrds->crd[d] + sw;
	}
      } else {
	MPI_Recv(iw, 6, MPI_INT, n, TAG_OFF_DIMS, format->obj.comm, MPI_STATUS_IGNORE);
	for (int d = 0; d < 3; d++) {
	  recv_crds[d] = calloc(im[d], sizeof(*recv_crds[d]));
	  MPI_Recv(recv_crds[d], im[d], MPI_FLOAT, n, TAG_CRDX + d, format->obj.comm,
		   MPI_STATUS_IGNORE);
	}
      }

      for (int d = 0; d < 3; d++) {
	for (int i = 0; i < im[d]; i++) {
	  lcrds->crd[d][i + ib[d]] = recv_crds[d][i];
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
communicate_fld(struct diag_format *format, struct mrc_f3 *gfld, int m, float scale,
		struct mrc_f3 *lfld)
{
  int sw = gfld->sw;

  struct mrc_f3 send_fld;
  mrc_domain_f3_alloc(gfld->domain, &send_fld, 1, SW_0);
  copy_and_scale(&send_fld, 0, gfld, m, scale, sw);

  if (format->rank != 0) {
    int iw[6], *ib = iw, *im = iw + 3;
    mrc_domain_get_local_offset_dims(send_fld.domain, ib, im);
    MPI_Send(iw, 6, MPI_INT, 0, TAG_OFF_DIMS, format->obj.comm);
    MPI_Send(send_fld.arr, send_fld.len, MPI_FLOAT, 0, TAG_DATA, format->obj.comm);
  } else { // format->rank == 0
    int size;
    MPI_Comm_size(format->obj.comm, &size);
    for (int n = 0; n < size; n++) {
      struct mrc_f3 recv_fld;
      if (n == 0) {
	mrc_f3_alloc_with_array(&recv_fld, send_fld.ib, send_fld.im, 1, send_fld.arr);
      } else {
	int iw[6], *ib = iw, *im = iw + 3;
	MPI_Recv(iw, 6, MPI_INT, n, TAG_OFF_DIMS, format->obj.comm, MPI_STATUS_IGNORE);
	mrc_f3_alloc(&recv_fld, ib, im, 1);
	MPI_Recv(recv_fld.arr, recv_fld.len, MPI_FLOAT, n, TAG_DATA, format->obj.comm,
		 MPI_STATUS_IGNORE);
      }
      
      mrc_f3_foreach(&recv_fld, ix,iy,iz, 0, 0) {
	MRC_F3(lfld,0, ix,iy,iz) = MRC_F3(&recv_fld,0, ix,iy,iz);
      } mrc_f3_foreach_end;
      
      mrc_f3_free(&recv_fld);
    }
  }

  mrc_f3_free(&send_fld);
}

static void
ds_xdmf_to_one_write_field(struct diag_format *format, float scale, struct mrc_f3 *fld,
			   int m)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(format);

#ifndef NDEBUG
  MPI_Barrier(format->obj.comm);
#endif
  int gdims[3];
  struct mrc_domain *ldomain = NULL;
  struct mrc_f3 lfld;

  if (format->rank == 0) {
    mrc_domain_get_global_dims(fld->domain, gdims);
    struct mrc_domain_simple_params lpar = {
      .ldims    = { gdims[0], gdims[1], gdims[2] },
      .nr_procs = { 1, 1, 1 },
    };
    ldomain = mrc_domain_create(MPI_COMM_SELF, "simple");
    mrc_domain_simple_set_params(ldomain, &lpar);
    struct mrc_crds *crds = mrc_domain_get_crds(ldomain);
    mrc_crds_set_type(crds, "rectilinear");
    mrc_domain_setup(ldomain);
    float *crd[3];
    for (int d = 0; d < 3; d++) {
      crd[d] = calloc(gdims[d], sizeof(*crd[d]));
    }
    // FIXME, need something to just alloc space
    mrc_crds_set_values(crds, crd[0], gdims[0], crd[1], gdims[1], crd[2], gdims[2]);
    for (int d = 0; d < 3; d++) {
      free(crd[d]);
    }
    mrc_domain_f3_alloc(ldomain, &lfld, 1, SW_0);
  }

  communicate_fld(format, fld, m, scale, &lfld);

  if (!hdf5->crds_done) {
    communicate_crds(format, fld, &lfld);
    hdf5->crds_done = true;
  }

  if (format->rank != 0) {
    return;
  }

  char path[100]; make_path(path, "3df");

  struct xdmf_spatial *xs = xdmf_spatial_find(format, "3df", -1);
  if (!xs) {
    xs = xdmf_spatial_create_3d(format, gdims, 1);
    hid_t group0 = H5Gcreate(hdf5->file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group0);
    hdf5_write_crds(format, gdims, ldomain, fld->sw);
  }
  save_fld_info(xs, strdup(fld->name[m]), false);

  hid_t group0 = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  hsize_t hdims[3] = { lfld.im[2], lfld.im[1], lfld.im[0] };
  hid_t group = H5Gcreate(group0, fld->name[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTmake_dataset_float(group, "3d", 3, hdims, lfld.arr);
  H5Gclose(group);
  H5Gclose(group0);

  mrc_f3_free(&lfld);
}


static struct diag_format_ops ds_xdmf_to_one_ops = {
  .name          = "xdmf_to_one",
  .size          = sizeof(struct diag_hdf5),
  .parallel      = true,
  .create        = ds_xdmf_create,
  .destroy       = ds_xdmf_destroy,
  .open          = ds_xdmf_to_one_open,
  .close         = ds_xdmf_to_one_close,
  .write_field   = ds_xdmf_to_one_write_field,
};

#ifdef H5_HAVE_PARALLEL

// ======================================================================
// xdmf_parallel

static void
ds_xdmf_parallel_open(struct diag_format *format, float sheet, int outtype, int step) 
{
#ifndef NDEBUG
  MPI_Barrier(format->obj.comm);
#endif
  
  struct diag_hdf5 *hdf5 = diag_hdf5(format);
  hdf5->step = step;
  hdf5->crds_done = false;

  assert(outtype == DIAG_TYPE_3D);
  struct diag_info *diag_info = &format->diag_info;
  char filename[strlen(diag_info->outdir) + strlen(diag_info->run) + 30];
  sprintf(filename, "%s/%s.3df.%06d_p%06d.h5", diag_info->outdir, diag_info->run,
	  step, 0);

  hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist, format->obj.comm, MPI_INFO_NULL);
  hdf5->file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
  H5Pclose(plist);

  hdf5->group_crd = H5Gcreate(hdf5->file, "crd", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for (int d = 0; d < 3; d++) {
    hdf5->crd_written[d] = false;
  }

  INIT_LIST_HEAD(&hdf5->xdmf_spatial_list);
}

static void
ds_xdmf_parallel_close(struct diag_format *format)
{
#ifndef NDEBUG
  MPI_Barrier(format->obj.comm);
#endif
  
  struct diag_hdf5 *hdf5 = diag_hdf5(format);

  H5Gclose(hdf5->group_crd);
  hdf5_close(format);

  if (format->rank != 0)
    return;

  while (!list_empty(&hdf5->xdmf_spatial_list)) {
    struct xdmf_spatial *xs = list_entry(hdf5->xdmf_spatial_list.next, typeof(*xs), list);

    xdmf_write_spatial_collection(format, hdf5->step, xs);
    xdmf_write_temporal_collection(format, xs->sfx, hdf5->step);
    xdmf_spatial_destroy(xs);
  }
}

static void
hdf5_write_crds_parallel(struct diag_format *format, struct mrc_f3 *fld)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(format);

  const char *xyz[3] = { "x", "y", "z" };

  int gdims[3], off[3], ldims[3], nr_procs[3], idx[3];
  mrc_domain_get_global_dims(fld->domain, gdims);
  mrc_domain_get_local_offset_dims(fld->domain, off, ldims);
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
      float *crd = crds->crd[d] + fld->sw;
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

static void
ds_xdmf_parallel_write_field(struct diag_format *format, float scale, struct mrc_f3 *fld,
			   int m)
{
  struct diag_hdf5 *hdf5 = diag_hdf5(format);

#ifndef NDEBUG
  MPI_Barrier(format->obj.comm);
#endif

  int gdims[3], off[3], ldims[3];
  mrc_domain_get_global_dims(fld->domain, gdims);
  mrc_domain_get_local_offset_dims(fld->domain, off, ldims);

  struct mrc_f3 lfld; // strip boundary, could be done through hyperslab, but
                      // still have to scale, anyway
  mrc_f3_alloc(&lfld, NULL, ldims, 1);
  copy_and_scale(&lfld, 0, fld, m, scale, fld->sw);

  char path[100]; make_path(path, "3df");

  if (!hdf5->crds_done) {
    hid_t group0 = H5Gcreate(hdf5->file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group0);
    hdf5_write_crds_parallel(format, fld);
    hdf5->crds_done = true;
  }

  if (format->rank == 0) {
    struct xdmf_spatial *xs = xdmf_spatial_find(format, "3df", -1);
    if (!xs) {
      xs = xdmf_spatial_create_3d(format, gdims, 1);
    }
    save_fld_info(xs, strdup(fld->name[m]), false);
  }

  hid_t group0 = H5Gopen(hdf5->file, path, H5P_DEFAULT);
  hsize_t hgdims[3] = { gdims[2], gdims[1], gdims[0] };
  hsize_t hldims[3] = { ldims[2], ldims[1], ldims[0] };
  hsize_t hoff[3]   = { off[2], off[1], off[0] };
  hid_t group = H5Gcreate(group0, fld->name[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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
  H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, dxpl, lfld.arr);
  H5Dclose(dset);
  H5Sclose(filespace);
  H5Sclose(memspace);
  H5Pclose(dxpl);

  H5Gclose(group);
  H5Gclose(group0);

  mrc_f3_free(&lfld);
}

static struct diag_format_ops ds_xdmf_parallel_ops = {
  .name          = "xdmf_parallel",
  .size          = sizeof(struct diag_hdf5),
  .param_descr   = diag_hdf5_params_descr,
  .param_offset  = offsetof(struct diag_hdf5, par),
  .parallel      = true,
  .create        = ds_xdmf_create,
  .destroy       = ds_xdmf_destroy,
  .open          = ds_xdmf_parallel_open,
  .close         = ds_xdmf_parallel_close,
  .write_field   = ds_xdmf_parallel_write_field,
};

#endif

// ======================================================================

void
libmrc_diag_hdf5_register()
{
  libmrc_diag_register_format(&ds_xdmf_ops);
  libmrc_diag_register_format(&ds_xdmf_serial_ops);
  libmrc_diag_register_format(&ds_xdmf_to_one_ops);
#ifdef H5_HAVE_PARALLEL
  libmrc_diag_register_format(&ds_xdmf_parallel_ops);
#endif
}

#endif

