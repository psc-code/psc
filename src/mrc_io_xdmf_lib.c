
#include "mrc_io_xdmf_lib.h"

#include <mrc_io_private.h>
#include <mrc_fld.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

// ======================================================================
// xdmf writing

static void
xdmf_write_header(FILE *f)
{
  fprintf(f, "<?xml version='1.0' ?>\n");
  fprintf(f, "<Xdmf xmlns:xi='http://www.w3.org/2001/XInclude' Version='2.0'>\n");
}

static void
xdmf_write_topology_m1(FILE *f, int im[3], const char *filename, int p)
{
  // FIXME crd0 hardcoded, should use mrc_m1_name()
  fprintf(f, "     <Topology TopologyType=\"3DRectMesh\" Dimensions=\"%d %d %d\"/>\n",
	  1, 1, im[0]);

  fprintf(f, "     <Geometry GeometryType=\"VXVYVZ\">\n");
  fprintf(f, "     <DataItem Name=\"VX\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[0]);
  fprintf(f, "        %s:/crd0/p%d/1d\n", filename, p);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VY\" DataType=\"Float\" Dimensions=\"%d\" Format=\"XML\">\n", 1);
  fprintf(f, "        0\n");
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VZ\" DataType=\"Float\" Dimensions=\"%d\" Format=\"XML\">\n", 1);
  fprintf(f, "        0\n");
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     </Geometry>\n");
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Node\">\n",
	  "crdx", "Scalar");
  fprintf(f, "       <DataItem Dimensions=\"%d %d %d%s\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
	  1, 1, im[0], "");
  fprintf(f, "        %s:/crd0/p%d/1d\n", filename, p);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
  fprintf(f, "\n");
}

static void
xdmf_write_topology_uniform_m3(FILE *f, int im[3], float xl[3], float dx[3], int _sw)
{
  int sw[3];
  for (int d = 0; d < 3; d++) {
    sw[d] = im[d] > 1 ? _sw : 0;
  }
  fprintf(f, "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d\"/>\n",
	  im[2] + 1 + 2*sw[2], im[1] + 1 + 2*sw[1], im[0] + 1 + 2*sw[0]);

  fprintf(f, "     <Geometry GeometryType=\"Origin_DxDyDz\">\n");
  fprintf(f, "     <DataItem Name=\"Origin\" DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">\n");
  fprintf(f, "        %g %g %g\n", xl[2] - sw[2] * dx[2], xl[1] - sw[1] * dx[1], xl[0] - sw[0] * dx[0]);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"DxDyDz\" DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">\n");
  fprintf(f, "        %g %g %g\n", dx[2], dx[1], dx[0]);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     </Geometry>\n");
  fprintf(f, "\n");
}

static void
xdmf_write_topology_m3(FILE *f, int im[3], const char *filename, const char *crd_nc_path[3], int p)
{
  // FIXME crd[012] hardcoded, should use mrc_m1_name()
  fprintf(f, "     <Topology TopologyType=\"3DRectMesh\" Dimensions=\"%d %d %d\"/>\n",
	  im[2] + 1, im[1] + 1, im[0] + 1);

  fprintf(f, "     <Geometry GeometryType=\"VXVYVZ\">\n");
  fprintf(f, "     <DataItem Name=\"VX\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[0] + 1);
  fprintf(f, "        %s:/%s/crd_nc[0]/p%d/1d\n", filename, crd_nc_path[0], p);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VY\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[1] + 1);
  fprintf(f, "        %s:/%s/crd_nc[1]/p%d/1d\n", filename, crd_nc_path[1], p);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     <DataItem Name=\"VZ\" DataType=\"Float\" Dimensions=\"%d\" Format=\"HDF\">\n", im[2] + 1);
  fprintf(f, "        %s:/%s/crd_nc[2]/p%d/1d\n", filename, crd_nc_path[2], p);
  fprintf(f, "     </DataItem>\n");
  fprintf(f, "     </Geometry>\n");
  fprintf(f, "\n");
}

static void
xdmf_write_fld_m1(FILE *f, struct xdmf_fld_info *fld_info, int im[3],
		  const char *filename, int p)
{
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Node\">\n",
	  fld_info->name, fld_info->is_vec ? "Vector" : "Scalar");
  fprintf(f, "       <DataItem Dimensions=\"%d %d %d%s\" NumberType=\"%s\" Precision=\"%d\" Format=\"HDF\">\n",
	  1, 1, im[0], fld_info->is_vec ? " 3" : "", fld_info->dtype, fld_info->precision);
  fprintf(f, "        %s:/%s/%s/p%d/1d\n", filename, fld_info->path, fld_info->name, p);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
}

static void
xdmf_write_fld_m3(FILE *f, struct xdmf_fld_info *fld_info, int im[3],
		  const char *filename, int p)
{
  fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Cell\">\n",
	  fld_info->name, fld_info->is_vec ? "Vector" : "Scalar");
  fprintf(f, "       <DataItem Dimensions=\"%d %d %d%s\" NumberType=\"%s\" Precision=\"%d\" Format=\"HDF\">\n",
	  im[2], im[1], im[0], fld_info->is_vec ? " 3" : "", fld_info->dtype, fld_info->precision);
  fprintf(f, "        %s:/%s/%s/p%d/3d\n", filename, fld_info->path, fld_info->name, p);
  fprintf(f, "       </DataItem>\n");
  fprintf(f, "     </Attribute>\n");
}

// ======================================================================
// xdmf_temporal

struct xdmf_temporal *
xdmf_temporal_create(const char *filename)
{
  struct xdmf_temporal *xt = calloc(1, sizeof(*xt));
  INIT_LIST_HEAD(&xt->timesteps);
  xt->filename = strdup(filename);
  return xt;
}

void
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

void
xdmf_temporal_append(struct xdmf_temporal *xt, const char *fname_spatial)
{
  struct xdmf_temporal_step *xt_step =
    malloc(sizeof(*xt_step) + strlen(fname_spatial) + 1);
  strcpy(xt_step->filename, fname_spatial);
  list_add_tail(&xt_step->entry, &xt->timesteps);
}

void
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
  __list_for_each_entry(xt_step, &xt->timesteps, entry, struct xdmf_temporal_step) {
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

void
xdmf_spatial_open(list_t *xdmf_spatial_list)
{
  INIT_LIST_HEAD(xdmf_spatial_list);
}

struct xdmf_spatial *
xdmf_spatial_find(list_t *xdmf_spatial_list, const char *name)
{
  struct xdmf_spatial *xs;
  __list_for_each_entry(xs, xdmf_spatial_list, entry, struct xdmf_spatial) {
    if (strcmp(xs->name, name) == 0) {
      return xs;
    }
  }
  return NULL;
}

struct xdmf_spatial *
xdmf_spatial_create_f1(list_t *xdmf_spatial_list, const char *name, 
		       struct mrc_domain *domain)
{
  // OPT, we could skip this on procs which aren't writing xdmf

  struct xdmf_spatial *xs = calloc(1, sizeof(*xs));
  xs->name = strdup(name);
  xs->dim = 1;
  xs->nr_global_patches = 1; // FIXME wrong when parallel
  xs->patch_infos = calloc(xs->nr_global_patches, sizeof(*xs->patch_infos));

  for (int gp = 0; gp < xs->nr_global_patches; gp++) {
    mrc_domain_get_local_patch_info(domain, gp, &xs->patch_infos[gp]);
    //    mrc_domain_get_global_patch_info(domain, gp, &xs->patch_infos[gp]);
  }

  list_add_tail(&xs->entry, xdmf_spatial_list);
  return xs;
}

struct xdmf_spatial *
xdmf_spatial_create_m3(list_t *xdmf_spatial_list, const char *name, 
		       struct mrc_domain *domain, struct mrc_io *io, int sw)
{
  // OPT, we could skip this on procs which aren't writing xdmf

  struct xdmf_spatial *xs = calloc(1, sizeof(*xs));
  xs->name = strdup(name);
  xs->dim = 3;
  xs->sw = sw;
  mrc_domain_get_nr_global_patches(domain, &xs->nr_global_patches);
  xs->patch_infos = calloc(xs->nr_global_patches, sizeof(*xs->patch_infos));
  for (int d = 0; d < 3; d++) {
    xs->xl[d] = calloc(xs->nr_global_patches, sizeof(*xs->xl[d]));
    xs->dx[d] = calloc(xs->nr_global_patches, sizeof(*xs->dx[d]));
  }

  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  double xl[3];
  double dx[3];
  if (strcmp(mrc_crds_type(crds), "amr_uniform") == 0 ||
      strcmp(mrc_crds_type(crds), "uniform") == 0) {
    xs->uniform = true;
    mrc_crds_get_param_double3(crds, "l", xl);
    mrc_crds_get_dx_base(crds, dx);
  } else {
    for (int d = 0; d < 3; d++) {
      xs->crd_nc_path[d] = mrc_io_obj_path(io, crds->crd_nc[d]);
    }
  }

  for (int gp = 0; gp < xs->nr_global_patches; gp++) {
    mrc_domain_get_global_patch_info(domain, gp, &xs->patch_infos[gp]);
    if (xs->uniform) {
      int level = xs->patch_infos[gp].level;
      for (int d = 0; d < 3; d++) {
	float refine = 1.f;
	if (xs->patch_infos[gp].ldims[d] > 1) {
	  refine = 1.f / (1 << level);
	}
	xs->xl[d][gp] = xl[d] + xs->patch_infos[gp].off[d] * dx[d] * refine;
	xs->dx[d][gp] = dx[d] * refine;
      }
    }
  }

  list_add_tail(&xs->entry, xdmf_spatial_list);
  return xs;
}

struct xdmf_spatial *
xdmf_spatial_create_m3_parallel(list_t *xdmf_spatial_list, const char *name, 
				struct mrc_domain *domain, int slab_off[3], 
				int slab_dims[3], struct mrc_io *io)
{
  // OPT, we could skip this on procs which aren't writing xdmf

  struct xdmf_spatial *xs = calloc(1, sizeof(*xs));
  xs->name = strdup(name);
  xs->dim = 3;
  xs->nr_global_patches = 1;
  xs->patch_infos = calloc(xs->nr_global_patches, sizeof(*xs->patch_infos));

  for (int d = 0; d < 3; d++) {
    xs->patch_infos[0].ldims[d] = slab_dims[d];
    xs->xl[d] = calloc(xs->nr_global_patches, sizeof(*xs->xl[d]));
    xs->dx[d] = calloc(xs->nr_global_patches, sizeof(*xs->dx[d]));
  }

  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  if (strcmp(mrc_crds_type(crds), "amr_uniform") == 0 ||
      strcmp(mrc_crds_type(crds), "uniform") == 0) {
    xs->uniform = true;
    const double *xl = mrc_crds_lo(crds);
    double dx[3];
    mrc_crds_get_dx(crds, 0, dx);

    for (int d = 0; d < 3; d++) {
      double xnorm = crds->xnorm;
      xs->xl[d][0] = (xl[d] + slab_off[d] * dx[d]) * xnorm;
      xs->dx[d][0] = dx[d] * xnorm;
    }
  } else {
    for (int d = 0; d < 3; d++) {
      xs->crd_nc_path[d] = strdup(mrc_io_obj_path(io, crds->crd_nc[d]));
    }
  }

  list_add_tail(&xs->entry, xdmf_spatial_list);
  return xs;
}

void
xdmf_spatial_save_fld_info(struct xdmf_spatial *xs, char *fld_name,
			   char *path, bool is_vec, int mrc_dtype)
{
  assert(xs->nr_fld_info < MAX_XDMF_FLD_INFO);
  struct xdmf_fld_info *fld_info = &xs->fld_info[xs->nr_fld_info++];
  fld_info->name = fld_name;
  fld_info->path = path;
  fld_info->is_vec = is_vec;
  switch (mrc_dtype) {
  case MRC_NT_FLOAT: 
  {
    fld_info->dtype = strdup("Float");
    fld_info->precision = 4;
    break;
  }
  case MRC_NT_DOUBLE: 
  {
    fld_info->dtype = strdup("Float");
    fld_info->precision = 8;
    break;
  }
  case MRC_NT_INT: 
  {
    fld_info->dtype = strdup("Int");
    fld_info->precision = 4;
    break;
  }
  default: assert(0);
  }
}


void
xdmf_spatial_write(struct xdmf_spatial *xs, const char *filename,
		   struct mrc_io *io)
{
  if (io->rank != 0)
    return;

  FILE *f = fopen(filename, "w");
  xdmf_write_header(f);
  fprintf(f, "<Domain>\n");
  if (xs->nr_global_patches > 1) {
    fprintf(f, "<Grid GridType=\"Collection\" CollectionType=\"Spatial\">\n");
    fprintf(f, "<Time Type=\"Single\" Value=\"%g\" />\n", io->time);
  }
  for (int s = 0; s < xs->nr_global_patches; s++) {
    fprintf(f, "   <Grid Name=\"patch-%s-%d\" GridType=\"Uniform\">\n", xs->name, s);
    if (xs->nr_global_patches == 1) {
      fprintf(f, "<Time Type=\"Single\" Value=\"%g\" />\n", io->time);
    }
    
    char fname[strlen(io->par.outdir) + strlen(io->par.basename) + 20];
    int rank = xs->patch_infos[s].rank;
    int patch = xs->patch_infos[s].patch;
    sprintf(fname, "%s.%06d_p%06d.h5", io->par.basename, io->step, rank);
    int *ldims = xs->patch_infos[s].ldims;
    if (xs->dim == 1) {
      if (xs->uniform) {
	assert(0);
      } else {
	xdmf_write_topology_m1(f, ldims, fname, patch);
      }
    } else if (xs->dim == 3) {
      if (xs->uniform) {
	float xl[3] = { xs->xl[0][s], xs->xl[1][s], xs->xl[2][s] };
	float dx[3] = { xs->dx[0][s], xs->dx[1][s], xs->dx[2][s] };
	xdmf_write_topology_uniform_m3(f, ldims, xl, dx, xs->sw);
      } else {
	xdmf_write_topology_m3(f, ldims, fname, xs->crd_nc_path, patch);
      }
    } else {
      assert(0);
    }
    
    for (int m = 0; m < xs->nr_fld_info; m++) {
      if (xs->dim == 1) {
	xdmf_write_fld_m1(f, &xs->fld_info[m], ldims, fname, patch);
      } else if (xs->dim == 3) {
	xdmf_write_fld_m3(f, &xs->fld_info[m], ldims, fname, patch);
      } else {
	assert(0);
      }
    }
    fprintf(f, "   </Grid>\n");
  }
  if (xs->nr_global_patches > 1) {
    fprintf(f, "</Grid>\n");
  }
  fprintf(f, "</Domain>\n");
  fprintf(f, "</Xdmf>\n");
  fclose(f);
}

void
xdmf_spatial_destroy(struct xdmf_spatial *xs)
{
  for (int m = 0; m < xs->nr_fld_info; m++) {
    struct xdmf_fld_info *fld_info = &xs->fld_info[m];
    free(fld_info->name);
    free(fld_info->path);
    free(fld_info->dtype);
  }
  list_del(&xs->entry);
  free(xs->name);
  free(xs->patch_infos);
  for (int d = 0; d < 3; d++) {
    free(xs->xl[d]);
    free(xs->dx[d]);
  }
  free(xs);
}

void
xdmf_spatial_close(list_t *xdmf_spatial_list, struct mrc_io *io,
		   struct xdmf_temporal *xt)
{
  while (!list_empty(xdmf_spatial_list)) {
    struct xdmf_spatial *xs = list_entry(xdmf_spatial_list->next, struct xdmf_spatial, entry);
    
    char fname_spatial[strlen(io->par.outdir) + strlen(io->par.basename) + 20];
    sprintf(fname_spatial, "%s/%s.%06d.xdmf", io->par.outdir, io->par.basename, io->step);
    xdmf_spatial_write(xs, fname_spatial, io);
    xdmf_spatial_destroy(xs);

    if (io->rank == 0) {
      sprintf(fname_spatial, "%s.%06d.xdmf", io->par.basename, io->step);
      xdmf_temporal_append(xt, fname_spatial);
      xdmf_temporal_write(xt);
    }
  }
}

