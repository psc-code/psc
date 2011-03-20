
#include "mrc_io_xdmf_lib.h"

#include <mrc_io_private.h>

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

void
xdmf_spatial_open(list_t *xdmf_spatial_list)
{
  INIT_LIST_HEAD(xdmf_spatial_list);
}

struct xdmf_spatial *
xdmf_spatial_find(list_t *xdmf_spatial_list, const char *name)
{
  struct xdmf_spatial *xs;
  list_for_each_entry(xs, xdmf_spatial_list, entry) {
    if (strcmp(xs->name, name) == 0) {
      return xs;
    }
  }
  return NULL;
}

struct xdmf_spatial *
xdmf_spatial_create_m3(list_t *xdmf_spatial_list, const char *name, 
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

  list_add_tail(&xs->entry, xdmf_spatial_list);
  return xs;
}

void
xdmf_spatial_save_fld_info(struct xdmf_spatial *xs, char *fld_name,
			   char *path, bool is_vec)
{
  assert(xs->nr_fld_info < MAX_XDMF_FLD_INFO);
  struct xdmf_fld_info *fld_info = &xs->fld_info[xs->nr_fld_info++];
  fld_info->name = fld_name;
  fld_info->path = path;
  fld_info->is_vec = is_vec;
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

void
xdmf_spatial_destroy(struct xdmf_spatial *xs)
{
  list_del(&xs->entry);
  free(xs->patch_infos);
  free(xs);
}

void
xdmf_spatial_close(list_t *xdmf_spatial_list, struct mrc_io *io,
		   struct xdmf_temporal *xt)
{
  while (!list_empty(xdmf_spatial_list)) {
    struct xdmf_spatial *xs = list_entry(xdmf_spatial_list->next, typeof(*xs), entry);
    
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

