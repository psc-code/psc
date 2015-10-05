
#include "mrc_domain_private.h"
#include "mrc_params.h"
#include "mrc_ddc.h"
#include "mrc_io.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define TAG_SCAN_OFF (1000)

static int
sfc_level_idx3_to_idx(struct mrc_domain_amr *amr, int l, int idx3[3])
{
  int i3[3];
  for (int d = 0; d < 3; d++) {
    if (amr->gdims[d] == 1) {
      i3[d] = 0;
    } else {
      i3[d] = idx3[d] << (amr->nr_levels - l);
    }
  }
  return sfc_idx3_to_idx(&amr->sfc, i3);
}

static void
sfc_level_idx_to_idx3(struct mrc_domain_amr *amr, int l, int sfc_idx, int idx3[3])
{
  sfc_idx_to_idx3(&amr->sfc, sfc_idx, idx3);
  for (int d = 0; d < 3; d++) {
    if (amr->gdims[d] != 1) {
      idx3[d] >>= (amr->nr_levels - l);
    }
  }
}

// ======================================================================
// map
// 
// maps between global patch index (contiguous) and 1D SFC idx
// (potentially non-contiguous)

static int
compare_level_sfc_idx(const void *p1, const void *p2)
{
  const struct mrc_amr_level_sfc_idx *s1 = p1, *s2 = p2;
  if (s1->sfc_idx < s2->sfc_idx) {
    return -1;
  } else if (s1->sfc_idx > s2->sfc_idx) {
    return 1;
  } else { // sfc_idx equal
    if (s1->l < s2->l) {
      return -1;
    } else if (s1->l > s2->l) {
      return 1;
    } else { // duplicates are not allowd
      assert(0);
    }
  }
}

static void
map_create(struct mrc_domain *domain)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  amr->map_gpatch_to_sfc = calloc(amr->nr_global_patches,
				  sizeof(*amr->map_gpatch_to_sfc));

  int gp = 0;
  while (!list_empty(&amr->global_patch_list)) {
    struct mrc_amr_patch *patch = list_entry(amr->global_patch_list.next,
					     struct mrc_amr_patch, entry);
    amr->map_gpatch_to_sfc[gp].l = patch->l;
    amr->map_gpatch_to_sfc[gp].sfc_idx = sfc_level_idx3_to_idx(amr, patch->l, patch->idx3);
    gp++;

    list_del(&patch->entry);
    free(patch);
  }
  qsort(amr->map_gpatch_to_sfc, amr->nr_global_patches,
	sizeof(*amr->map_gpatch_to_sfc), compare_level_sfc_idx);
}

static void
map_destroy(struct mrc_domain *domain)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  free(amr->map_gpatch_to_sfc);
}

static int
map_sfc_idx_to_gpatch(struct mrc_domain *domain, int l, int sfc_idx)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  // FIXME, do binary search
  for (int gp = 0; gp < amr->nr_global_patches; gp++) {
    if (amr->map_gpatch_to_sfc[gp].l == l &&
	amr->map_gpatch_to_sfc[gp].sfc_idx == sfc_idx) {
      return gp;
    }
  }
  return -1;
}

static void
map_gpatch_to_sfc_idx(struct mrc_domain *domain, int gpatch, int *l, int *sfc_idx)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  *l       = amr->map_gpatch_to_sfc[gpatch].l;
  *sfc_idx = amr->map_gpatch_to_sfc[gpatch].sfc_idx;
}

// ======================================================================

static void
gpatch_to_rank_patch(struct mrc_domain *domain, int gpatch,
		     int *rank, int *patch)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  // FIXME, this can be done much more efficiently using binary search...
  for (int i = 0; i < domain->size; i++) {
    if (gpatch < amr->gpatch_off_all[i+1]) {
      *rank = i;
      *patch = gpatch - amr->gpatch_off_all[i];
      break;
    }
  }
}

// ======================================================================
// mrc_domain_amr

// ----------------------------------------------------------------------
// mrc_domain_amr_get_global_patch_info

static void
mrc_domain_amr_get_global_patch_info(struct mrc_domain *domain, int gpatch,
				     struct mrc_patch_info *info)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  assert(gpatch < amr->nr_global_patches);
  info->global_patch = gpatch;
  gpatch_to_rank_patch(domain, gpatch, &info->rank, &info->patch);
  int sfc_idx;
  map_gpatch_to_sfc_idx(domain, gpatch, &info->level, &sfc_idx);
  
  assert(info->rank >= 0);
  
  int p3[3];
  sfc_level_idx_to_idx3(amr, info->level, sfc_idx, p3);
  for (int d = 0; d < 3; d++) {
    info->ldims[d] = amr->gdims[d];
    info->off[d] = p3[d] * amr->gdims[d];
    info->idx3[d] = p3[d];
  }
}

// ----------------------------------------------------------------------
// mrc_domain_amr_get_local_patch_info

static void
mrc_domain_amr_get_local_patch_info(struct mrc_domain *domain, int patch,
				      struct mrc_patch_info *info)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  mrc_domain_amr_get_global_patch_info(domain, amr->gpatch_off + patch,
				       info);
}

// ----------------------------------------------------------------------
// mrc_domain_amr_create

static void
mrc_domain_amr_create(struct mrc_domain *domain)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  INIT_LIST_HEAD(&amr->global_patch_list);
  amr->nr_levels = -1;
  // FIXME, this doesn't make too much sense, since a proper ddc is
  // problem / field specific
  mrc_ddc_set_type(domain->ddc, "amr");
  mrc_ddc_set_domain(domain->ddc, domain);
}

// ----------------------------------------------------------------------
// mrc_domain_amr_add_patch
//
// needs to be called collectively

static void
mrc_domain_amr_add_patch(struct mrc_domain *domain, int l, int idx3[3])
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  struct mrc_amr_patch *patch = calloc(1, sizeof(*patch));
  patch->l = l;
  for (int d = 0; d < 3; d++) {
    patch->idx3[d] = idx3[d];
  }
  if (l > amr->nr_levels) {
    amr->nr_levels = l;
  }

  list_add_tail(&patch->entry, &amr->global_patch_list);
  amr->nr_global_patches++;
}

static void
setup_gpatch_off_all(struct mrc_domain *domain)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  amr->gpatch_off_all = calloc(domain->size + 1, sizeof(*amr->gpatch_off_all));
  int nr_global_patches = amr->nr_global_patches;

  if (amr->nr_patches >= 0) {
    // prescribed mapping patch <-> proc
    MPI_Comm comm = mrc_domain_comm(domain);
    int *nr_patches_all = calloc(domain->size, sizeof(*nr_patches_all));
    MPI_Gather(&amr->nr_patches, 1, MPI_INT, nr_patches_all, 1, MPI_INT,
	       0, comm);
    MPI_Bcast(nr_patches_all, domain->size, MPI_INT, 0, comm);

    for (int i = 1; i <= domain->size; i++) {
      int nr_patches = nr_patches_all[i-1];
      amr->gpatch_off_all[i] = amr->gpatch_off_all[i-1] + nr_patches;
    }
    free(nr_patches_all);
  } else {
    // map patch <-> proc uniformly (roughly)
    int patches_per_proc = nr_global_patches / domain->size;
    int patches_per_proc_rmndr = nr_global_patches % domain->size;
    for (int i = 1; i <= domain->size; i++) {
      int nr_patches = patches_per_proc + ((i-1) < patches_per_proc_rmndr);
      amr->gpatch_off_all[i] = amr->gpatch_off_all[i-1] + nr_patches;
    }
  }
}

// ----------------------------------------------------------------------
// mrc_domain_amr_setup

static void
mrc_domain_amr_setup(struct mrc_domain *domain)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  MPI_Comm comm = mrc_domain_comm(domain);
  MPI_Comm_rank(comm, &domain->rank);
  MPI_Comm_size(comm, &domain->size);

  int sfc_np[3];
  for (int d = 0; d < 3; d++) {
    if (amr->gdims[d] == 1) {
      sfc_np[d] = 1;
    } else {
      sfc_np[d] = 1 << amr->nr_levels;
    }
  }
  sfc_setup(&amr->sfc, sfc_np);
  map_create(domain);

  setup_gpatch_off_all(domain);

  amr->gpatch_off = amr->gpatch_off_all[domain->rank];
  amr->nr_patches = amr->gpatch_off_all[domain->rank+1] - amr->gpatch_off;

  amr->patches = calloc(amr->nr_patches, sizeof(*amr->patches));
  for (int p = 0; p < amr->nr_patches; p++) {
    struct mrc_patch_info info;
    mrc_domain_amr_get_global_patch_info(domain, p + amr->gpatch_off, &info);
    for (int d = 0; d < 3; d++) {
      amr->patches[p].ldims[d] = info.ldims[d];
      amr->patches[p].off[d] = info.off[d];
    }
  }

  mrc_domain_setup_super(domain);
}

// ----------------------------------------------------------------------
// mrc_domain_amr_destroy

static void
mrc_domain_amr_destroy(struct mrc_domain *domain)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  free(amr->gpatch_off_all);
  free(amr->patches);
  map_destroy(domain);
}

// ----------------------------------------------------------------------
// mrc_domain_amr_view

static void
mrc_domain_amr_view(struct mrc_domain *domain)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  if (domain->rank == 0) {
    mprintf("nr_levels = %d\n", amr->nr_levels);
    mprintf("nr_global_patches = %d\n", amr->nr_global_patches);
    for (int gp = 0; gp < amr->nr_global_patches; gp++) {
      mprintf("map: gp = %d l = %d sfc_idx = %d\n", gp,
	      amr->map_gpatch_to_sfc[gp].l,
	      amr->map_gpatch_to_sfc[gp].sfc_idx);
    }
  }

  if (domain->rank == 0) {
    for (int gp = 0; gp < amr->nr_global_patches; gp++) {
      struct mrc_patch_info info;
      mrc_domain_get_global_patch_info(domain, gp, &info);
      mprintf("info: rank = %d\n", info.rank);
      mprintf("      patch = %d\n", info.patch);
      mprintf("      global_patch = %d\n", info.global_patch);
      mprintf("      l = %d\n", info.level);
      mprintf("      idx3 = %d,%d,%d\n", info.idx3[0], info.idx3[1], info.idx3[2]);
    }
  }

#if 0
  for (int proc = 0; proc < domain->size; proc++) {
    if (domain->rank == proc) {
      for (int p = 0; p < amr->nr_patches; p++) {
	struct mrc_patch *patch = &amr->patches[p];
	mprintf("patch %d: ldims %dx%dx%d off %dx%dx%d\n", p,
		patch->ldims[0], patch->ldims[1], patch->ldims[2],
	      patch->off[0], patch->off[1], patch->off[2]);
      }
    }
    MPI_Barrier(mrc_domain_comm(domain));
  }
#endif
}

// ----------------------------------------------------------------------
// mrc_domain_amr_get_patches

static struct mrc_patch *
mrc_domain_amr_get_patches(struct mrc_domain *domain, int *nr_patches)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);
  if (nr_patches) {
    *nr_patches = amr->nr_patches;
  }
  return amr->patches;
}

// ----------------------------------------------------------------------
// mrc_domain_amr_get_global_dims

static void
mrc_domain_amr_get_global_dims(struct mrc_domain *domain, int *dims)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  for (int d = 0; d < 3; d++) {
    dims[d] = amr->gdims[d];
  }
}

// ----------------------------------------------------------------------
// mrc_domain_amr_get_nr_global_patches

static void
mrc_domain_amr_get_nr_global_patches(struct mrc_domain *domain, int *nr_global_patches)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  *nr_global_patches = amr->nr_global_patches;
}

// ----------------------------------------------------------------------
// mrc_domain_amr_get_level_idx3_patch_info

static void
mrc_domain_amr_get_level_idx3_patch_info(struct mrc_domain *domain, int level,
					 int idx3[3], struct mrc_patch_info *info)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  info->rank = -1;
  info->patch = -1;
  info->global_patch = -1;
  for (int d = 0; d < 3; d++) {
    if (idx3[d] < 0 || idx3[d] >= (1 << level))
      return;
  }
  
  int sfc_idx = sfc_level_idx3_to_idx(amr, level, idx3);
  int gpatch = map_sfc_idx_to_gpatch(domain, level, sfc_idx);
  if (gpatch < 0) {
    return;
  }
  mrc_domain_amr_get_global_patch_info(domain, gpatch, info);
}

// ----------------------------------------------------------------------
// mrc_domain_amr_get_neighbor_rank_patch

extern void
mrc_domain_get_neighbor_patch_same(struct mrc_domain *domain, int p,
				   int dx[3], int *p_nei);

static void
mrc_domain_amr_get_neighbor_rank_patch(struct mrc_domain *domain, int p, int dir[3],
				       int *nei_rank, int *nei_patch)
{
  mrc_domain_get_neighbor_patch_same(domain, p, dir, nei_patch);
  if (*nei_patch < 0) {
    *nei_rank = -1;
  } else {
    *nei_rank = 0;
  }
}

// ----------------------------------------------------------------------
// mrc_domain_amr_get_nr_levels

static void
mrc_domain_amr_get_nr_levels(struct mrc_domain *domain, int *p_nr_levels)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  *p_nr_levels = amr->nr_levels + 1;
}

// ----------------------------------------------------------------------
// mrc_domain_amr_write

static void
mrc_domain_amr_write(struct mrc_domain *domain, struct mrc_io *io)
{
  int nr_global_patches;
  mrc_domain_amr_get_nr_global_patches(domain, &nr_global_patches);
  mrc_io_write_int(io, domain, "nr_global_patches", nr_global_patches);

  for (int gp = 0; gp < nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_amr_get_global_patch_info(domain, gp, &info);
    char path[strlen(mrc_io_obj_path(io, domain)) + 10];
    sprintf(path, "%s/p%d", mrc_io_obj_path(io, domain), gp);
    mrc_io_write_attr_int3(io, path, "ldims", info.ldims);
    mrc_io_write_attr_int3(io, path, "off", info.off);
    mrc_io_write_attr_int3(io, path, "idx3", info.idx3);
    int l, sfc_idx;
    map_gpatch_to_sfc_idx(domain, gp, &l, &sfc_idx);
    mrc_io_write_attr_int(io, path, "level", l);
    mrc_io_write_attr_int(io, path, "sfc_idx", sfc_idx);
  }
}

// ----------------------------------------------------------------------
// mrc_domain_amr_plot

static void
mrc_domain_amr_plot(struct mrc_domain *domain)
{
  struct mrc_domain_amr *amr = mrc_domain_amr(domain);

  if (domain->rank != 0) {
    return;
  }

  if (amr->gdims[2] != 1) {
    mprintf("WARNING: plotting only x-y domain\n");
  }

  const char *name = mrc_domain_name(domain);
  char filename[strlen(name) + 20];
  sprintf(filename, "%s-patches.asc", name);
  FILE *file = fopen(filename, "w");
  sprintf(filename, "%s-curve.asc", name);
  FILE *file_curve = fopen(filename, "w");

  int nr_global_patches;
  mrc_domain_get_nr_global_patches(domain, &nr_global_patches);
  
  for (int gp = 0; gp < nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_get_global_patch_info(domain, gp, &info);
    float xb[3], xe[3];
    for (int d = 0; d < 3; d++) {
      xb[d] = (float) info.off[d] / (1 << info.level);
      xe[d] = (float) (info.off[d] + info.ldims[d]) / (1 << info.level);
    }
    fprintf(file, "%g %g %d %d\n", xb[0], xb[1], info.rank, info.level);
    fprintf(file, "%g %g %d %d\n", xe[0], xb[1], info.rank, info.level);
    fprintf(file, "\n");
    fprintf(file, "%g %g %d %d\n", xb[0], xe[1], info.rank, info.level);
    fprintf(file, "%g %g %d %d\n", xe[0], xe[1], info.rank, info.level);
    fprintf(file, "\n");
    fprintf(file, "\n");

    fprintf(file_curve, "%g %g %d %d\n", .5*(xb[0] + xe[0]), .5*(xb[1] + xe[1]),
	    info.rank, info.level);
  }

  fclose(file);
  fclose(file_curve);
}

static struct mrc_ddc *
mrc_domain_amr_create_ddc(struct mrc_domain *domain)
{
  struct mrc_ddc *ddc = mrc_ddc_create(domain->obj.comm);
  mrc_ddc_set_type(ddc, "amr");
  mrc_ddc_set_domain(ddc, domain);
  return ddc;
}

static struct mrc_param_select curve_descr[] = {
  { .val = CURVE_BYDIM   , .str = "bydim"    },
  { .val = CURVE_MORTON  , .str = "morton"   },
  { .val = CURVE_HILBERT , .str = "hilbert"  },
  {},
};

#define VAR(x) (void *)offsetof(struct mrc_domain_amr, x)
static struct param mrc_domain_amr_params_descr[] = {
  { "m"               , VAR(gdims)           , PARAM_INT3(32, 32, 32) },
  { "curve_type"      , VAR(sfc.curve_type)  , PARAM_SELECT(CURVE_BYDIM,
							    curve_descr) },
  { "nr_patches"      , VAR(nr_patches)      , PARAM_INT(-1) },
  {},
};
#undef VAR

struct mrc_domain_ops mrc_domain_amr_ops = {
  .name                      = "amr",
  .size                      = sizeof(struct mrc_domain_amr),
  .param_descr               = mrc_domain_amr_params_descr,
  .create                    = mrc_domain_amr_create,
  .setup                     = mrc_domain_amr_setup,
  .view                      = mrc_domain_amr_view,
  .write                     = mrc_domain_amr_write,
  .destroy                   = mrc_domain_amr_destroy,
  .add_patch                 = mrc_domain_amr_add_patch,
  .get_patches               = mrc_domain_amr_get_patches,
  .get_nr_global_patches     = mrc_domain_amr_get_nr_global_patches,
  .get_global_dims           = mrc_domain_amr_get_global_dims,
  .get_global_patch_info     = mrc_domain_amr_get_global_patch_info,
  .get_local_patch_info      = mrc_domain_amr_get_local_patch_info,
  .get_level_idx3_patch_info = mrc_domain_amr_get_level_idx3_patch_info,
  .get_nr_levels             = mrc_domain_amr_get_nr_levels,
  .get_neighbor_rank_patch   = mrc_domain_amr_get_neighbor_rank_patch,
  .plot                      = mrc_domain_amr_plot,
  .create_ddc                = mrc_domain_amr_create_ddc,
};

