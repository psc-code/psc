
#include "mrc_domain_private.h"
#include "mrc_params.h"
#include "mrc_ddc.h"
#include "mrc_io.h"
#include "mrc_io_private.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define TAG_SCAN_OFF (1000)

static inline struct mrc_domain_multi *
mrc_domain_multi(struct mrc_domain *domain)
{
  return domain->obj.subctx;
}

// ======================================================================
// map
// 
// maps between global patch index (contiguous) and 1D SFC idx
// (potentially non-contiguous)
// if sfc_indices is NULL, the map will be the indentity

static void
map_create(struct mrc_domain *domain, int *sfc_indices, int nr_gpatches)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  if (!sfc_indices) {
    return;
  }
  multi->gp = malloc(sizeof(int) * multi->nr_global_patches);
  for (int i = 0; i < nr_gpatches; i++) {
    multi->gp[i] = sfc_indices[i];
  }

  //Create the bintree for performant searching
  int vals[nr_gpatches];
  for (int i = 0; i < nr_gpatches; i++) {
    vals[i] = i;
  }
  bintree_create_from_ordered_list(&multi->g_patches, sfc_indices, vals, nr_gpatches);
}

static void
map_destroy(struct mrc_domain *domain)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  if (!multi->gp) {
    return;
  }
  free(multi->gp);
  bintree_destroy(&multi->g_patches);
}

static int
map_sfc_idx_to_gpatch(struct mrc_domain *domain, int sfc_idx)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  if (!multi->gp) {
    return sfc_idx;
  }
  int retval;
  int rc = bintree_get(&multi->g_patches, sfc_idx, &retval);
  if (rc == 0) {
    return -1;
  }
  return retval;
}

static int
map_gpatch_to_sfc_idx(struct mrc_domain *domain, int gpatch)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  if (!multi->gp) {
    return gpatch;
  }
  return multi->gp[gpatch];
}

// ======================================================================

static void
gpatch_to_rank_patch(struct mrc_domain *domain, int gpatch,
		     int *rank, int *patch)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  // FIXME, this can be done much more efficiently using binary search...
  for (int i = 0; i < domain->size; i++) {
    if (gpatch < multi->gpatch_off_all[i+1]) {
      *rank = i;
      *patch = gpatch - multi->gpatch_off_all[i];
      return;
    }
  }
  assert(0);
}

// ======================================================================

static void
mrc_domain_multi_view(struct mrc_domain *domain)
{
#if 0
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int proc = 0; proc < domain->size; proc++) {
    if (domain->rank == proc) {
      for (int p = 0; p < multi->nr_patches; p++) {
	struct mrc_patch *patch = &multi->patches[p];
	mprintf("patch %d: ldims %dx%dx%d off %dx%dx%d\n", p,
		patch->ldims[0], patch->ldims[1], patch->ldims[2],
	      patch->off[0], patch->off[1], patch->off[2]);
      }
    }
    MPI_Barrier(obj->comm);
  }
#endif
}

static void
mrc_domain_multi_get_global_patch_info(struct mrc_domain *domain, int gpatch,
				       struct mrc_patch_info *info)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  assert(gpatch < multi->nr_global_patches);
  info->global_patch = gpatch;
  gpatch_to_rank_patch(domain, gpatch, &info->rank, &info->patch);
  
  assert(info->rank >= 0);
  
  int sfc_idx = map_gpatch_to_sfc_idx(domain, gpatch);
  int p3[3];
  sfc_idx_to_idx3(&multi->sfc, sfc_idx, p3);
  info->level = 0;
  for (int d = 0; d < 3; d++) {
    info->ldims[d] = multi->ldims[d][p3[d]];
    info->off[d] = multi->off[d][p3[d]];
    info->idx3[d] = p3[d];
  }
}

static void
mrc_domain_multi_get_local_patch_info(struct mrc_domain *domain, int patch,
				      struct mrc_patch_info *info)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  mrc_domain_multi_get_global_patch_info(domain, multi->gpatch_off + patch,
					 info);
}

static void
setup_gpatch_off_all(struct mrc_domain *domain)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  multi->gpatch_off_all = calloc(domain->size + 1, sizeof(*multi->gpatch_off_all));
  int nr_global_patches = multi->nr_global_patches;

  if (multi->nr_patches >= 0) {
    // prescribed mapping patch <-> proc
    MPI_Comm comm = mrc_domain_comm(domain);
    int *nr_patches_all = calloc(domain->size, sizeof(*nr_patches_all));
    MPI_Gather(&multi->nr_patches, 1, MPI_INT, nr_patches_all, 1, MPI_INT,
	       0, comm);
    MPI_Bcast(nr_patches_all, domain->size, MPI_INT, 0, comm);

    for (int i = 1; i <= domain->size; i++) {
      int nr_patches = nr_patches_all[i-1];
      multi->gpatch_off_all[i] = multi->gpatch_off_all[i-1] + nr_patches;
    }
    free(nr_patches_all);
  } else {
    // map patch <-> proc uniformly (roughly)
    int patches_per_proc = nr_global_patches / domain->size;
    int patches_per_proc_rmndr = nr_global_patches % domain->size;
    for (int i = 1; i <= domain->size; i++) {
      int nr_patches = patches_per_proc + ((i-1) < patches_per_proc_rmndr);
      multi->gpatch_off_all[i] = multi->gpatch_off_all[i-1] + nr_patches;
    }
  }
}

static void
mrc_domain_multi_setup_map(struct mrc_domain *domain)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  if (!multi->have_activepatches) {
    map_create(domain, NULL, -1);
  } else {
    int sfc_indices[multi->nr_global_patches];
  
    int npatches = 0;
  
    //TODO Find a smarter way than iterating over all possible patches
    int npt = multi->np[0] * multi->np[1] * multi->np[2];
    for(int i = 0; i < npt; i++) {
      int idx[3];
      sfc_idx_to_idx3(&multi->sfc, i, idx);
      if(bitfield3d_isset(&multi->activepatches, idx)) {
	//Register the patch
	sfc_indices[npatches] = i;
	npatches++;
      }
    }
    map_create(domain, sfc_indices, multi->nr_global_patches);
  }
}

static void
mrc_domain_multi_do_setup(struct mrc_domain *domain)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  MPI_Comm comm = mrc_domain_comm(domain);
  MPI_Comm_rank(comm, &domain->rank);
  MPI_Comm_size(comm, &domain->size);

  int *np = multi->np;
  if (multi->p_activepatches) {
    //Copy the activepatch-list
    bitfield3d_copy(&multi->activepatches, multi->p_activepatches);
    multi->have_activepatches = true;
    multi->nr_global_patches = bitfield3d_count_bits_set(&multi->activepatches);
  } else {
    multi->nr_global_patches = np[0] * np[1] * np[2];
  }
  
  for (int d = 0; d < 3; d++) {
    int ldims[3], rmndr[3];
    ldims[d] = multi->gdims[d] / np[d];
    rmndr[d] = multi->gdims[d] % np[d];
    assert(rmndr[d] == 0); // we used to support this, but not anymore

    multi->ldims[d] = calloc(np[d], sizeof(*multi->ldims[d]));
    multi->off[d] = calloc(np[d], sizeof(*multi->off[d]));
    for (int i = 0; i < np[d]; i++) {
      multi->ldims[d][i] = ldims[d] + (i < rmndr[d]);
      if (i > 0) {
	multi->off[d][i] = multi->off[d][i-1] + multi->ldims[d][i-1];
      }
    }
  }

  sfc_setup(&multi->sfc, multi->np);
  mrc_domain_multi_setup_map(domain);

  setup_gpatch_off_all(domain);

  multi->gpatch_off = multi->gpatch_off_all[domain->rank];
  multi->nr_patches = multi->gpatch_off_all[domain->rank+1] - multi->gpatch_off;

  multi->patches = calloc(multi->nr_patches, sizeof(*multi->patches));
  for (int p = 0; p < multi->nr_patches; p++) {
    struct mrc_patch_info info;
    mrc_domain_multi_get_local_patch_info(domain, p, &info);
    for (int d = 0; d < 3; d++) {
      multi->patches[p].ldims[d] = info.ldims[d];
      multi->patches[p].off[d] = info.off[d];
    }
  }
}

static void
mrc_domain_multi_setup(struct mrc_domain *domain)
{
  mrc_domain_multi_do_setup(domain);
  mrc_domain_setup_super(domain);
}

static void
mrc_domain_multi_destroy(struct mrc_domain *domain)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int d = 0; d < 3; d++) {
    free(multi->ldims[d]);
    free(multi->off[d]);
  }
  free(multi->gpatch_off_all);
  free(multi->patches);
  map_destroy(domain);
  if (multi->have_activepatches) {
    bitfield3d_destroy(&multi->activepatches);
  }
  sfc_destroy(&multi->sfc);
}

static struct mrc_patch *
mrc_domain_multi_get_patches(struct mrc_domain *domain, int *nr_patches)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);
  if (nr_patches) {
    *nr_patches = multi->nr_patches;
  }
  return multi->patches;
}

static void
mrc_domain_multi_get_global_dims(struct mrc_domain *domain, int *dims)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int d = 0; d < 3; d++) {
    dims[d] = multi->gdims[d];
  }
}

static void
mrc_domain_multi_get_nr_procs(struct mrc_domain *domain, int *nr_procs)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  for (int d = 0; d < 3; d++) {
    nr_procs[d] = multi->np[d];
  }
}

static void
mrc_domain_multi_get_nr_levels(struct mrc_domain *domain, int *p_nr_levels)
{
  *p_nr_levels = 1;
}

static void
mrc_domain_multi_get_nr_global_patches(struct mrc_domain *domain, int *nr_global_patches)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  *nr_global_patches = multi->nr_global_patches;
}

static void
mrc_domain_multi_get_level_idx3_patch_info(struct mrc_domain *domain, int level,
					   int idx[3], struct mrc_patch_info *info)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  assert(level == 0);
  //Check if the patch is active
  if (multi->have_activepatches &&
      !bitfield3d_isset(&multi->activepatches, idx)) {
    info->rank = -1;
    info->patch = -1;
    info->global_patch = -1;
    for(int d = 0; d < 3; d++) {
      info->ldims[d] = multi->ldims[d][idx[d]];
      info->off[d] = multi->off[d][idx[d]];
      info->idx3[d] = idx[d];
    }
    return;
  }

  int sfc_idx = sfc_idx3_to_idx(&multi->sfc, idx);
  int gpatch = map_sfc_idx_to_gpatch(domain, sfc_idx);
  mrc_domain_multi_get_global_patch_info(domain, gpatch, info);
}

static void
mrc_domain_multi_write(struct mrc_domain *domain, struct mrc_io *io)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  int nr_global_patches;
  mrc_domain_multi_get_nr_global_patches(domain, &nr_global_patches);
  mrc_io_write_int(io, domain, "nr_global_patches", nr_global_patches);

  struct mrc_io_ops *io_ops = (struct mrc_io_ops *) io->obj.ops;
  if (io_ops->parallel) {
    return;
  }

  int mpi_rank;
  MPI_Comm_rank(mrc_domain_comm(domain), &mpi_rank);

  char path[strlen(mrc_io_obj_path(io, domain)) + 10];
  sprintf(path, "%s/rank_%d", mrc_io_obj_path(io, domain), mpi_rank);
  mrc_io_write_attr_int(io, path, "nr_patches", multi->nr_patches);

#if 0
  int mpi_size;
  MPI_Comm_size(mrc_domain_comm(domain), &mpi_size);

  char path[strlen(mrc_io_obj_path(io, domain)) + 10];
  int last_rank = 0, nr_patches_in_rank = 0;
  for (int gp = 0; gp < nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_multi_get_global_patch_info(domain, gp, &info);
    sprintf(path, "%s/p%d", mrc_io_obj_path(io, domain), gp);
    mrc_io_write_attr_int3(io, path, "ldims", info.ldims);
    mrc_io_write_attr_int3(io, path, "off", info.off);
    mrc_io_write_attr_int3(io, path, "idx3", info.idx3);
    int sfc_idx = map_gpatch_to_sfc_idx(domain, gp);
    mrc_io_write_attr_int(io, path, "sfc_idx", sfc_idx);

    while (last_rank < info.rank) {
      sprintf(path, "%s/rank_%d", mrc_io_obj_path(io, domain), last_rank);
      mrc_io_write_attr_int(io, path, "nr_patches", nr_patches_in_rank);
      last_rank++;
      nr_patches_in_rank = 0;
    }
    nr_patches_in_rank++;
  }
  while (last_rank < mpi_size) {
    sprintf(path, "%s/rank_%d", mrc_io_obj_path(io, domain), last_rank);
    mrc_io_write_attr_int(io, path, "nr_patches", nr_patches_in_rank);
    last_rank++;
    nr_patches_in_rank = 0;
  }
#endif
}

static void
mrc_domain_multi_read(struct mrc_domain *domain, struct mrc_io *io)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  // This isn't a collective value, so we better
  // don't take what's been read from the file
  int mpi_rank, mpi_size;
  MPI_Comm_rank(mrc_domain_comm(domain), &mpi_rank);
  MPI_Comm_size(mrc_domain_comm(domain), &mpi_size);

  struct mrc_io_ops *io_ops = (struct mrc_io_ops *) io->obj.ops;

  // FIXME: With the addition of **do_setup() below it may no longer
  // be necessary to have this read/write of local patch number
  // for non-parallel writers, as the structure will be regenerated
  // anyway.
  if(!io_ops->parallel) {
    const char *path = mrc_io_obj_path(io, domain);
    char path2[strlen(path) + 20];
    sprintf(path2, "%s/rank_%d", path, mpi_rank);
    mrc_io_read_attr_int(io, path2, "nr_patches", &multi->nr_patches);
  } else {
    multi->nr_patches = -1;
  }

#if 0
  // need to read attributes collectively :(
  // FIXME, maybe we should just write an array?
  for (int i = 0; i < mpi_size; i++) {
  char path2[strlen(path) + 20];
    sprintf(path2, "%s/rank_%d", path, i);
    int nr_patches;
    mrc_io_read_attr_int(io, path2, "nr_patches", &nr_patches);
    if (i == mpi_rank) {
      multi->nr_patches = nr_patches;
    }
  }
#endif

  // FIXME? We're mostly redoing things from scratch, rather
  // than reading them...
  // WARNING: This may break PSC checkpointing. (but I think setup was getting run anyway)
  mrc_domain_multi_do_setup(domain);
  mrc_domain_read_super(domain, io);
}

static void
mrc_domain_multi_plot(struct mrc_domain *domain)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  if (domain->rank != 0) {
    return;
  }

  if (multi->np[2] != 1) {
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
    fprintf(file, "%d %d %d\n", info.off[0], info.off[1], info.rank);
    fprintf(file, "%d %d %d\n", info.off[0] + info.ldims[0], info.off[1], info.rank);
    fprintf(file, "\n");
    fprintf(file, "%d %d %d\n", info.off[0], info.off[1] + info.ldims[1], info.rank);
    fprintf(file, "%d %d %d\n", info.off[0] + info.ldims[0], info.off[1] + info.ldims[1],
	    info.rank);
    fprintf(file, "\n");
    fprintf(file, "\n");

    fprintf(file_curve, "%g %g %d\n", info.off[0] + .5 * info.ldims[0],
	    info.off[1] + .5 * info.ldims[1], info.rank);
  }

  fclose(file);
  fclose(file_curve);
}

// ----------------------------------------------------------------------
// mrc_domain_multi_get_neighbor_rank_patch

static void
mrc_domain_multi_get_neighbor_rank_patch(struct mrc_domain *domain, int p, int dir[3],
					 int *nei_rank, int *nei_patch)
{
  struct mrc_domain_multi *multi = mrc_domain_multi(domain);

  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(domain, p, &info);
  int patch_idx_nei[3];
  for (int d = 0; d < 3; d++) {
    patch_idx_nei[d] = info.idx3[d] + dir[d];
    if (domain->bc[d] == BC_PERIODIC) {
      if (patch_idx_nei[d] < 0) {
	patch_idx_nei[d] += multi->np[d];
      }
      if (patch_idx_nei[d] >= multi->np[d]) {
	patch_idx_nei[d] -= multi->np[d];
      }
    }
    if (patch_idx_nei[d] < 0 || patch_idx_nei[d] >= multi->np[d]) {
      *nei_rank = -1;
      *nei_patch = -1;
      return;
    }
  }
  mrc_domain_get_level_idx3_patch_info(domain, 0, patch_idx_nei, &info);
  *nei_rank = info.rank;
  *nei_patch = info.patch;
}

static struct mrc_ddc *
mrc_domain_multi_create_ddc(struct mrc_domain *domain)
{
  struct mrc_ddc *ddc = mrc_ddc_create(domain->obj.comm);
  mrc_ddc_set_type(ddc, "multi");
  mrc_ddc_set_domain(ddc, domain);
  return ddc;
}

static struct mrc_param_select curve_descr[] = {
  { .val = CURVE_BYDIM   , .str = "bydim"    },
  { .val = CURVE_MORTON  , .str = "morton"   },
  { .val = CURVE_HILBERT , .str = "hilbert"  },
  {},
};

#define VAR(x) (void *)offsetof(struct mrc_domain_multi, x)
static struct param mrc_domain_multi_params_descr[] = {
  { "m"               , VAR(gdims)           , PARAM_INT3(32, 32, 32),
    .help = "Global dimensions of the domain (as number of cells per direction)." },
  { "np"              , VAR(np)              , PARAM_INT3(1, 1, 1),
    .help = "Number of patches to divide global domain into (per direction)" },
  { "curve_type"      , VAR(sfc.curve_type)  , PARAM_SELECT(CURVE_BYDIM,
							    curve_descr),
    .help = "Type of spacing filling curve to use for distributing patches." },
  { "nr_patches"      , VAR(nr_patches)      , PARAM_INT(-1) },
  { "activepatches"   , VAR(p_activepatches) , PARAM_PTR(NULL) },
  {},
};
#undef VAR

struct mrc_domain_ops mrc_domain_multi_ops = {
  .name                      = "multi",
  .size                      = sizeof(struct mrc_domain_multi),
  .param_descr               = mrc_domain_multi_params_descr,
  .setup                     = mrc_domain_multi_setup,
  .view                      = mrc_domain_multi_view,
  .write                     = mrc_domain_multi_write,
  .read                      = mrc_domain_multi_read,
  .destroy                   = mrc_domain_multi_destroy,
  .get_patches               = mrc_domain_multi_get_patches,
  .get_global_dims           = mrc_domain_multi_get_global_dims,
  .get_nr_procs              = mrc_domain_multi_get_nr_procs,
  .get_nr_levels             = mrc_domain_multi_get_nr_levels,
  .get_nr_global_patches     = mrc_domain_multi_get_nr_global_patches,
  .get_global_patch_info     = mrc_domain_multi_get_global_patch_info,
  .get_local_patch_info      = mrc_domain_multi_get_local_patch_info,
  .get_level_idx3_patch_info = mrc_domain_multi_get_level_idx3_patch_info,
  .get_neighbor_rank_patch   = mrc_domain_multi_get_neighbor_rank_patch,
  .plot                      = mrc_domain_multi_plot,
  .create_ddc                = mrc_domain_multi_create_ddc,
};

