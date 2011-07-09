
#include "mrc_domain_private.h"
#include "mrc_params.h"
#include "mrc_ddc.h"
#include "mrc_io.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define TAG_SCAN_OFF (1000)

static inline struct mrc_domain_dynamic *
mrc_domain_dynamic(struct mrc_domain *domain)
{
  return domain->obj.subctx;
}

// ======================================================================
// map
// 
// maps between global patch index (contiguous) and 1D SFC idx
// (potentially non-contiguous)

static void
map_create(struct mrc_domain *domain, int *sfc_indices, int nr_gpatches)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);

  multi->gp = malloc(sizeof(int) * multi->nr_gpatches);
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
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);

  free(multi->gp);
  bintree_destroy(&multi->g_patches);
}

static int
map_sfc_idx_to_gpatch(struct mrc_domain *domain, int sfc_idx)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);

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
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);

  return multi->gp[gpatch];
}

// ======================================================================

static void
sfc_idx_to_rank_patch(struct mrc_domain *domain, int sfc_idx,
		      int *rank, int *patch)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);
  //Get gp->array index
  int gpatch = map_sfc_idx_to_gpatch(domain, sfc_idx);
  if (gpatch < 0) {
    *rank = -1;
    *patch = -1;
    return;
  }
  
  // FIXME, multi can be done much more efficiently using binary search...
  for (int i = 0; i < domain->size; i++) {
    if (gpatch < multi->gpatch_off_all[i+1]) {
      *rank = i;
      *patch = gpatch - multi->gpatch_off_all[i];
      break;
    }
  }
}

// ======================================================================

static void
mrc_domain_dynamic_view(struct mrc_domain *domain)
{
}

// FIXME, get rid of multi one always use idx3 one?
static void
mrc_domain_dynamic_get_global_patch_info(struct mrc_domain *domain, int gpatch,
				       struct mrc_patch_info *info)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);

  assert(gpatch < multi->nr_gpatches);
  info->global_patch = gpatch;
  int sfc_idx = map_gpatch_to_sfc_idx(domain, gpatch);
  sfc_idx_to_rank_patch(domain, sfc_idx, &info->rank, &info->patch);
  
  assert(info->rank >= 0);

  int p3[3];
  sfc_idx_to_idx3(&multi->sfc, sfc_idx, p3);
  for (int d = 0; d < 3; d++) {
    info->ldims[d] = multi->_ldims[d];
    info->off[d] = p3[d] * multi->_ldims[d];
    info->idx3[d] = p3[d];
  }
}

static void
mrc_domain_dynamic_get_local_patch_info(struct mrc_domain *domain, int patch,
				      struct mrc_patch_info *info)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);

  mrc_domain_dynamic_get_global_patch_info(domain, multi->gpatch_off + patch, info);
}

static void
mrc_domain_dynamic_setup_patches(struct mrc_domain *domain, int *nr_patches_all)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);

  int sfc_indices[multi->nr_gpatches];
  
  multi->patches = malloc(sizeof(*multi->patches) * multi->nr_patches);
  
  multi->gpatch_off_all = calloc(domain->size + 1, sizeof(*multi->gpatch_off_all));
  for (int i = 1; i <= domain->size; i++) {
    multi->gpatch_off_all[i] = multi->gpatch_off_all[i-1] + nr_patches_all[i-1];
  }
 
  multi->gpatch_off = multi->gpatch_off_all[domain->rank];

  int activerank = 0;
  int npatches = 0;
  
  //TODO Find a smarter way than iterating over all possible patches
  int npt = multi->np[0] * multi->np[1] * multi->np[2];
  for(int i = 0; i < npt; i++) {
    int idx[3];
    sfc_idx_to_idx3(&multi->sfc, i, idx);
    if(bitfield3d_isset(&multi->activepatches, idx)) {
      //Calculate rank
      if (npatches >= multi->gpatch_off_all[activerank+1]) {
	activerank++;
      }
      
      //Register the patch
      sfc_indices[npatches] = i;

      if (activerank == domain->rank) { // Create the patch on multi processor
	//Setup patches[lpatch]
	int lpatch = npatches - multi->gpatch_off;
	for(int d = 0; d < 3; d++) {
	  multi->patches[lpatch].off[d] = idx[d] * multi->_ldims[d];
	  multi->patches[lpatch].ldims[d] = multi->_ldims[d];
	}
      }
      npatches++;
    }
  }
  
  map_create(domain, sfc_indices, multi->nr_gpatches);
}

static void
mrc_domain_dynamic_setup(struct mrc_domain *domain)
{
  assert(!domain->is_setup);
  domain->is_setup = true;

  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);
  	
  //Assert some properties we rely on
  for(int d=0; d<3; ++d) assert((multi->gdims[d] % multi->np[d]) == 0);	//The domain-size must be a mutliple of the patch-size
  
  //Copy the activepatch-list
  bitfield3d_copy(&multi->activepatches, multi->p_activepatches);
  
  multi->nr_gpatches = bitfield3d_count_bits_set(&multi->activepatches);
  
  MPI_Comm comm = mrc_domain_comm(domain);
  MPI_Comm_rank(domain->obj.comm, &domain->rank);
  MPI_Comm_size(domain->obj.comm, &domain->size);
  
  //If nr_patches is not set, try to distribute patches evenly
  if(multi->nr_patches < 0) {
    multi->nr_patches = multi->nr_gpatches / domain->size +
      (domain->rank < (multi->nr_gpatches % domain->size) ? 1 : 0);
  }
  
  int *nr_patches_all = calloc(domain->size, sizeof(*nr_patches_all));
  MPI_Gather(&multi->nr_patches, 1, MPI_INT, nr_patches_all, 1, MPI_INT, 0, comm);
  MPI_Bcast(nr_patches_all, domain->size, MPI_INT, 0, comm);
  
  //Adjust for differences between the requested amount of patches and the actual amount
  //Multi might happen because the load-balancer doesn't know about now-to-be-created patches
  int nr_gpatches = 0;
  mprintf("nr_patches(from lb): %d\n", multi->nr_patches);
  if(domain->rank == 0)	{
    printf("\npatches on host: \n");
    for(int i = 0; i < domain->size; ++i) {
      nr_gpatches += nr_patches_all[i];
      mprintf("\t%d: %d\n", i, nr_patches_all[i]);
    }
    assert(nr_gpatches == multi->nr_gpatches);
  }

  for(int d = 0; d < 3; d++) {
    multi->_ldims[d] = multi->gdims[d] / multi->np[d];
  }
  
  //Create list of patches
  sfc_setup(&multi->sfc, multi->np);
  mrc_domain_dynamic_setup_patches(domain, nr_patches_all);
}

static void
mrc_domain_dynamic_destroy(struct mrc_domain *domain)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);
  
  free(multi->gpatch_off_all);
  free(multi->patches);
  bitfield3d_destroy(&multi->activepatches);
  map_destroy(domain);
}

static struct mrc_patch *
mrc_domain_dynamic_get_patches(struct mrc_domain *domain, int *nr_patches)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);
  if (nr_patches) {
    *nr_patches = multi->nr_patches;
  }
  return multi->patches;
}

static void
mrc_domain_dynamic_get_global_dims(struct mrc_domain *domain, int *dims)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);

  for (int d = 0; d < 3; d++) {
    dims[d] = multi->gdims[d];
  }
}

static void
mrc_domain_dynamic_get_nr_procs(struct mrc_domain *domain, int *nr_procs)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);

  for (int d = 0; d < 3; d++) {
    nr_procs[d] = multi->np[d];
  }
}

static void
mrc_domain_dynamic_get_bc(struct mrc_domain *domain, int *bc)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);

  for (int d = 0; d < 3; d++) {
    bc[d] = multi->bc[d];
  }
}

static void
mrc_domain_dynamic_get_nr_global_patches(struct mrc_domain *domain, int *nr_global_patches)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);

  *nr_global_patches = multi->nr_gpatches;
}

static void
mrc_domain_dynamic_get_idx3_patch_info(struct mrc_domain *domain, int idx[3],
				     struct mrc_patch_info *info)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);
  //Check if the patch is active
  if (!bitfield3d_isset(&multi->activepatches, idx)) {
    info->rank = -1;
    info->patch = -1;
    info->global_patch = -1;
    for(int d = 0; d < 3; d++) {
      info->ldims[d] = multi->_ldims[d];
      info->off[d] = idx[d] * multi->_ldims[d];
      info->idx3[d] = idx[d];
    }
    return;
  }

  int sfc_idx = sfc_idx3_to_idx(&multi->sfc, idx);
  int gpatch = map_sfc_idx_to_gpatch(domain, sfc_idx);
  mrc_domain_dynamic_get_global_patch_info(domain, gpatch, info);
}

static void
mrc_domain_dynamic_write(struct mrc_domain *domain, struct mrc_io *io)
{
  struct mrc_domain_dynamic* multi = mrc_domain_dynamic(domain);
  int nr_global_patches;
  mrc_domain_dynamic_get_nr_global_patches(domain, &nr_global_patches);
  mrc_io_write_attr_int(io, mrc_domain_name(domain), "nr_global_patches", multi->nr_gpatches);
  
  //Iterate over all global patches
  for (int i = 0; i < nr_global_patches; i++) {
    char path[strlen(mrc_domain_name(domain)) + 10];
    sprintf(path, "%s/p%d", mrc_domain_name(domain), i);
    mrc_io_write_attr_int3(io, path, "ldims", multi->_ldims);
    int sfc_idx = map_gpatch_to_sfc_idx(domain, i);
    mrc_io_write_attr_int(io, path, "sfc_idx", sfc_idx);
    int idx[3];
    sfc_idx_to_idx3(&multi->sfc, sfc_idx, idx);
    for(int d=0; d<3; ++d) idx[d] *= multi->_ldims[d];
    mrc_io_write_attr_int3(io, path, "off", idx);
  }
}

static void
mrc_domain_dynamic_plot(struct mrc_domain *domain)
{
}

static struct mrc_ddc *
mrc_domain_dynamic_create_ddc(struct mrc_domain *domain)
{
  struct mrc_ddc *ddc = mrc_ddc_create(domain->obj.comm);
  mrc_ddc_set_type(ddc, "multi");
  mrc_ddc_set_domain(ddc, domain);
  return ddc;
}

static struct mrc_param_select bc_descr[] = {
  { .val = BC_NONE       , .str = "none"     },
  { .val = BC_PERIODIC   , .str = "periodic" },
  {},
};

static struct mrc_param_select curve_descr[] = {
  { .val = CURVE_BYDIM   , .str = "bydim"    },
  { .val = CURVE_MORTON  , .str = "morton"   },
  { .val = CURVE_HILBERT , .str = "hilbert"  },
  {},
};

#define VAR(x) (void *)offsetof(struct mrc_domain_dynamic, x)
static struct param mrc_domain_dynamic_params_descr[] = {
  { "m"               , VAR(gdims)           , PARAM_INT3(32, 32, 32) },
  { "np"              , VAR(np)              , PARAM_INT3(1, 1, 1)    },
  { "bcx"             , VAR(bc[0])           , PARAM_SELECT(BC_NONE,
							    bc_descr) },
  { "bcy"             , VAR(bc[1])           , PARAM_SELECT(BC_NONE,
							    bc_descr) },
  { "bcz"             , VAR(bc[2])           , PARAM_SELECT(BC_NONE,
							    bc_descr) },
  { "curve_type"      , VAR(sfc.curve_type)  , PARAM_SELECT(CURVE_BYDIM,
							    curve_descr) },
  { "nr_patches"      , VAR(nr_patches)      , PARAM_INT(-1) },
  { "activepatches"   , VAR(p_activepatches) , PARAM_PTR(NULL) },
  {},
};
#undef VAR

struct mrc_domain_ops mrc_domain_dynamic_ops = {
  .name                  = "dynamic",
  .size                  = sizeof(struct mrc_domain_dynamic),
  .param_descr           = mrc_domain_dynamic_params_descr,
  .setup                 = mrc_domain_dynamic_setup,
  .view                  = mrc_domain_dynamic_view,
  .write                 = mrc_domain_dynamic_write,
  .destroy               = mrc_domain_dynamic_destroy,
  .get_patches           = mrc_domain_dynamic_get_patches,
  .get_global_dims       = mrc_domain_dynamic_get_global_dims,
  .get_nr_procs          = mrc_domain_dynamic_get_nr_procs,
  .get_bc                = mrc_domain_dynamic_get_bc,
  .get_nr_global_patches = mrc_domain_dynamic_get_nr_global_patches,
  .get_global_patch_info = mrc_domain_dynamic_get_global_patch_info,
  .get_local_patch_info  = mrc_domain_dynamic_get_local_patch_info,
  .get_idx3_patch_info   = mrc_domain_dynamic_get_idx3_patch_info,
  .plot                  = mrc_domain_dynamic_plot,
  .create_ddc            = mrc_domain_dynamic_create_ddc,
};
