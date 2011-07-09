
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

// ----------------------------------------------------------------------
// bydim space filling curve

static void
sfc_bydim_setup(struct mrc_domain_dynamic *multi)
{
}

static int
sfc_bydim_idx3_to_idx(struct mrc_domain_dynamic *multi, const int p[3])
{
  int *np = multi->np;
  return (p[2] * np[1] + p[1]) * np[0] + p[0];
}

static void
sfc_bydim_idx_to_idx3(struct mrc_domain_dynamic *multi, int idx, int p[3])
{
  int *np = multi->np;
  p[0] = idx % np[0]; idx /= np[0];
  p[1] = idx % np[1]; idx /= np[1];
  p[2] = idx;
}

// ----------------------------------------------------------------------
// morton space filling curve

static void
sfc_morton_setup(struct mrc_domain_dynamic *multi)
{
  int *np = multi->np;
  int *nbits = multi->nbits;
  for (int d = 0; d < 3; d++) {
    int n = np[d];
    nbits[d] = 0;
    while (n > 1) {
      n >>= 1;
      nbits[d]++;
    }
    // each dim must be power of 2
    assert(np[d] == 1 << nbits[d]);
  }

  multi->nbits_max = nbits[0];
  if (nbits[1] > multi->nbits_max) multi->nbits_max = nbits[1];
  if (nbits[2] > multi->nbits_max) multi->nbits_max = nbits[2];
}

static int
sfc_morton_idx3_to_idx(struct mrc_domain_dynamic *multi, const int p[3])
{
  int *nbits = multi->nbits;
  int nbits_max = multi->nbits_max;

  int pos = 0;
  int idx = 0;
  for (int b = 0; b < nbits_max; b++) {
    for (int d = 0; d < 3; d++) {
      if (b >= nbits[d])
	continue;

      if (p[d] & (1 << b)) {
	idx |= (1 << pos);
      }
      pos++;
    }
  }

  return idx;
}

static void
sfc_morton_idx_to_idx3(struct mrc_domain_dynamic *multi, int idx, int p[3])
{
  int *nbits = multi->nbits;
  int nbits_max = multi->nbits_max;

  for (int d = 0; d < 3; d++) {
    p[d] = 0;
  }

  int pos = 0;
  for (int b = 0; b < nbits_max; b++) {
    for (int d = 0; d < 3; d++) {
      if (b >= nbits[d])
	continue;

      if (idx & (1 << pos)) {
	p[d] |= (1 << b);
      }
      pos++;
    }
  }
}

// ----------------------------------------------------------------------
// hilbert space filling curve

#include "hilbert.h"

static void
sfc_hilbert_setup(struct mrc_domain_dynamic *multi)
{
  int *np = multi->np;
  int *nbits = multi->nbits;
  for (int d = 0; d < 3; d++) {
    int n = np[d];
    nbits[d] = 0;
    while (n > 1) {
      n >>= 1;
      nbits[d]++;
    }
    // each dim must be power of 2
    assert(np[d] == 1 << nbits[d]);
  }

  multi->nbits_max = nbits[0];
  if (nbits[1] > multi->nbits_max) multi->nbits_max = nbits[1];
  if (nbits[2] > multi->nbits_max) multi->nbits_max = nbits[2];
  multi->hilbert_nr_dims = 0;
  for (int d = 0; d < 3; d++) {
    if (nbits[d] == 0)
      continue;

    multi->hilbert_dim[multi->hilbert_nr_dims] = d;
    multi->hilbert_nr_dims++;
  }
  // all not invariant dimensions must be equal
  int d0 = multi->hilbert_dim[0];
  for (int i = 0; i < multi->hilbert_nr_dims; i++) {
    int d = multi->hilbert_dim[i];
    assert(nbits[d] == nbits[d0]);
  }
}

static int
sfc_hilbert_idx3_to_idx(struct mrc_domain_dynamic *multi, const int p[3])
{
  if (multi->hilbert_nr_dims == 0)
    return 0;

  int nbits_max = multi->nbits_max;

  bitmask_t p_bm[3];
  for (int i = 0; i < multi->hilbert_nr_dims; i++) {
    int d = multi->hilbert_dim[i];
    p_bm[i] = p[d];
  }
  return hilbert_c2i(multi->hilbert_nr_dims, nbits_max, p_bm);
}

static void
sfc_hilbert_idx_to_idx3(struct mrc_domain_dynamic *multi, int idx, int p[3])
{
  int nbits_max = multi->nbits_max;

  bitmask_t p_bm[3];
  hilbert_i2c(multi->hilbert_nr_dims, nbits_max, idx, p_bm);
  for (int d = 0; d < 3; d++) {
    p[d] = 0;
  }
  for (int i = 0; i < multi->hilbert_nr_dims; i++) {
    int d = multi->hilbert_dim[i];
    p[d] = p_bm[i];
  }
}

// ----------------------------------------------------------------------
// space filling curve

static void
sfc_setup(struct mrc_domain_dynamic *multi)
{
  switch (multi->curve_type) {
  case CURVE_BYDIM: return sfc_bydim_setup(multi);
  case CURVE_MORTON: return sfc_morton_setup(multi);
  case CURVE_HILBERT: return sfc_hilbert_setup(multi);
  default: assert(0);
  }
}

static int
sfc_idx3_to_idx(struct mrc_domain_dynamic *multi, const int p[3])
{
  switch (multi->curve_type) {
  case CURVE_BYDIM: return sfc_bydim_idx3_to_idx(multi, p);
  case CURVE_MORTON: return sfc_morton_idx3_to_idx(multi, p);
  case CURVE_HILBERT: return sfc_hilbert_idx3_to_idx(multi, p);
  default: assert(0);
  }
}

static void
sfc_idx_to_idx3(struct mrc_domain_dynamic *multi, int idx, int p[3])
{
  switch (multi->curve_type) {
  case CURVE_BYDIM: return sfc_bydim_idx_to_idx3(multi, idx, p);
  case CURVE_MORTON: return sfc_morton_idx_to_idx3(multi, idx, p);
  case CURVE_HILBERT: return sfc_hilbert_idx_to_idx3(multi, idx, p);
  default: assert(0);
  }
}

// ======================================================================
// map
// 
// maps between global patch index (contiguous) and 1D SFC idx
// (potentially non-contiguous)

static void
map_create(struct mrc_domain *domain, int *sfc_indices, int nr_gpatches)
{
  struct mrc_domain_dynamic *this = mrc_domain_dynamic(domain);

  this->gp = malloc(sizeof(int) * this->nr_gpatches);
  for (int i = 0; i < nr_gpatches; i++) {
    this->gp[i] = sfc_indices[i];
  }

  //Create the bintree for performant searching
  int vals[nr_gpatches];
  for (int i = 0; i < nr_gpatches; i++) {
    vals[i] = i;
  }
  bintree_create_from_ordered_list(&this->g_patches, sfc_indices, vals, nr_gpatches);
}

static void
map_destroy(struct mrc_domain *domain)
{
  struct mrc_domain_dynamic *this = mrc_domain_dynamic(domain);

  free(this->gp);
  bintree_destroy(&this->g_patches);
}

static int
map_sfc_idx_to_gpatch(struct mrc_domain *domain, int sfc_idx)
{
  struct mrc_domain_dynamic *this = mrc_domain_dynamic(domain);

  int retval;
  int rc = bintree_get(&this->g_patches, sfc_idx, &retval);
  if (rc == 0) {
    return -1;
  }
  return retval;
}

static int
map_gpatch_to_sfc_idx(struct mrc_domain *domain, int gpatch)
{
  struct mrc_domain_dynamic *this = mrc_domain_dynamic(domain);

  return this->gp[gpatch];
}

// ======================================================================

static void
sfc_idx_to_rank_patch(struct mrc_domain *domain, int sfc_idx,
		      int *rank, int *patch)
{
  struct mrc_domain_dynamic *this = mrc_domain_dynamic(domain);
  //Get gp->array index
  int gpatch = map_sfc_idx_to_gpatch(domain, sfc_idx);
  if (gpatch < 0) {
    *rank = -1;
    *patch = -1;
    return;
  }
  
  *rank = this->rank[gpatch];
  *patch = this->patch[gpatch];
}

// ======================================================================

static void
mrc_domain_dynamic_view(struct mrc_domain *domain)
{
}

// FIXME, get rid of this one always use idx3 one?
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
  sfc_idx_to_idx3(multi, sfc_idx, p3);
  for (int d = 0; d < 3; d++) {
    info->ldims[d] = multi->ldims[d];
    info->off[d] = p3[d] * multi->ldims[d];
    info->idx3[d] = p3[d];
  }
}

static void
mrc_domain_dynamic_get_local_patch_info(struct mrc_domain *domain, int patch,
				      struct mrc_patch_info *info)
{
  struct mrc_domain_dynamic *multi = mrc_domain_dynamic(domain);

  int sfc_idx = multi->gpatch[patch]; // FIXME, there must be an easier way
  int gpatch = map_sfc_idx_to_gpatch(domain, sfc_idx);
  mrc_domain_dynamic_get_global_patch_info(domain, gpatch, info);
}

static void mrc_domain_dynamic_setup_patches(struct mrc_domain *domain, int firstpatch)
{
  struct mrc_domain_dynamic *this = mrc_domain_dynamic(domain);
  
  int ngp = this->np[0] * this->np[1] * this->np[2];
  
  int gpatchkeys[this->nr_gpatches];
  
  this->rank = malloc(sizeof(int) * this->nr_gpatches);
  this->patch = malloc(sizeof(int) * this->nr_gpatches);
  this->gpatch = malloc(sizeof(int) * this->nr_patches);
  this->patches = malloc(sizeof(*this->patches) * this->nr_patches);
  
  //Gather firstpatch from all ranks
  int *firstpatch_all = calloc(domain->size, sizeof(*firstpatch_all));
  MPI_Comm comm = mrc_domain_comm(domain);
  MPI_Gather(&firstpatch, 1, MPI_INT, firstpatch_all, 1, MPI_INT, 0, comm);
  MPI_Bcast(firstpatch_all, domain->size, MPI_INT, 0, comm);
 
  int activerank = 0;
  int npatches = 0;
  
  //assert(this->nr_gpatches == bitfield3d_count_bits_set(&this->activepatches));
  
  //TODO Find a smarter way than iterating over all possible patches
  for(int i=0; i<ngp; ++i)
  {
    int idx[3];
    sfc_idx_to_idx3(this, i, idx);
    if(bitfield3d_isset(&this->activepatches, idx))
    {
      //Calculate rank
      if(activerank < ( domain->size - 1 ) && firstpatch_all[activerank+1] <= npatches) activerank++;
      
      //Register the patch
      gpatchkeys[npatches] = i;

      this->rank[npatches] = activerank;
      int lpatch = npatches - firstpatch_all[activerank];
      this->patch[npatches] = lpatch;
      
      if(activerank == domain->rank) //Create the patch on this processor
	  {
		  this->gpatch[lpatch] = i;
		  
		  //Setup patches[lpatch]
		  for(int d=0; d<3; ++d)
		  {
			  this->patches[lpatch].off[d] = idx[d] * this->ldims[d];
			  this->patches[lpatch].ldims[d] = this->ldims[d];
		  }
      }
      npatches+=1;
    }
  }
  
  map_create(domain, gpatchkeys, this->nr_gpatches);
}

static void
mrc_domain_dynamic_setup(struct mrc_domain *domain)
{
  assert(!domain->is_setup);
  domain->is_setup = true;

  struct mrc_domain_dynamic *this = mrc_domain_dynamic(domain);
  	
  //Assert some properties we rely on
  for(int d=0; d<3; ++d) assert((this->gdims[d] % this->np[d]) == 0);	//The domain-size must be a mutliple of the patch-size
  
  //Copy the activepatch-list
  bitfield3d_copy(&this->activepatches, this->p_activepatches);
  
  this->nr_gpatches = bitfield3d_count_bits_set(&this->activepatches);
  
  MPI_Comm comm = mrc_domain_comm(domain);
  MPI_Comm_rank(domain->obj.comm, &domain->rank);
  MPI_Comm_size(domain->obj.comm, &domain->size);
  
  //If nr_patches is not set, try to distribute patches evenly
  if(this->nr_patches < 0)
  {
    this->nr_patches = this->nr_gpatches / domain->size + (domain->rank < (this->nr_gpatches % domain->size) ? 1 : 0);
  }
  
  int *nr_patches_all = calloc(domain->size, sizeof(*nr_patches_all));
  MPI_Gather(&this->nr_patches, 1, MPI_INT, nr_patches_all, 1, MPI_INT, 0, comm);
  MPI_Bcast(nr_patches_all, domain->size, MPI_INT, 0, comm);
  
  //Adjust for differences between the requested amount of patches and the actual amount
  //This might happen because the load-balancer doesn't know about now-to-be-created patches
  int nrequestedpatches = 0;
	printf("Rank %d, nr_patches(from lb): %d\n", domain->rank, this->nr_patches);
	if(domain->rank == 0)
	{
		printf("\npatches on host: \n");
		for(int i=0; i<domain->size; ++i)
		{
			nrequestedpatches += nr_patches_all[i];
			printf("\t%d: %d\n", i, nr_patches_all[i]);
		}
	}
	
	/* This is not neccessary anymore after modifying the balancer
	 * I Left it here for now as I'm not convinced yet this was the best way
	
	int excess = (this->nr_gpatches - nrequestedpatches);
	
	if(excess != 0)
	{
		assert(0); //This should not happen anymore
		this->nr_patches += (excess / domain->size) + (excess % domain->size > domain->rank ? 1 : 0);
	}
	else if(excess < 0)
	{
		this->nr_patches += (excess / domain->size) + (domain->rank < (-excess % domain->size) ? -1 : 0);
	}
	
	//Gather new # of patches
	MPI_Gather(&this->nr_patches, 1, MPI_INT, nr_patches_all, 1, MPI_INT, 0, comm);
	MPI_Bcast(nr_patches_all, domain->size, MPI_INT, 0, comm);
	*/
	
	//and calculate first patch to process on this domain
	int firstpatch = 0;	//The first global patch to be owned by THIS processor
	for(int i=0; i<domain->rank; ++i)
	{
		firstpatch += nr_patches_all[i];
	}
	
	for(int d=0; d<3; ++d)
	{
		this->ldims[d] = this->gdims[d] / this->np[d];
	}
  
	//Create list of patches
	sfc_setup(this);
	mrc_domain_dynamic_setup_patches(domain, firstpatch);
}

static void
mrc_domain_dynamic_destroy(struct mrc_domain *domain)
{
  struct mrc_domain_dynamic *this = mrc_domain_dynamic(domain);
  
  free(this->patches);
  
  free(this->rank);
  free(this->patch);
  free(this->gpatch);
  //
  bitfield3d_destroy(&this->activepatches);
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
  struct mrc_domain_dynamic *this = mrc_domain_dynamic(domain);
  //Check if the patch is active
  if(!bitfield3d_isset(&this->activepatches, idx))
  {
    info->rank = -1;
    info->patch = -1;
    info->global_patch = -1;
    for(int d=0; d<3; ++d)
    {
      info->ldims[d] = this->ldims[d];
      info->off[d] = idx[d] * this->ldims[d];
      info->idx3[d] = idx[d];
    }
    return;
  }

  int sfc_idx = sfc_idx3_to_idx(this, idx);
  int gpatch = map_sfc_idx_to_gpatch(domain, sfc_idx);
  mrc_domain_dynamic_get_global_patch_info(domain, gpatch, info);
}

static void
mrc_domain_dynamic_write(struct mrc_domain *domain, struct mrc_io *io)
{
  struct mrc_domain_dynamic* this = mrc_domain_dynamic(domain);
  int nr_global_patches;
  mrc_domain_dynamic_get_nr_global_patches(domain, &nr_global_patches);
  mrc_io_write_attr_int(io, mrc_domain_name(domain), "nr_global_patches", this->nr_gpatches);
  
  //Iterate over all global patches
  for (int i = 0; i < nr_global_patches; i++) {
    char path[strlen(mrc_domain_name(domain)) + 10];
    sprintf(path, "%s/p%d", mrc_domain_name(domain), i);
    mrc_io_write_attr_int3(io, path, "ldims", this->ldims);
    int sfc_idx = map_gpatch_to_sfc_idx(domain, i);
    mrc_io_write_attr_int(io, path, "sfc_idx", sfc_idx);
    int idx[3];
    sfc_idx_to_idx3(this, sfc_idx, idx);
    for(int d=0; d<3; ++d) idx[d] *= this->ldims[d];
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
  { "curve_type"      , VAR(curve_type)      , PARAM_SELECT(CURVE_BYDIM,
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
