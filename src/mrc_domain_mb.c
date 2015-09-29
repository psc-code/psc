
#include <mrc_domain_private.h>
#include <mrc_fld.h>
#include <mrc_block_factory.h>
#include <mrc_crds.h>
#include <mrc_params.h>
#include <mrc_ddc.h>
#include <mrc_io.h>

#include <stdlib.h>
#include <string.h>

#define mrc_domain_mb(domain) mrc_to_subobj(domain, struct mrc_domain_mb)
#define CE assert(ierr == 0)

const char *face2str[NR_FACES] = {
  "FACE_LEFT",
  "FACE_RIGHT",
  "FACE_BOTTOM",
  "FACE_TOP",
  "FACE_BACK",
  "FACE_FRONT",
};

// ----------------------------------------------------------------------
// MB_Create


static void
_mrc_domain_mb_create(struct mrc_domain *mb)
{
  mrc_crds_set_type(mb->crds, "mb");
  // I would prefer this to be marked as a member object in the params
  // struct, but there doesn't seem to be any easy way to recover
  // via the params interface. So we do it this hacky way.
  struct mrc_domain_mb *sub = mrc_domain_mb(mb);
  sub->_block_factory = mrc_block_factory_create(mrc_domain_comm(mb));
  mrc_domain_add_child(mb, (struct mrc_obj *) sub->_block_factory);
  
}

static void
_mrc_domain_mb_destroy(struct mrc_domain *mb)
{

  struct mrc_domain_mb *sub = mrc_domain_mb(mb);

  // free the blocks created by the factory.
  free(sub->mb_blocks);
  free(sub->patches);
  free(sub->patch_info);

}  

static void
_mrc_domain_mb_view(struct mrc_domain *domain)
{
  struct mrc_domain_mb *sub = mrc_domain_mb(domain);
}

// ----------------------------------------------------------------------
// MB_Setup

static int
block_index_to_patch_index(struct mrc_domain *mb, int b, int k[3])
{
  struct mrc_domain_mb *sub = mrc_domain_mb(mb);
  return b * sub->mb_ppb_total + 
    (k[2]*sub->ppb[1] + k[1])*sub->ppb[0] + k[0];
}

static void
mrc_domain_mb_get_blocks(struct mrc_domain *domain, struct MB_block **blocks, int *nr_blocks)
{
  struct mrc_domain_mb *sub = mrc_domain_mb(domain);
  *nr_blocks = sub->nr_blocks;
  *blocks = sub->mb_blocks;
}

static void
mrc_domain_mb_get_global_patch_info(struct mrc_domain *domain, int patch,
					struct mrc_patch_info *info)
{
  struct mrc_domain_mb *mb = mrc_domain_mb(domain);
  *info = mb->patch_info[patch];
  for (int d = 0; d < 3; d++) {
    info->ldims[d] = mb->patches[patch].ldims[d];
    info->off[d] = mb->patches[patch].off[d];
  }
  info->patch = info->global_patch - mb->gpatch_off;
}


static void
mrc_domain_mb_get_local_patch_info(struct mrc_domain *domain, int patch,
				      struct mrc_patch_info *info)
{
  struct mrc_domain_mb *mb = mrc_domain_mb(domain);

  mrc_domain_mb_get_global_patch_info(domain, mb->gpatch_off + patch,
					 info);
}


// ----------------------------------------------------------------------
// subclass setup

static void
mrc_domain_mb_do_setup(struct mrc_domain *mb)
{
  MPI_Comm comm = mrc_domain_comm(mb);
  MPI_Comm_size(comm, &mb->size);
  MPI_Comm_rank(comm, &mb->rank);
  struct mrc_domain_mb *sub = mrc_domain_mb(mb);

  // Blocks may already exist (if we read the domain from a file),
  // otherwise we need to run a factory.
  if ( !sub->mb_blocks ) {
    mrc_block_factory_run(sub->_block_factory, mb);
  }
  // global list of all patches and related info
  sub->mb_ppb_total = sub->ppb[0]*sub->ppb[1]*sub->ppb[2];
  sub->nr_global_patches = sub->nr_blocks * sub->mb_ppb_total;
  sub->patches = (struct mrc_patch *) calloc(sub->nr_global_patches, sizeof(struct mrc_patch));
  sub->patch_info = (struct mrc_patch_info *) calloc(sub->nr_global_patches, sizeof(struct mrc_patch_info));

  // list of patches which will be on the local processor
  int rank = mb->rank, size = mb->size;
  sub->nr_patches = sub->nr_global_patches / size + (rank < sub->nr_global_patches % size);
  int start = 0, j = 0;
  for (int i = 0; i < size; i++) {
    int nr_loc_patches = sub->nr_global_patches / size + (i < sub->nr_global_patches % size);
    for (int p = 0; p < nr_loc_patches; p++) {
      sub->patch_info[j + p].rank = i;
    }
    if (i == rank) {
      sub->gpatch_off = j;
      start = j;
    }
    j += nr_loc_patches;
  }

  int k = 0;
  MB_foreach_block(mb, block) {
#if 0
    ASSERT(block->mx[0] % mb->mb_ppb[0] == 0);
    ASSERT(block->mx[1] % mb->mb_ppb[1] == 0);
    ASSERT(block->mx[2] % mb->mb_ppb[2] == 0);
#else
    assert(block->mx[0] % sub->ppb[0] == 0);
    assert(block->mx[1] % sub->ppb[1] == 0);
    assert(block->mx[2] % sub->ppb[2] == 0);
#endif
    int kk[3];
    for (kk[2] = 0; kk[2] < sub->ppb[2]; kk[2]++) {
      for (kk[1] = 0; kk[1] < sub->ppb[1]; kk[1]++) {
	for (kk[0] = 0; kk[0] < sub->ppb[0]; kk[0]++) {
	  struct mrc_patch *patch = &sub->patches[k];
	  struct mrc_patch_info *info = &sub->patch_info[k];
	  info->global_patch = k;
	  info->p_block = block->nr_block;
	  for (int d = 0; d < 3; d++) {
	    patch->ldims[d] = block->mx[d] / sub->ppb[d];
	    info->p_ix[d]  = kk[d] * patch->ldims[d];
	    // FIXME: This just lines all the blocks up in X !!
	    patch->off[d] = info->p_ix[d];
	    if (d == 0) {
	      for (int b = block->nr_block - 1; b >= 0; b--) {
		patch->off[d] += sub->mb_blocks[b].mx[d];
	      }
	    }
	    info->p_bmx[d] = block->mx[d];
	  }

	  for (int f = 0; f < NR_FACES; f++) {
	    int dir = face2dir(f), bnd = face2bnd(f);
	    if ((bnd == 0 && kk[dir] == 0) ||
		(bnd == 1 && kk[dir] == sub->ppb[dir]-1)) {
	      // use b.c. from block
	      info->p_pface[f].pf_face  = block->faces[f].face;
	      info->p_pface[f].pf_btype = block->faces[f].btype;
	      int nblock = block->faces[f].block;
	      if (nblock < 0) {
		info->p_pface[f].pf_patch = nblock;
		continue;
	      }
	      for (int d = 0; d < 3; d++) {
		info->p_pface[f].pf_map[d] = block->faces[f].map[d];
	      }
	      int nface = block->faces[f].face;
	      int *map = block->faces[f].map;
	      int ndir = face2dir(nface), nbnd = face2bnd(nface); 
	      int kb[3];
	      kb[ndir] = nbnd ? sub->ppb[ndir] - 1 : 0;
	      for (int d = 0; d < 3; d++) {
		if (map[d] != 0) {
		  if (map[d] & MB_R) {
		    kb[map2dir(map[d])] = sub->ppb[d]-1-kk[d];
		  } else {
		    kb[map2dir(map[d])] = kk[d];
		  }
		}
	      }
	      info->p_pface[f].pf_patch = block_index_to_patch_index(mb, nblock, kb);
	    } else {
	      // create new internal b.c.
	      int kb[3];
	      for (int d = 0; d < 3; d++) {
		kb[d] = kk[d];
		info->p_pface[f].pf_map[d] = d + 1;
	      }
	      kb[dir] += bnd ? 1 : -1;
	      info->p_pface[f].pf_face = f ^ 1;
	      info->p_pface[f].pf_map[dir] = 0;
	      info->p_pface[f].pf_patch = block_index_to_patch_index(mb, block->nr_block, kb );
	      info->p_pface[f].pf_btype = BTYPE_PATCH; 	
	    }
	  }
	  k++;
	}
      }
    }
  } MB_foreach_block_end;
  
};

static void
_mrc_domain_mb_setup(struct mrc_domain *mb) 
{

  mrc_domain_mb_do_setup(mb);

  if (mb->ddc) {
    mrc_ddc_set_type(mb->ddc, "mb");
    mrc_ddc_set_domain(mb->ddc, mb);
  }

  mrc_domain_setup_super(mb);
 
}

static void
mrc_domain_mb_get_global_dims(struct mrc_domain *domain, int *dims)
{
  // This isn't well defined for multi-block domains
  assert(0);

}

static void
mrc_domain_mb_get_nr_procs(struct mrc_domain *domain, int *nr_procs)
{
  // FIXME: the 3d nr procs is meaningless for multiblock domains,
  // so this is just a hacky work around so fld writers don't 
  // choke
  struct mrc_domain_mb *sub = mrc_domain_mb(domain);
  nr_procs[0] = sub->nr_blocks * sub->ppb[0];
  for (int d = 1; d < 3; d++) {
    nr_procs[d] = sub->ppb[d];
  }
}


static struct mrc_patch *
mrc_domain_mb_get_patches(struct mrc_domain *domain, int *nr_patches)
{
  struct mrc_domain_mb *mb = mrc_domain_mb(domain);
  if (nr_patches) {
    *nr_patches = mb->nr_patches;
  }
  return &mb->patches[mb->gpatch_off];
}

static void
mrc_domain_mb_get_nr_global_patches(struct mrc_domain *domain, int *nr_global_patches)
{
  struct mrc_domain_mb *mb = mrc_domain_mb(domain);

  *nr_global_patches = mb->nr_global_patches;
}

static struct mrc_ddc *
mrc_domain_mb_create_ddc(struct mrc_domain *domain)
{
  struct mrc_ddc *ddc = mrc_ddc_create(domain->obj.comm);
  mrc_ddc_set_type(ddc, "mb");
  mrc_ddc_set_domain(ddc, domain);
  return ddc;
}

static void
mrc_domain_mb_write(struct mrc_domain *domain, struct mrc_io *io)
{

  struct mrc_domain_mb *sub = mrc_domain_mb(domain);
  mrc_io_write_int(io, domain, "nr_blocks", sub->nr_blocks);

  const char *path = mrc_io_obj_path(io, domain);  
  char path2[strlen(path) + 10];

  for (int b = 0; b < sub->nr_blocks; b++) {
    sprintf(path2, "%s/block-%d", path, b);
    struct MB_block *block = &sub->mb_blocks[b];
    mrc_io_write_attr_int3(io, path2, "mx", block->mx);
    mrc_io_write_attr_int(io, path2, "nr_block", block->nr_block);
    mrc_io_write_attr_double3(io, path2, "xl", block->xl);
    mrc_io_write_attr_double3(io, path2, "xh", block->xh);
    // FIXME: write coord_gens?

    for (int f = 0; f < NR_FACES; f++) {
      struct MB_face *face = &block->faces[f];
      char face_path[strlen(path2) + 10];
      sprintf(face_path, "%s/face-%d", path2, f);
      mrc_io_write_attr_int(io, face_path, "block", face->block);
      mrc_io_write_attr_int(io, face_path, "face", face->face);
      mrc_io_write_attr_int3(io, face_path, "map", face->map);
      mrc_io_write_attr_int(io, face_path, "btype", face->btype);
    }
  }
}


static void
mrc_domain_mb_read(struct mrc_domain *domain, struct mrc_io *io)
{
#warning MB initialization orders are all sorts of messed up here
  struct mrc_domain_mb *sub = mrc_domain_mb(domain);
  mrc_io_read_int(io, domain, "nr_blocks", &sub->nr_blocks);
  sub->mb_blocks = (struct MB_block *) calloc(sub->nr_blocks, sizeof(struct MB_block));

  const char *path = mrc_io_obj_path(io, domain);  
  char path2[strlen(path) + 10];
  
  for (int b = 0; b < sub->nr_blocks; b++) {
    sprintf(path2, "%s/block-%d", path, b);
    struct MB_block *block = &sub->mb_blocks[b];
    mrc_io_read_attr_int3(io, path2, "mx", &block->mx);
    mrc_io_read_attr_int(io, path2, "nr_block", &block->nr_block);
    mrc_io_read_attr_double3(io, path2, "xl", &block->xl);
    mrc_io_read_attr_double3(io, path2, "xh", &block->xh);
    
    for (int f = 0; f < NR_FACES; f++) {
      struct MB_face *face = &block->faces[f];
      char face_path[strlen(path2) + 10];
      sprintf(face_path, "%s/face-%d", path2, f);
      mrc_io_read_attr_int(io, face_path, "block", &face->block);
      mrc_io_read_attr_int(io, face_path, "face", &face->face);
      mrc_io_read_attr_int3(io, face_path, "map", &face->map);
      mrc_io_read_attr_int(io, face_path, "btype", &face->btype);
    }
  }


  mrc_domain_mb_do_setup(domain);
  mrc_domain_read_super(domain, io);
}


#define VAR(x) (void *)offsetof(struct mrc_domain_mb, x)
static struct param mrc_domain_mb_params_descr[] = {
  { "block_factory"   , VAR(_block_factory)   , PARAM_OBJ(mrc_block_factory) },
  { "ppb"             , VAR(ppb)             , PARAM_INT3(1,1,1),
    .help = " number of patches per block in x, y, z direction" },
  {},
};
#undef VAR


static struct mrc_obj_method mrc_domain_mb_methods[] = {
  MRC_OBJ_METHOD("block_idx_to_patch_idx", block_index_to_patch_index),
  MRC_OBJ_METHOD("get_blocks", mrc_domain_mb_get_blocks),
};

struct mrc_domain_ops mrc_domain_mb_ops = {
  .name                      = "mb",
  .methods                   = mrc_domain_mb_methods,
  .size                      = sizeof(struct mrc_domain_mb),
  .param_descr               = mrc_domain_mb_params_descr,
  .create                    = _mrc_domain_mb_create,
  .setup                     = _mrc_domain_mb_setup,
  .destroy                   = _mrc_domain_mb_destroy,
  .view                      = _mrc_domain_mb_view,
  .write                     = mrc_domain_mb_write,
  .read                      = mrc_domain_mb_read,
  .get_global_dims           = mrc_domain_mb_get_global_dims,
  .get_nr_procs              = mrc_domain_mb_get_nr_procs,
  .get_global_patch_info     = mrc_domain_mb_get_global_patch_info,
  .get_local_patch_info      = mrc_domain_mb_get_local_patch_info,
  .get_nr_global_patches     = mrc_domain_mb_get_nr_global_patches,
  .get_patches               = mrc_domain_mb_get_patches,
  .create_ddc                = mrc_domain_mb_create_ddc,
};



// M3d functions that need to be in this file for linking reasons
// will be deprecated
/* int */
/* M3d_WriteASCII2(struct mrc_fld_x *f, const char *filename, struct MB_bc *bc) */
/* { */
/*   int ierr; */
  
/*   pfb; */
/*   struct mrc_fld_x *l = mrc_fld_x_create(mrc_fld_x_comm(f)); */
/*   ierr = M3d_GetLocal(f, l); CE; */
/*   ierr = M3d_FillGhostCells(f, l, bc); CE; */
/*   ierr = M3d_WriteASCII(l, filename, f->_domain->mb_trafo->_cc); CE; */
/*   mrc_fld_x_destroy(l); */
/*   pfr; */
/* } */

// ======================================================================
// M3d_ReadH5

// FIXME, don't pass mb, find it (cache...)

/* int */
/* M3d_ReadH5(struct mrc_fld_x *x, hid_t loc, const char *name, struct mrc_domain *mb) */
/* { */
/*   struct mrc_domain_mb *sub = mrc_domain_mb(mb); */
/*   int ierr; */
/*   int bs, sw; */
  
/*   pfb; */
/*   hid_t g_id = MRC_H5Gopen(loc, name, H5P_DEFAULT); */
/*   ierr = MRC_H5Aread_int(g_id, "sw", &sw); CE; */
/*   ierr = MRC_H5Aread_int(g_id, "bs", &bs); CE; */
/*   ASSERT(bs > 0); */

/*   ierr = MB_GetFld(mb, bs, sw, x); */

/*   for (int b = 0; b < sub->nr_blocks; b++) { */
/*     struct MB_block *blk = &sub->mb_blocks[b]; */
/*     char s[20]; */
/*     sprintf(s, "block-%d", b); */
/*     hid_t b_id = 0, dataset = 0, filespace = 0; */
/*     if (mb->rank == 0) { */
/*       b_id = H5Gopen(g_id, s, H5P_DEFAULT); */

/*       dataset = H5Dopen(b_id, "v", H5P_DEFAULT); */
/*       ASSERT(dataset >= 0); */

/*       filespace = H5Dget_space(dataset); */
/*     } */

/*     // all patches that make up this block */
/*     for (int p = 0; p < sub->nr_global_patches; p++) { */
/*       struct mrc_patch *patch = &sub->patches[p]; */
/*       struct mrc_patch_info info; */
/*       mrc_domain_get_global_patch_info(mb, p, &info); */
/*       if (info.p_block != b) */
/* 	continue; */
      
/*       int bl[3] = {}, bh[3] = {}; */
/*       for (int d = 0; d < 3; d++) { */
/* 	if (info.p_ix[d] == 0) { */
/* 	  bl[d] = sw; */
/* 	}  */
/* 	if (info.p_ix[d] + patch->ldims[d] == blk->mx[d]) { */
/* 	  bh[d] = sw; */
/* 	} */
/*       } */
/*       hsize_t mdims[] = { patch->ldims[2] + 2*sw, */
/* 			  patch->ldims[1] + 2*sw, */
/* 			  patch->ldims[0] + 2*sw, */
/* 			  mrc_fld_m3d_nr_comps(x)}; */
/*       int msize = mdims[0] * mdims[1] * mdims[2] * mdims[3]; */

/*       if (mb->rank == 0) { */
/* 	double *data = NULL; */
/* 	ierr = PetscMalloc(msize * sizeof(double), &data); CE; */
/* 	hsize_t offset[] = { info.p_ix[2], */
/* 			     info.p_ix[1], */
/* 			     info.p_ix[0], */
/* 			     0 }; */
/* 	ierr = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, */
/* 				   mdims, NULL); CE; */
/* 	hid_t memspace = H5Screate_simple(4, mdims, NULL); */
	
/* 	hid_t plist = H5Pcreate(H5P_DATASET_XFER); */
/* 	//	ierr = H5Pset_dxpl_mpio(plist, H5FD_MPIO_INDEPENDENT); CE; */
/* 	ierr = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, plist, */
/* 		       data); CE; */
/* 	ierr = H5Pclose(plist); CE; */
	
/* 	if (info.rank == 0) {  */
/* 	  memcpy(&MRC_D5(x,0, -sw,-sw,-sw,info.patch), data, */
/* 		 msize * sizeof(double)); CE; */
/* 	} else { */
/* 	  ierr = MPI_Send(data, msize, MPI_DOUBLE, info.rank, 0x3000, */
/* 			  mrc_domain_comm(mb)); CE; */
/* 	} */
/* 	ierr = PetscFree(data); CE; */
/*       } else { */
/* 	if (info.rank == mb->rank) { */
/* 	  ierr = MPI_Recv(&MRC_D5(x,0, -sw,-sw,-sw,info.patch), msize, MPI_DOUBLE, */
/* 			  0, 0x3000, mrc_domain_comm(mb), MPI_STATUS_IGNORE); CE; */
/* 	} */
/*       } */
/*     } */

/*     if (mb->rank == 0) { */
/*       ierr = H5Dclose(dataset); CE; */
/*       ierr = H5Sclose(filespace); CE; */
/*       ierr = H5Gclose(b_id); CE; */
/*     } */
/*   } */

/*   ierr = MRC_H5Gclose(g_id); CE; */
/*   pfr; */
/* } */
 


// =====================================================================

/* // k is GLOBAL patch index */
/* int */
/* MB_writeASCII_patch(struct mrc_domain *mb, int k, struct mrc_fld_x *f, FILE *file, struct mrc_fld_x *coord) */
/* { */

/*   struct mrc_domain_mb *sub = mrc_domain_mb(mb); */
/*   struct mrc_patch *patch = &sub->patches[k]; */
/*   struct mrc_patch_info info; */
/*   mrc_domain_get_global_patch_info(mb, k, &info); */

/*   pfb; */
/*   int jxs, jys, jxe, jye; */
/*   int sw = mrc_fld_m3d_sw(f); */
/*   jxs = -sw; jxe = info.ldims[0] + sw; */
/*   jys = -sw; jye = info.ldims[1] + sw; */
/*   int jz = 0; */
  
/*   struct MB_coord *crds = mrc_domain_get_crds(mb); */
/*   for (int jy = jys; jy < jye; jy++) { */
/*     for (int jx = jxs; jx < jxe; jx++) { */
/*       int jjx = jx + info.p_ix[0] + sub->mb_blocks[info.p_block].lxs[0]; */
/*       int jjy = jy + info.p_ix[1] + sub->mb_blocks[info.p_block].lxs[1]; */
/*       fprintf(file, "%d %d %g %g", jjx, jjy, XI0(crds, jx, info.patch),  */
/* 	      XI1(crds, jy, info.patch)); */
/*       if (coord) { */
/* 	fprintf(file, " % g %g", MRC_D5(coord,0, jx,jy,jz,info.patch), */
/* 		MRC_D5(coord,1, jx,jy,jz),info.patch); */
/*       } */
/*       for (int m = 0; m < mrc_fld_m3d_nr_comps(f); m++) { */
/* 	fprintf(file, " %.15g", MRC_D5(f, m, jx,jy,jz),info.patch); */
/*       } */
/*       fprintf(file, "\n"); */
/*     } */
/*     fprintf(file, "\n"); */
/*   } */
/*   pfr; */
/* } */

/* int */
/* MB_WriteASCII(struct mrc_domain *mb, struct mrc_fld_x *f, const char *filename, struct mrc_fld_x *coord) */
/* { */
/*   int ierr; */
  
/*   pfb; */
/*   char *fn; */
/*   ierr = PetscMalloc(strlen(filename)+10, &fn); CE; */

/*   sprintf(fn, "%s.asc", filename); */
/*   printf("Writing '%s'\n", fn); */
/*   fflush(stdout); */
/*   FILE *file = fopen(fn, "w"); */
/*   MB_foreach_patch(mb, patch) { */
/*     MB_writeASCII_patch(mb, mrc_domain_mb(mb)->gpatch_off + patch, f, file, coord); */
/*     fprintf(file, "\n"); */
/*   } MB_foreach_patch_end; */
/*   fclose(file); */

/*   MB_foreach_patch(mb, patch) { */
/*     sprintf(fn, "%s-%d.asc", filename, mrc_domain_mb(mb)->gpatch_off + patch); */
/*     FILE *file = fopen(fn, "w"); */
/*     MB_writeASCII_patch(mb, mrc_domain_mb(mb)->gpatch_off + patch, f, file, coord); */
/*     fclose(file); */
/*   } MB_foreach_patch_end; */

/*   ierr = PetscFree(fn); CE; */
/*   pfr; */
/* } */


/* // ====================================================================== */
/* // MB_WriteH5 */

// ======================================================================
// MB_ReadH5


