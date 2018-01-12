
#include "mrc_ddc_private.h"

#include <mrc_params.h>
#include <mrc_domain_private.h>
#include <petsc.h>

#include <stdlib.h>
#include <assert.h>

#define mrc_domain_mb(domain) mrc_to_subobj(domain, struct mrc_domain_mb)
#define CE assert(ierr == 0)
#define pfb 
#define pfr return(0)

struct mrc_ddc_mb {
  struct mrc_domain *domain;
  VecScatter mb_gtol[MAX_SW][MAX_BS];
  VecScatter mb_ltol_edge[MAX_SW][MAX_BS];
  // FIXME: Patch per block should probably go in here
};

// Some typedefs to make the method fetching a little easier
typedef Vec (*fgp_t)(struct mrc_fld *);
typedef void (*fpp_t)(struct mrc_fld *, Vec *);
typedef void (*fsp_t)(struct mrc_fld *, Vec);

#define global_size(mb, sw) (((mb)->patches[0].ldims[0] + 2 * (sw[0])) *	\
			     ((mb)->patches[0].ldims[1] + 2 * (sw[1])) *	\
			     ((mb)->patches[0].ldims[2] + 2 * (sw[2])) *	\
			     (mb)->nr_global_patches)

#define local_size(mb, sw) (((mb)->patches[0].ldims[0] + 2 * (sw[0])) *  \
			    ((mb)->patches[0].ldims[1] + 2 * (sw[1])) *  \
			    ((mb)->patches[0].ldims[2] + 2 * (sw[2])) *  \
			    (mb)->nr_patches)

#define to_mrc_ddc_mb(ddc) ((struct mrc_ddc_mb *) (ddc)->obj.subctx)

static int
mrc_patch_index(struct mrc_domain *domain, int b, int k[3])
{
  int (*blk_idx_p_idx)(struct mrc_domain *mb, int b, int k[3]);
  blk_idx_p_idx = (int (*)(struct mrc_domain *, int, int*)) 
    mrc_domain_get_method(domain, "block_idx_to_patch_idx");
  
  assert(blk_idx_p_idx);
  return blk_idx_p_idx(domain, b, k);
}

// Get the (single int) index of a three + patch index in
// the worldwide vector based on the number of ghost
// points given. Replaces __I3_patch macro
// gpatch is a global patch number (ie, old patch->nr)
static inline int
get_world_array_index(struct mrc_domain *mb, 
		 int gpatch, int s[3], 
		 int jx, int jy, int jz)
{
  
  struct mrc_domain_mb *sub = mrc_domain_mb(mb);
  struct mrc_patch *pmb = &sub->patches[gpatch];
  // p_off_shift (formerly a patch attribute) is the constant
  // part of the mapping from coordinate indicies to vector index.
  // Comes out if you expand out all the index
  // mappings.
  int p_off_shift = (((sub)->patches[0].ldims[0] + 2 * (s[0])) *
		     ((sub)->patches[0].ldims[1] + 2 * (s[1])) *
		     ((sub)->patches[0].ldims[2] + 2 * (s[2])) * 
		     (gpatch)) 
    - (((-s[2])*(pmb->ldims[1] + 2*s[1])+(-s[1]))*(pmb->ldims[0] + 2*s[0])+(-s[0]));

  assert(jx >= -s[0] && jx < (pmb)->ldims[0]+s[0]);
  assert(jy >= -s[1] && jy < (pmb)->ldims[1]+s[1]);
  assert(jz >= -s[2] && jz < (pmb)->ldims[2]+s[2]);
  return p_off_shift + 
    ((jz)*(pmb->ldims[1] + 2*s[1])+(jy))*(pmb->ldims[0] + 2*s[0])+(jx);
}

static void
_mrc_ddc_mb_destroy(struct mrc_ddc *ddc)
{
  struct mrc_ddc_mb *sub = to_mrc_ddc_mb(ddc);

  int ierr;
  for (int s = 0; s < MAX_SW; s++) {
    for (int bs = 0; bs < MAX_BS; bs++) {
      if (sub->mb_gtol[s][bs]) {
	ierr = VecScatterDestroy(&(sub->mb_gtol[s][bs])); CE;
      }
      if (sub->mb_ltol_edge[s][bs]) {
	ierr = VecScatterDestroy(&(sub->mb_ltol_edge[s][bs])); CE;
      }
    }
  }
  
}
// ----------------------------------------------------------------------
// mrc_ddc_multi_set_domain

static void
mrc_ddc_mb_set_domain(struct mrc_ddc *ddc, struct mrc_domain *domain)
{
  struct mrc_ddc_mb *sub = to_mrc_ddc_mb(ddc);

  sub->domain = domain;
}


struct mrc_domain *
mrc_ddc_mb_get_domain(struct mrc_ddc *ddc)
{
  struct mrc_ddc_mb *sub = to_mrc_ddc_mb(ddc);

  return sub->domain;
}



// ----------------------------------------------------------------------
// MB_find_face
//
// currently, this is only called for exterior points

static int
MB_find_face(struct mrc_domain *mb, int k, const int i[3])
{
  struct mrc_domain_mb *sub = mrc_domain_mb(mb);
  struct mrc_patch *patch = &sub->patches[k];
  int f = -1;
  for (int d = 0; d < 3; d++) {
    if (i[d] < 0) {
      f = 2*d;
      break;
    }
    if (i[d] >= patch->ldims[d]) {
      f = 2*d + 1;
      break;
    }
  }
  return f;
}

static int
map_to_interior(struct mrc_domain *mb, int bg, const int ig[3], int *pbl, int il[3])
{
  struct mrc_domain_mb *sub = mrc_domain_mb(mb);
  struct MB_block *gblock = &sub->mb_blocks[bg];
  // find face to neighbor
  int f = -1, in = 0;
  int is_exterior = 0;
  for (int d = 0; d < 3; d++) {
    if (ig[d] < 0) {
      is_exterior = 1;
      if (gblock->faces[2*d].block < 0)
	continue;
      f = 2*d;
      in = -ig[d]-1;
      break;
    }
    if (ig[d] >= gblock->mx[d]) {
      is_exterior = 1;
      if (gblock->faces[2*d+1].block < 0)
	continue;
      f = 2*d + 1;
      in = ig[d] - gblock->mx[d];
      break;
    }
  }
  // all interior or exterior but no connectivity info that direction
  if (f < 0) {
    if (is_exterior) {
      return -1;
    }
    il[0] = ig[0]; il[1] = ig[1]; il[2] = ig[2];
    *pbl = gblock->nr_block;
    return 0;
  }
  struct MB_face *face = &gblock->faces[f];
  int nface = face->face;
  assert(face->block >= 0);
  struct MB_block *nblock = &sub->mb_blocks[face->block];
  int ii[3]; // map to next block
  for (int d = 0; d < 3; d++) {
    if (face->map[d] == 0) {
      ii[face2dir(nface)] = 
	(face2bnd(nface) == 1) ? nblock->mx[face2dir(nface)]-1 - in: in;
    } else {
      ii[map2dir(face->map[d])] = (face->map[d] & MB_R) ? gblock->mx[d]-1 - ig[d] : ig[d];
    }
  }
  return map_to_interior(mb, face->block, ii, pbl, il);
}

static int
mb_patch_map_to_interior(struct mrc_domain *mb, int kg, const int ig[3], int *pkl, int il[3])
{
  int bg, _ig[3], bl, _il[3];

  struct mrc_domain_mb *sub = mrc_domain_mb(mb);
  struct mrc_patch_info ginfo;
  mrc_domain_get_global_patch_info(mb, kg, &ginfo);
  bg = ginfo.p_block;
  for (int d = 0; d < 3; d++) {
    _ig[d] = ig[d] + ginfo.p_ix[d];
  }
  if (map_to_interior(mb, bg, _ig, &bl, _il) < 0)
    return -1;

  int kl[3];
  for (int d = 0; d < 3; d++) {
    kl[d] = _il[d] / (sub->mb_blocks[bl].mx[d]/sub->ppb[d]);
  }
  *pkl = mrc_patch_index(mb, bl, kl);
  for (int d = 0; d < 3; d++) {
    il[d] = _il[d] - sub->patch_info[*pkl].p_ix[d];
  }
  return 0;
}


// This struct defines the mapping from logical coordinates 0,1,2
// (0 being normal to the face) to actual coordinates 

struct face_map {
  int idx[3]; // maps logical coord # to target coordinate #
  int bnd;    // which face boundary (lo/hi) are we on?
  int dim[3]; // target dimensions for 0 (normal) direction
};



static void
map_setup(struct face_map *map, struct mrc_domain *mb, int k, int f)
{
  struct mrc_domain_mb *sub = mrc_domain_mb(mb);
  map->idx[0] = face2dir(f);
  map->bnd = face2bnd(f);
  map->dim[0] = sub->patches[k].ldims[face2dir(f)];
}


static inline void
lmap_apply(const struct face_map *map, const int ix[3], int il[3])
{ 
  il[map->idx[0]] = (ix[0] < 0) ? ix[0] : map->dim[0] + ix[0];
  for (int d = 1; d < 3; d++) {
    il[map->idx[d]] = ix[d];
  }
}

static inline void
fill_ghost(struct mrc_domain *mb, int sw[3], int kl, const int il[3], int *idxl, int *idxg, int *pcnt)
{
  int ig[3], kg;
  int rc = mb_patch_map_to_interior(mb, kl, il, &kg, ig);
  if (rc) { // no point --> physical b.c.
    return;
  }
  idxl[*pcnt] = get_world_array_index(mb,kl,sw, il[0],il[1],il[2]);
  idxg[*pcnt] = get_world_array_index(mb,kg, (int[3]){0,0,0}, ig[0],ig[1],ig[2]);
  (*pcnt)++;
}

static void
map_and_fill_ghosts(struct mrc_domain *mb, int sw[3], const struct face_map *gmap, int k,
		    int *idxl, int *idxg, int *pcnt,
		    const int ixs[3], const int ixe[3])
{
  int ix[3];
  for (ix[0] = ixs[0]; ix[0] < ixe[0]; ix[0]++) {
    for (ix[1] = ixs[1]; ix[1] < ixe[1]; ix[1]++) {
      for (ix[2] = ixs[2]; ix[2] < ixe[2]; ix[2]++) {
	int il[3];
	lmap_apply(gmap, ix, il);
	fill_ghost(mb, sw, k, il, idxl, idxg, pcnt);
      }
    }
  }
}


static void
fill_face(struct mrc_domain *mb, int sw[3], struct face_map *lmap, int lk, int lf, int *idxl, int *idxg, 
	  int *pcnt)
{
  struct mrc_domain_mb *sub = mrc_domain_mb(mb);
  struct mrc_patch_info linfo;
  mrc_domain_get_global_patch_info(mb, lk, &linfo);
  struct MB_pface *lpface = &linfo.p_pface[lf];
  int gpatch_nr = lpface->pf_patch;
  struct mrc_patch *gpatch = &sub->patches[gpatch_nr];
  struct face_map smap;
  int rev_map[3];    // [1,2]: whether we have to reverse directions
  map_setup(&smap, mb, gpatch_nr, lpface->pf_face);
  int i = 1;
  for (int d = 0; d < 3; d++) {
    if (lpface->pf_map[d] != 0) {
      smap.idx[i] = map2dir(lpface->pf_map[d]);
      rev_map[i]  = lpface->pf_map[d] & MB_R;
      smap.dim[i] = gpatch->ldims[smap.idx[i]];
      assert(lmap->dim[i] == smap.dim[i]);
      i++;
    }
  }
  assert(i == 3);

  int ix[3];
  int ixs, ixe;
  // FIXME: Steve made this change, but he doesn't understand this function...
  // set the number of ghost points we need in the normal (dim[0]) direction
  int nr_ghosts = sw[lmap->idx[0]];
  if (lmap->bnd == 1) {
    ixs = 0; ixe = nr_ghosts;
  } else {
    ixs = -nr_ghosts; ixe = 0;
  }
  int mx[3] = { lmap->dim[0], lmap->dim[1], lmap->dim[2] };
  int cnt = *pcnt;
  for (ix[0] = ixs; ix[0] < ixe; ix[0]++) {
    for (ix[1] = 0; ix[1] < mx[1]; ix[1]++) {
      for (ix[2] = 0; ix[2] < mx[2]; ix[2]++) {
	int ig[3], is[3];
	lmap_apply(lmap, ix, ig);
	
	int isx0 = (lmap->bnd != smap.bnd) ? ix[0] : -ix[0]-1;
	is[smap.idx[0]] = (isx0 < 0) ? smap.dim[0] + isx0 : isx0;
	
	for (int d = 1; d < 3; d++) {
	  is[smap.idx[d]] = rev_map[d] ? smap.dim[d]-1 - ix[d]: ix[d];
	}
	
	idxl[cnt] = get_world_array_index(mb,lk,sw, ig[0],ig[1],ig[2]);
	assert(idxl[cnt] >= 0);
	idxg[cnt] = get_world_array_index(mb,gpatch_nr,(int[3]){0,0,0} , is[0],is[1],is[2]);
	assert(idxg[cnt] >= 0);
	cnt++;
      }
    }
  }
  *pcnt = cnt;
}

// ----------------------------------------------------------------------
// MB_GetVector

static int
MB_GetVector(struct mrc_domain *mb, int bs, int sw[3], Vec *pv)
{
  int ierr;
  Vec v;

  struct mrc_domain_mb *sub = mrc_domain_mb(mb);
  pfb;
  int nr_ghosts = MAX(sw[0], sw[1]);
  nr_ghosts = MAX(nr_ghosts, sw[2]);
  assert(nr_ghosts < MAX_SW);
  if (!mrc_domain_is_setup(mb)) {
    mrc_domain_setup(mb);
  }
  ierr = VecCreate(mrc_domain_comm(mb), &v); CE;
  //  ierr = VecSetSizes(v, bs*sub->mb_loc_N[sw], bs*sub->mb_N[sw]); CE;
  int loc_size = local_size(sub, sw);
  ierr = VecSetSizes(v, bs * loc_size, PETSC_DECIDE); CE;
  ierr = VecSetBlockSize(v, bs); CE;
  ierr = VecSetUp(v); CE;

 *pv = v;
  pfr;
}


// ----------------------------------------------------------------------
// MB_RestoreVector

static int
MB_RestoreVector(struct mrc_domain *mb, Vec *pv)
{
  int ierr;

  pfb;
  ierr = VecDestroy(pv); CE;
  pfr;
}



static int
getGtoL(struct mrc_ddc *ddc, int bs, int sw[3], VecScatter *pgtol)
{
  int ierr;

  struct mrc_domain *mb = mrc_ddc_mb_get_domain(ddc);
  struct mrc_domain_mb *sub = mrc_domain_mb(mb);
  struct mrc_ddc_mb *ddc_sub = to_mrc_ddc_mb(ddc);

  int nr_ghosts = MAX(sw[0],sw[1]);
  nr_ghosts = MAX(nr_ghosts, sw[2]);

  pfb;
  *pgtol = NULL;
  assert(bs-1 < MAX_BS);
  if (ddc_sub->mb_gtol[nr_ghosts][bs-1]) {
    *pgtol = ddc_sub->mb_gtol[nr_ghosts][bs-1];
    pfr;
  }

  int *idxg, *idxl, cnt=0;
  // FIXME, double the size shouldn't be necessary...
  int global_size = global_size(sub, sw);
  ierr = PetscMalloc(sizeof(*idxg) * global_size * 2, &idxg); CE;
  ierr = PetscMalloc(sizeof(*idxl) * global_size * 2, &idxl); CE;
  // OPT we don't need to recreate the index array for every bs
  MB_foreach_patch(mb, patch) {
    // interior
    mrc_patch_foreach(mb, patch, jx,jy,jz, 0, 0) {      
      idxl[cnt] = get_world_array_index(mb,sub->gpatch_off + patch,sw, jx,jy,jz);
      idxg[cnt] = get_world_array_index(mb,sub->gpatch_off + patch,(int[3]){0,0,0}, jx,jy,jz);
      cnt++;
    } mrc_patch_foreach_end;

    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mb, patch, &info);

    // ghost cells
    for (int f = 0; f < NR_FACES; f++) {
      struct MB_pface *pface = &info.p_pface[f];
      if (pface->pf_patch < 0) {
	continue;
      }
      struct face_map lmap;
      map_setup(&lmap, mb, info.global_patch, f);
      int f_sw[3]; // The ghost points mapped into the face normal frame.
      f_sw[0] = sw[lmap.idx[0]];
      int i = 1;
      for (int d = 0; d < 3; d++) {
	if (pface->pf_map[d] != 0) {
	  lmap.idx[i] = d;
	  lmap.dim[i] = info.ldims[lmap.idx[i]];
	  f_sw[i] = sw[lmap.idx[i]];
	  i++;
	}
      }
      assert(i == 3);
      fill_face(mb, sw, &lmap, info.global_patch, f, idxl, idxg, &cnt);

      int ixs, ixe;
      if (lmap.bnd == 1) {
	ixs = 0; ixe = f_sw[0];
      } else {
	ixs = -f_sw[0]; ixe = 0;
      }
      int mx[3] = { lmap.dim[0], lmap.dim[1], lmap.dim[2] };
      // edges
      map_and_fill_ghosts(mb, sw, &lmap, info.global_patch, idxl, idxg, &cnt,
			  (int[3]) { ixs, -f_sw[1]  ,  0        },
			  (int[3]) { ixe,  0      ,  mx[2]    });
      map_and_fill_ghosts(mb, sw, &lmap, info.global_patch, idxl, idxg, &cnt,
			  (int[3]) { ixs, mx[1]   ,  0        },
			  (int[3]) { ixe, mx[1]+f_sw[1], mx[2]    });
      map_and_fill_ghosts(mb, sw, &lmap, info.global_patch, idxl, idxg, &cnt,
			  (int[3]) { ixs, 0      ,  -f_sw[2]  },
			  (int[3]) { ixe, mx[1]  ,  0         });
      map_and_fill_ghosts(mb, sw, &lmap, info.global_patch, idxl, idxg, &cnt,
			  (int[3]) { ixs, 0      ,  mx[2]     },
			  (int[3]) { ixe, mx[1]  ,  mx[2]+f_sw[2] });
      // corners
      map_and_fill_ghosts(mb, sw, &lmap, info.global_patch, idxl, idxg, &cnt,
			  (int[3]) { ixs, -f_sw[1] ,  -f_sw[2]   },
			  (int[3]) { ixe, 0       ,  0        });
      map_and_fill_ghosts(mb, sw, &lmap, info.global_patch, idxl, idxg, &cnt,
			  (int[3]) { ixs, -f_sw[1] ,  mx[2]    },
			  (int[3]) { ixe, 0       ,  mx[2]+f_sw[2] });
      map_and_fill_ghosts(mb, sw, &lmap, info.global_patch, idxl, idxg, &cnt,
			  (int[3]) { ixs, mx[1]   ,  -f_sw[2]    },
			  (int[3]) { ixe, mx[1]+f_sw[1],  0        });
      map_and_fill_ghosts(mb, sw, &lmap, info.global_patch, idxl, idxg, &cnt,
			  (int[3]) { ixs, mx[1]   ,  mx[2]    },
			  (int[3]) { ixe, mx[1]+f_sw[1],  mx[2]+f_sw[2] });
    }
    //    ASSERT(cnt <= mb->gN); FIXME!!!
  } MB_foreach_patch_end;

  /*********************************
   * With petsc > 3.3 the indexing is relative to block, not element
   * so this extra step is not needed
   * for (int i = 0; i < cnt; i++) {
   * idxl[i] *= bs;
   * idxg[i] *= bs;
   * }
  */

  IS isl, isg;
  ierr = ISCreateBlock(mrc_domain_comm(mb), bs, cnt, idxl, PETSC_COPY_VALUES, &isl); CE;
  ierr = ISCreateBlock(mrc_domain_comm(mb), bs, cnt, idxg, PETSC_COPY_VALUES, &isg); CE;
  ierr = PetscFree(idxl); CE;
  ierr = PetscFree(idxg); CE;

  Vec lvec, gvec;
  ierr = MB_GetVector(mb, bs, (int[3]){0,0,0}, &gvec); CE;
  ierr = MB_GetVector(mb, bs, sw, &lvec); CE;
  ierr = VecScatterCreate(gvec, isg, lvec, isl, &ddc_sub->mb_gtol[nr_ghosts][bs-1]); CE;
  ierr = MB_RestoreVector(mb, &gvec); CE;
  ierr = MB_RestoreVector(mb, &lvec); CE;
  ierr = ISDestroy(&isl); CE;
  ierr = ISDestroy(&isg); CE;

  *pgtol = ddc_sub->mb_gtol[nr_ghosts][bs-1];
  pfr;
}

// ======================================================================
// ddc_GlobalToLocal
// FIXME: We should figure out a way to impliment this for the rest of the
// ddc classes too. Shouldn't actually be that hard...

static void
mrc_ddc_mb_global_to_local_fld(struct mrc_ddc *ddc, struct mrc_fld *gfld, struct mrc_fld *lfld)
{
  int ierr;
  assert(strcmp(mrc_ddc_type(ddc), "mb") == 0);
  // (I don't know that this comment pertained to, but it was probably important)
  // but it doesn't exist yet, and I don't feel like writing
  // it in

  // FIXME: I shouldn't have to access private attributes here.
  assert(gfld->_nr_ghosts == 0);
  VecScatter gtol;
  // FIXME: I can't make this an op because I can't include a petsc file at the top level
  // in mrc_fld. 

  fgp_t fld_get_petsc = (fgp_t) mrc_fld_get_method(gfld, "get_petsc_vec");
  Vec gfld_v = fld_get_petsc(gfld);
  Vec lfld_v = fld_get_petsc(lfld);


  // we're going to access private mrc_fld stuff, because I'm to lazy to write in 
  // a get_param_int_array
  assert(lfld->_sw.nr_vals == 5);
  int sw[3];
  for(int d = 0; d < 3; d++) {
    sw[d] = lfld->_sw.vals[d];
  }

  ierr = getGtoL(ddc, mrc_fld_nr_comps(lfld), sw, &gtol); CE;

  ierr = VecScatterBegin(gtol, gfld_v, lfld_v, 
			 INSERT_VALUES, SCATTER_FORWARD); CE;
  ierr = VecScatterEnd  (gtol, gfld_v, lfld_v, 
			 INSERT_VALUES, SCATTER_FORWARD); CE;

  fpp_t fld_put_petsc = (fpp_t) mrc_fld_get_method(gfld, "put_petsc_vec");
  fld_put_petsc(gfld, &gfld_v);
  fld_put_petsc(lfld, &lfld_v);

}



// p is GLOBAL patch number
static void
fill_bnd_edge(struct mrc_domain *mb, int jxs[3], int jxe[3], int p, int sw[3],
	      int *idxg, int *idxl, int *pcnt)
{
  struct mrc_domain_mb *sub = mrc_domain_mb(mb);
  struct mrc_patch *patch = &sub->patches[p];
  struct mrc_patch_info info;
  mrc_domain_get_global_patch_info(mb, p, &info);
  int *mx = patch->ldims;
  int jx[3];
  for (jx[2] = jxs[2]; jx[2] < jxe[2]; jx[2]++) {
    for (jx[1] = jxs[1]; jx[1] < jxe[1]; jx[1]++) {
      for (jx[0] = jxs[0]; jx[0] < jxe[0]; jx[0]++) {
	for (int dir = 0; dir < 3; dir++) {
	  int f, jn;
	  if (jx[dir] < 0) {
	    f = 0;
	    jn = -jx[dir] - 1;
	  } else if (jx[dir] >= mx[dir]) {
	    f = 1;
	    jn = jx[dir] - mx[dir];
	  } else {
	    continue;
	  }
	  int npatch = info.p_pface[2*dir+f].pf_patch;
	  if (npatch < 0) // external?
	    continue;

	  // I would prefer to use get_global_patch_info here, but I'm not sure it's worth it.
	  int nface = sub->patch_info[npatch].p_pface[2*dir+f].pf_face;
	  //	  printf("patch %d npatch %d\n", patch->p_nr, npatch);
	  //	  printf("face %d, nface %d\n", 2*dir+f, nface);
	  assert(nface == 2*dir + (1-f)); 
	  // otherwise, we have to remap indices, directions

	  int ix[3];
	  int *nmx = sub->patches[npatch].ldims;
	  for (int d = 0; d < 3; d++) {
	    if (d == dir) {
	      ix[d] = face2bnd(nface) == 1 ? nmx[d] - jn - 1 : jn;
	    } else {
	      ix[d] = jx[d];
	    }
	  }
	  idxl[*pcnt] = get_world_array_index(mb,p,sw, jx[0],jx[1],jx[2]);
	  idxg[*pcnt] = get_world_array_index(mb,npatch,sw, ix[0],ix[1],ix[2]);
	  (*pcnt)++;
	}
      }
    }
  }
}
static int
getLtoL_Edge(struct mrc_ddc *ddc, int bs, int sw[3], VecScatter *pltol)
{
  int ierr;
  struct mrc_domain *mb = mrc_ddc_mb_get_domain(ddc);
  struct mrc_domain_mb *sub = mrc_domain_mb(mb);
  struct mrc_ddc_mb *ddc_sub = to_mrc_ddc_mb(ddc);

  pfb;

  *pltol = NULL;
  assert(bs-1 < MAX_BS);

  int nr_ghosts = MAX(sw[0], sw[1]);
  nr_ghosts = MAX(nr_ghosts, sw[2]);

  if (ddc_sub->mb_ltol_edge[nr_ghosts][bs-1]) {
    *pltol = ddc_sub->mb_ltol_edge[nr_ghosts][bs-1];
    pfr;
  }

  int *idxg, *idxl, cnt=0;

 // FIXME: this should be more than it needs.. not a big deal
  int local_size = local_size(sub,((int[3]){1, 1, 1}) );
  ierr = PetscMalloc(sizeof(*idxg) * local_size, &idxg); CE;
  ierr = PetscMalloc(sizeof(*idxl) * local_size, &idxl); CE;
  // OPT we don't need to recreate the index array for every bs

  MB_foreach_patch(mb, patch) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(mb, patch, &info);
    int *mx = info.ldims;
    int jxs[3], jxe[3];
    // FIXME: I'm going to take a gamble here and assume that because we're using
    // patch ldims that dir represents regular x,y,z numbering.
    for (int dir = 0; dir < 1; dir++) {
      for (int f = 0; f < 2; f++) {
	if (info.p_pface[2*dir+f].pf_patch < 0) { // external?
	  if (f == 0) {
	    jxs[dir] = -sw[dir]; jxe[dir] = 0;
	  } else {
	    jxs[dir] = mx[dir]; jxe[dir] = mx[dir]+sw[dir];
	  }
	  
	  for (int dir1 = dir+1; dir1 < 3; dir1++) {
	    for (int f1 = 0; f1 < 2; f1++) {
	      if (info.p_pface[2*dir1+f1].pf_patch < 0)  // external?
		continue;
	      
	      if (f1 == 0) {
		jxs[dir1] = -sw[dir1]; jxe[dir1] = 0;
	      } else {
		jxs[dir1] = mx[dir1]; jxe[dir1] = mx[dir1]+sw[dir1];
	      }
	      for (int dir2 = 0; dir2 < 3; dir2++) {
		if (dir2 == dir || dir2 == dir1)
		  continue;
		
		jxs[dir2] = 0; jxe[dir2] = mx[dir2];
		fill_bnd_edge(mb, jxs, jxe, sub->gpatch_off + patch, sw, idxg, idxl, &cnt);
	      }
	    }
	  }
	}
      }
      
    }
  } MB_foreach_patch_end;

  /*********************************
   * With petsc > 3.3 the indexing is relative to block, not element
   * so this extra step is not needed
   *
   *  for (int i = 0; i < cnt; i++) {
   *    idxl[i] *= bs;
   *    idxg[i] *= bs;
   *  }
   */

  IS isl, isg;
  ierr = ISCreateBlock(mrc_domain_comm(mb), bs, cnt, idxl, PETSC_COPY_VALUES, &isl); CE;
  ierr = ISCreateBlock(mrc_domain_comm(mb), bs, cnt, idxg, PETSC_COPY_VALUES, &isg); CE;
  ierr = PetscFree(idxl); CE;
  ierr = PetscFree(idxg); CE;

  Vec lvec;
  ierr = MB_GetVector(mb, bs, sw, &lvec); CE;
  ierr = VecScatterCreate(lvec, isg, lvec, isl, 
			  &ddc_sub->mb_ltol_edge[nr_ghosts][bs-1]); CE;
  ierr = MB_RestoreVector(mb, &lvec); CE;
  ierr = ISDestroy(&isl); CE;
  ierr = ISDestroy(&isg); CE;

  *pltol = ddc_sub->mb_ltol_edge[nr_ghosts][bs-1];
  pfr;
}

// FIXME: ugly name
static void
mrc_ddc_mb_fill_ghost_edges_fld(struct mrc_ddc *ddc, int mb, int me, 
				struct mrc_fld *x)
#warning FIXME mb_fill_ghost_edges will not do single components, just full block
{
  int ierr;

  VecScatter ltol;
  // FIXME: shouldn't have to access private var to get nr_ghosts
  assert(x->_sw.nr_vals == 5);
  int sw[3];
  for(int d = 0; d < 3; d++) {
    sw[d] = x->_sw.vals[d];
  }  
  ierr = getLtoL_Edge(ddc, mrc_fld_nr_comps(x), sw, &ltol); CE;

  fgp_t fld_get_petsc = (fgp_t) mrc_fld_get_method(x, "get_petsc_vec");
  Vec x_v = fld_get_petsc(x);
  ierr = VecScatterBegin(ltol, x_v, x_v, 
			 INSERT_VALUES, SCATTER_FORWARD); CE;
  ierr = VecScatterEnd  (ltol, x_v, x_v, 
			 INSERT_VALUES, SCATTER_FORWARD); CE;
  fpp_t fld_put_petsc = (fpp_t) mrc_fld_get_method(x, "put_petsc_vec");
  fld_put_petsc(x, &x_v);
}


static void
mrc_ddc_mb_fill_ghosts_fld(struct mrc_ddc *ddc, int mb, int me,
			      struct mrc_fld *fld)
{
#warning FIXME mb fill_ghosts will not do single componenets, just full block
  // FIXME: Slow, has to create a field.

  int bs = mrc_fld_nr_comps(fld);
  struct mrc_domain *domain = fld->_domain;

  struct mrc_fld *xg = mrc_domain_fld_create(domain, 0, NULL);
  mrc_fld_set_type(xg, mrc_fld_type(fld));
  mrc_fld_set_param_int(xg, "nr_comps", bs);
  mrc_fld_setup(xg);  

  MB_foreach_patch(domain, patch) {
    mrc_patch_foreach(domain, patch, jx,jy,jz, 0,0) {
      for (int m = 0; m < bs; m++) {
	MRC_D5(xg, m, jx,jy,jz, patch) = MRC_D5(fld, m, jx,jy,jz, patch);
      }
    } mrc_patch_foreach_end;
  } MB_foreach_patch_end;
  
  mrc_ddc_global_to_local_fld(domain->ddc, xg, fld);

  // FIXME: not sure if I really need to call this. Think so..

  mrc_ddc_mb_fill_ghost_edges_fld(ddc, 0, bs, fld);
  mrc_fld_destroy(xg);


}

// ======================================================================
// Methods
// ====================
// These are some functions particular to the MB domain that map patches
// and indices across boundaries taking into account the block mappings.
// They could use a better name and API, and intrinsically rely on petsc
// at the moment.


// Maps from the local patch (lpatch) interior coordinates to the global
// domain coordinates. Note this is not field specific.
// Used to be called mrc_matrix_find_global_index
static int
map_patch_interior_index_to_global(struct mrc_domain *domain, 
         int lpatch, 
         int jx, int jy, int jz)
{
  struct mrc_domain_mb *mb = mrc_domain_mb(domain);
  return get_world_array_index(domain, mb->gpatch_off + lpatch,
             (int[3]){0,0,0}, jx, jy, jz);
}

// Maps from an arbitrary coordinate relative to global patch k1 to some interior
// point in the global array. Unlike the above function this one will map points
// exterior to a given patch into their appropriate place in the domain.
// Used to be called MB_get_index
static int
map_global_patch_index_to_global(struct mrc_domain *mb, int k1, int i1x, int i1y, int i1z)
{
  int i1[3] = { i1x, i1y, i1z }, k2, i2[3];
  int ierr = mb_patch_map_to_interior(mb, k1, i1, &k2, i2);
  if (ierr) {
    return -1;
  }
  return get_world_array_index(mb,k2,(int[3]){0,0,0}, i2[0], i2[1], i2[2]);
}

// ======================================================================
// mrc_ddc_multi_ops

#define VAR(x) (void *)offsetof(struct mrc_ddc_mb, x)
static struct param mrc_ddc_mb_descr[] = {
  { "domain"          , VAR(domain)       , PARAM_OBJ(mrc_domain) },
  {},
};
#undef VAR

static struct mrc_obj_method mrc_ddc_mb_methods[] = {
  MRC_OBJ_METHOD("map_patch_interior_index_to_global", map_patch_interior_index_to_global),
  MRC_OBJ_METHOD("map_global_patch_index_to_global", map_global_patch_index_to_global),
  {}
};

struct mrc_ddc_ops mrc_ddc_mb_ops = {
  .name                  = "mb",
  .size                  = sizeof(struct mrc_ddc_mb),
  .param_descr           = mrc_ddc_mb_descr,
  .methods               = mrc_ddc_mb_methods,
  .destroy               = _mrc_ddc_mb_destroy,
  .set_domain            = mrc_ddc_mb_set_domain,
  .get_domain            = mrc_ddc_mb_get_domain,
  .fill_ghosts_fld       = mrc_ddc_mb_fill_ghosts_fld,
  .global_to_local_fld   = mrc_ddc_mb_global_to_local_fld,
  .fill_ghost_edges_fld  = mrc_ddc_mb_fill_ghost_edges_fld,
};
