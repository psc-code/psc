#include "psc.h"
#include "psc_ppu.h"

#include <string.h>

// Built off Kai's counting sort 2 method, 
// which uses a counting sort with a seperate array of cell indices. 
// This will reduce what we need to transfer the the spes if I can ever
// cut the number of variables per particle down to 8.

#if PARTICLES_BASE == PARTICLES_FORTRAN
static inline in
__nint(real x)
{
  return (int)(x + real(10.5)) - 10;
}

static inline void
find_cell(real xi, real yi, real zi, int l[3])
{
  l[0] = __nint(xi / psc.dx[0]) - psc.ilo[0]; 
  l[1] = __nint(yi / psc.dx[1]) - psc.ilo[1];
  l[2] = __nint(zi / psc.dx[2]) - psc.ilo[2];
}


static void
cbe_find_cell_indices(particles_base_t *pp)
{
  int block_size[3] = spu_ctl.block_size;
  int blkgd[3] = spu_ctl.block_grid;

  for (int i = 0; i < pp->n_part; i++) {
    particle_base_t *p = particles_base_get_one(pp,i);

    int ci[3];
    find_cell(p->xi,p->yi,p->zi,ci);
    ci[0] /= block_size[0];
    ci[1] /= block_size[1];
    ci[2] /= block_size[2];
    p->cni = ((ci[2] 
	      * blkgd[1] + ci[1])
	      * blkgd[0] + ci[0]);
    assert(p->cni < spu_ctl.blocks);
  }
}

static void
cbe_countsort()
{
  static int pr; 
  if (!pr) {
    pr = prof_register("cbe_countsort",1.,0,0);
  }
  
  find_cell_indices(&psc.pp);
  
  int N = 2 * spu_ctl.nblocks; 

  // For offloading to the spes we need to save the offsets 
  // of each block in the particle array. 
  unsigned in *cnis = malloc(psc.pp.n_part * sizeof(*cnis));
  for (int i = 0; i < psc.pp.n_part; i++) {
    cnis[i] = get_cell_index(&psc.pp.particles[i]);
    assert(cnis[i] < N);
  }
  
  prof_start(pr);

  unsigned int *cnts;

  if (! spu_ctl.cnts) {
    cnts = malloc(N * sizeof(*cnts));
    spu_ctl.cnts = cnts;
  }
  else {
    cnts = spu_ctl.cnts;
  }

  memset(cnts, 0, N * sizeof(*cnts));

  // count
  for (int i = 0; i < psc.pp.n_part; i++) {
    unsigned int cni = cnis[i];
    cnts[cni]++;
  }

  // calc offsets
  int cur = 0;
  for(int i = 1; i < N; i++) {
    int n = cnts[i];
    cnts[i] = cur;
    cur += n;
  }
  assert(cur == psc.pp.n_part);

  // move into new position
  particle_base_t *particles2 = malloc(psc.pp.n_part * sizeof(*particles2));
  for(int i = 0; i < psc.pp.n_part; i++) {
    unsigned int cni = cnis[i];
    int n = 1;
    while (i+n < psc.pp.n_part && cnis[i+n] == cni) {
      n++;
    }
    memcpy(&particles2[cnts[cni]], &psc.pp.particles[i], n * sizeof(*particles2));
    cnts[cni] += n;
    i += n - 1;
  }

  // back to in-place
  memcpy(psc.pp.particles, particles2, psc.pp.n_part * sizeof(*particles2));

  free(particles2);
  free(cnts);

  // Need to assign particle beginning/ending
  // addresses to the blocks, which means we 
  // get to have a little fun. 

  psc_cell_block_t * curr = *blocks_list;
  for(int i = 1; i <= spu_ctl.nblocks; i++){
    curr->part_start = (unsigned long long) (psc.pp.particles + cnts[i]);
    curr->part_end = (unsigned long long) (psc.pp.particles + cnts[i+1]);
    cur++;
  }

  spu_ctl.blocks_ready = 1;
  
  prof_stop(pr);
}

struct psc_sort_ops psc_sort_ops_cbe = {
  .name = "cbesort",
  .sort = cbe_countsort,
};

/// \file cbe_sort.h 
///
/// The Cell, like Cuda, needs a bit of a customized
/// sorting call. The algorithm isn't particularly
/// unique, but a we have a little bit of information
/// which needs to be stored, and the grain of the 
/// sort is much coarser than what is needed for
/// the other methods. 
