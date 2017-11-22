
#ifndef PSC_ACCUMULATOR_OPS_H
#define PSC_ACCUMULATOR_OPS_H

template<class A>
struct PscAccumulatorOps {
  typedef A Accumulator;

#if 0
#define VOX(x,y,z) VOXEL(x,y,z, g->nx,g->ny,g->nz)

  enum { accumulators_n_block = 256 };

  // this version still kinda has the original pipeline handling,
  
  void
  clear_accumulator_( accumulator_array_t * RESTRICT aa )
  {
    int i0 = VOX(1,1,1);
    int n0      = VOX(aa->g->nx,aa->g->ny,aa->g->nz) - i0 + 1;
    int n_array = aa->n_pipeline + 1;
    int s_array = aa->stride;
    accumulator_t* a0 = aa->a + i0;
    for (int rank = 0; rank < n_array; rank++) {
      int n = n0;
      int i;
      DISTRIBUTE(n, accumulators_n_block, rank, aa->n_pipeline, i, n);
      accumulator_t *a = a0 + i;
      for(int cnt = n_array; cnt; cnt--, a += s_array ) {
	mprintf("rank %d clear a %ld to %ld\n", rank, a-aa->a, a-aa->a+n);
	CLEAR( a, n );
      }
    }
  }
#undef VOX

#endif
  
  void
  clear_accumulator(Accumulator& accumulator)
  {
    grid_t *g = accumulator.g;
    
    int n_array = accumulator.n_pipeline + 1;
    for(int arr = 0; arr < n_array; arr++) {
      accumulator_t *a_begin = &accumulator(arr, 1,1,1);
      accumulator_t *a_end = &accumulator(arr, g->nx, g->ny, g->nz);
      // FIXME, the + 1 in n0 doesn't really make sense to me.  And,
      // originally, this was extended to 128 byte boundaries, too,
      // which I dropped -- which is also a behavior change, which I
      // though shouldn't matter as it's outside the local domain, but
      // it might, too
      CLEAR(a_begin, a_end - a_begin + 1);
    }
  }
  void clear_accumulator_array(Accumulator *accumulator)
  {
    TIC clear_accumulator(*accumulator); TOC(clear_accumulators, 1);
  }
};


#endif
