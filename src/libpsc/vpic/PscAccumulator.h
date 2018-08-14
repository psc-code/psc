
#pragma once

template<typename MfieldsAccumulator, typename MfieldsState>
struct PscAccumulatorOps
{
  using Grid = typename MfieldsAccumulator::Grid;
  
  // ----------------------------------------------------------------------
  // clear

  static void clear(MfieldsAccumulator& acc)
  {
    auto g = acc.grid();
    int n_array = acc.n_pipeline() + 1;
    for (int arr = 0; arr < n_array; arr++) {
      auto a_begin = &acc(arr, 1,1,1);
      auto a_end = &acc(arr, g->nx, g->ny, g->nz);
      // FIXME, the + 1 in n0 doesn't really make sense to me.  And,
      // originally, this was extended to 128 byte boundaries, too,
      // which I dropped -- which is also a behavior change, which I
      // though shouldn't matter as it's outside the local domain, but
      // it might, too
      memset(a_begin, 0, (a_end - a_begin + 1) * sizeof(*a_begin));
    }
  }

  // ----------------------------------------------------------------------
  // reduce
  
  static void reduce(MfieldsAccumulator& acc)
  {
    auto g = acc.grid();
    int si = sizeof(typename MfieldsAccumulator::Element) / sizeof(float);
    int nr = acc.n_pipeline() + 1 - 1;
    int sr = si * acc.stride();

    auto a_begin = &acc(0,1,1,1);
    auto a_end = &acc(0, g->nx, g->ny, g->nz);
    int n = a_end - a_begin + 1;

    // a is broken into restricted rw and ro parts to allow the compiler
    // to do more aggresive optimizations

    float * RESTRICT a = reinterpret_cast<float *>(a_begin);
    const float * RESTRICT ALIGNED(16) b = a + sr;

    float f[si];

    for(int i = 0; i < n; i++) {
      int j = i*si;
      for (int m = 0; m < si; m++) {
	f[m] = a[j+m];
      }
      for (int r = 0; r < nr; r++) {
	int k = j + r*sr;
	for (int m = 0; m < si; m++) {
	  f[m]  += b[k+m];
	}
      }
      for (int m = 0; m < si; m++) {
	a[j+m] = f[m];
      }
    }
  }

  // ----------------------------------------------------------------------
  // unload
  
  static void unload( /*const*/ MfieldsAccumulator& acc, MfieldsState& mflds)
  {
    auto g = acc.grid();
    float cx = 0.25 * g->rdy * g->rdz / g->dt;
    float cy = 0.25 * g->rdz * g->rdx / g->dt;
    float cz = 0.25 * g->rdx * g->rdy / g->dt;
    
    auto& fa = mflds.getPatch(0);
    Field3D<typename MfieldsState::Patch> F(fa);
    Field3D<MfieldsAccumulator> A(acc);

    int nx = g->nx, ny = g->ny, nz = g->nz;
    // FIXME, these limits seem to go too far out compared to what we zeroed before
    for(int k = 1; k <= nz+1; k++) {
      for(int j = 1; j <= ny+1; j++) {
	for(int i = 1; i <= nx+1; i++) {
	  F(i,j,k).jfx += cx * (A(i,j,k).jx[0] + A(i,j-1,k).jx[1] + A(i,j,k-1).jx[2] + A(i,j-1,k-1).jx[3]);
	  F(i,j,k).jfy += cy * (A(i,j,k).jy[0] + A(i,j,k-1).jy[1] + A(i-1,j,k).jy[2] + A(i-1,j,k-1).jy[3]);
	  F(i,j,k).jfz += cz * (A(i,j,k).jz[0] + A(i-1,j,k).jz[1] + A(i,j-1,k).jz[2] + A(i-1,j-1,k).jz[3]);
	}
      }
    }
  }
};

