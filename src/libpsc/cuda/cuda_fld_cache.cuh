
#pragma once

// ======================================================================
// FldCache

template<typename BS>
struct FldCache
{
  const static int BLOCKSIZE_X = BS::x::value;
  const static int BLOCKSIZE_Y = BS::y::value;
  const static int BLOCKSIZE_Z = BS::z::value;

  using real_t = float;
  
  __device__ FldCache() = default;
  __device__ FldCache(const FldCache&) = delete;
  
  __device__ void load(DFields d_flds, int *ci0)
  {
    off_ = (-(ci0[2] - 2) * (BLOCKSIZE_Y + 4) +
	    -(ci0[1] - 2));
#if 0
    ci0y_ = ci0[1];
    ci0z_ = ci0[2];
#endif
    
    int ti = threadIdx.x;
    int n = BLOCKSIZE_X * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4);
    while (ti < n) {
      int tmp = ti;
      int jy = tmp % (BLOCKSIZE_Y + 4) - 2;
      tmp /= BLOCKSIZE_Y + 4;
      int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
      // OPT? currently it seems faster to do the loop rather than do m by threadidx
      for (int m = EX; m <= HZ; m++) {
#if 0
	float val = d_flds(m, 0,jy+ci0[1],jz+ci0[2]);
	//printf("C load %d %d: %g (%d)\n", jy+ci0[1], jz+ci0[2], val, m);
	if (!isfinite(val)) {
	  printf("CUDA_ERROR: load %g m %d jz %d %d ci0 %d %d\n", val, m, jy, jz, ci0[1], ci0[2]);
	}
#endif
	(*this)(m, 0, jy+ci0[1], jz+ci0[2]) = d_flds(m, 0,jy+ci0[1],jz+ci0[2]);
      }
      ti += THREADS_PER_BLOCK;
    }
  }

  __host__ __device__ float operator()(int m, int i, int j, int k) const
  {
    return data_[index(m, i,j,k)];
  }

private: // it's supposed to be a (read-only) cache, after all
  __host__ __device__ float& operator()(int m, int i, int j, int k)
  {
    return data_[index(m, i,j,k)];
  }

private:
  __host__ __device__ int index(int m, int i, int j, int k) const
  {
#if 0
    if (j < ci0y_ -2 || j >= ci0y_ + BLOCKSIZE_Y + 2 ||
	k < ci0z_ -2 || k >= ci0z_ + BLOCKSIZE_Y + 2) {
      printf("CUDA_ERROR: fld cache jk %d %d ci0 %d %d\n",
	     j, k, ci0y_, ci0z_);
    }
#endif
    return (((m - EX) * (BLOCKSIZE_Z + 4)
	     + k) * (BLOCKSIZE_Y + 4)
	    + j) + off_;
  }

  float data_[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
  int off_;
#if 0
  int ci0y_, ci0z_;
#endif
};

