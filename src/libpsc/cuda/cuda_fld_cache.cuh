
#pragma once

// ======================================================================
// FldCache

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
struct FldCache
{
  using real_t = float;
  
  __device__ FldCache() = default;
  __device__ FldCache(const FldCache&) = delete;
  
  __device__ void load(DFields d_flds, int *ci0)
  {
    off_ = (-(ci0[2] - 2) * (BLOCKSIZE_Y + 4) +
	    -(ci0[1] - 2));
    
    int ti = threadIdx.x;
    int n = BLOCKSIZE_X * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4);
    while (ti < n) {
      int tmp = ti;
      int jy = tmp % (BLOCKSIZE_Y + 4) - 2;
      tmp /= BLOCKSIZE_Y + 4;
      int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
      // OPT? currently it seems faster to do the loop rather than do m by threadidx
      for (int m = EX; m <= HZ; m++) {
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
    return (((m - EX) * (BLOCKSIZE_Z + 4)
	     + k) * (BLOCKSIZE_Y + 4)
	    + j) + off_;
  }

  float data_[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
  int off_;
};

