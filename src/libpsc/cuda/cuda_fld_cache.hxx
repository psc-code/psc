
#pragma once

// #define DEBUG_FLDCACHE

// ======================================================================
// FldCache

template <typename BS, typename DIM>
struct FldCache
{
  using dim = DIM;
  using value_type = float;

  const static int BLOCKSIZE_X = BS::x::value;
  const static int BLOCKSIZE_Y = BS::y::value;
  const static int BLOCKSIZE_Z = BS::z::value;

  FldCache() = default;
  __device__ FldCache(const FldCache&) = delete;

  template <typename E>
  __device__ void load(E&& gt, Int3 ib, int* ci0)
  {
    off_ =
      ((-(ci0[2] - 2) * (BLOCKSIZE_Y + 4) + -(ci0[1] - 2)) * (BLOCKSIZE_X + 4) +
       -(ci0[0] - 2));

#ifdef DEBUG_FLDCACHE
    ci0_[0] = ci0[0];
    ci0_[1] = ci0[1];
    ci0_[2] = ci0[2];
#endif

    auto _d_flds = make_Fields3d<dim_xyz>(gt, ib);
    int n = (BLOCKSIZE_X + 4) * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4);
    for (int ti = threadIdx.x; ti < n; ti += THREADS_PER_BLOCK) {
      int tmp = ti;
      int jx = tmp % (BLOCKSIZE_X + 4) - 2;
      tmp /= BLOCKSIZE_X + 4;
      int jy = tmp % (BLOCKSIZE_Y + 4) - 2;
      tmp /= BLOCKSIZE_Y + 4;
      int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
      for (int m = EX; m <= HZ; m++) {
#ifdef DEBUG_FLDCACHE
        float val = _d_flds(m, jx + ci0[0], jy + ci0[1], jz + ci0[2]);
        // printf("C load %d %d %d: %g (%d)\n", jx+ci0[0], jy+ci0[1], jz+ci0[2],
        // val, m);
        if (!isfinite(val)) {
          printf("CUDA_ERROR: load %g m %d jxyz %d %d %d ci0 %d %d %d\n", val,
                 m, jx, jy, jz, ci0[0], ci0[1], ci0[2]);
        }
#endif
        (*this)(m, jx + ci0[0], jy + ci0[1], jz + ci0[2]) =
          _d_flds(m, jx + ci0[0], jy + ci0[1], jz + ci0[2]);
      }
    }
  }

  __host__ __device__ float operator()(int m, int i, int j, int k) const
  {
    return data_[index(m, i, j, k)];
  }

private: // it's supposed to be a (read-only) cache, after all
  __host__ __device__ float& operator()(int m, int i, int j, int k)
  {
    return data_[index(m, i, j, k)];
  }

private:
  __host__ __device__ int index(int m, int i, int j, int k) const
  {
#ifdef DEBUG_FLDCACHE
    if (i < ci0_[0] - 2 || i >= ci0_[0] + BLOCKSIZE_X + 2 || j < ci0_[1] - 2 ||
        j >= ci0_[1] + BLOCKSIZE_Y + 2 || k < ci0_[2] - 2 ||
        k >= ci0_[2] + BLOCKSIZE_Z + 2) {
      printf("CUDA_ERROR: fld cache ijk %d %d %d ci0 %d %d %d\n", i, j, k,
             ci0_[0], ci0_[1], ci0_[2]);
    }
#endif
    return ((((m - EX) * (BLOCKSIZE_Z + 4) + k) * (BLOCKSIZE_Y + 4) + j) *
              (BLOCKSIZE_X + 4) +
            i) +
           off_;
  }

  float data_[6 * (BLOCKSIZE_X + 4) * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
  int off_;
#ifdef DEBUG_FLDCACHE
  int ci0_[3];
#endif
};

// ======================================================================
// FldCache dim_yz specialization

template <typename BS>
struct FldCache<BS, dim_yz>
{
  using dim = dim_yz;
  using value_type = float;

  static_assert(BS::x::value == 1, "FldCache: dim_yz needs BS::x == 1");
  const static int BLOCKSIZE_X = BS::x::value;
  const static int BLOCKSIZE_Y = BS::y::value;
  const static int BLOCKSIZE_Z = BS::z::value;

  FldCache() = default;
  __device__ FldCache(const FldCache&) = delete;

  template <typename E>
  __device__ void load(E&& gt, Int3 ib, int* ci0)
  {
    off_ = (-(ci0[2] - 2) * (BLOCKSIZE_Y + 4) + -(ci0[1] - 2));
#ifdef DEBUG_FLDCACHE
    ci0y_ = ci0[1];
    ci0z_ = ci0[2];
#endif

    auto _d_flds = make_Fields3d<dim_yz>(gt, ib);
    int ti = threadIdx.x;
    int n = BLOCKSIZE_X * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4);
    while (ti < n) {
      int tmp = ti;
      int jy = tmp % (BLOCKSIZE_Y + 4) - 2;
      tmp /= BLOCKSIZE_Y + 4;
      int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
      // OPT? currently it seems faster to do the loop rather than do m by
      // threadidx
      for (int m = EX; m <= HZ; m++) {
#ifdef DEBUG_FLDCACHE
        float val = _d_flds(m, 0, jy + ci0[1], jz + ci0[2]);
        // printf("C load %d %d: %g (%d)\n", jy+ci0[1], jz+ci0[2], val, m);
        if (!isfinite(val)) {
          printf("CUDA_ERROR: load %g m %d jz %d %d ci0 %d %d\n", val, m, jy,
                 jz, ci0[1], ci0[2]);
        }
#endif
        (*this)(m, 0, jy + ci0[1], jz + ci0[2]) =
          _d_flds(m, 0, jy + ci0[1], jz + ci0[2]);
      }
      ti += THREADS_PER_BLOCK;
    }
  }

  __host__ __device__ float operator()(int m, int i, int j, int k) const
  {
    return data_[index(m, i, j, k)];
  }

private: // it's supposed to be a (read-only) cache, after all
  __host__ __device__ float& operator()(int m, int i, int j, int k)
  {
    return data_[index(m, i, j, k)];
  }

private:
  __host__ __device__ int index(int m, int i, int j, int k) const
  {
#ifdef DEBUG_FLDCACHE
    if (j < ci0y_ - 2 || j >= ci0y_ + BLOCKSIZE_Y + 2 || k < ci0z_ - 2 ||
        k >= ci0z_ + BLOCKSIZE_Z + 2) {
      printf("CUDA_ERROR: fld cache jk %d %d ci0 %d %d\n", j, k, ci0y_, ci0z_);
    }
#endif
    return (((m - EX) * (BLOCKSIZE_Z + 4) + k) * (BLOCKSIZE_Y + 4) + j) + off_;
  }

  float data_[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
  int off_;
#ifdef DEBUG_FLDCACHE
  int ci0y_, ci0z_;
#endif
};
