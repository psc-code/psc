
#include "cuda_mfields.h"
#include "cuda_bits.h"

#include "fields.hxx"

#include <cstdio>
#include <cassert>

// ======================================================================
// cuda_mfields

// ----------------------------------------------------------------------
// ctor

cuda_mfields::cuda_mfields(const Grid_t& grid, int _n_fields, const Int3& ibn)
  : ib(-ibn),
    im(grid.ldims + 2 * ibn),
    n_patches(grid.n_patches()),
    n_fields(_n_fields),
    n_cells_per_patch(im[0] * im[1] * im[2]),
    n_cells(n_patches * n_cells_per_patch),
    d_flds_(n_fields * n_cells),
    grid_(grid)
{
  cuda_base_init();
}

// ----------------------------------------------------------------------
// to_json

mrc_json_t cuda_mfields::to_json()
{
  mrc_json_t json = mrc_json_object_new(9);
  mrc_json_object_push_integer(json, "n_patches", n_patches);
  mrc_json_object_push_integer(json, "n_fields", n_fields);
  mrc_json_object_push_integer(json, "n_cells_per_patch", n_cells_per_patch);
  mrc_json_object_push_integer(json, "n_cells", n_cells);

  mrc_json_object_push(json, "ib", mrc_json_integer_array_new(3, ib));
  mrc_json_object_push(json, "im", mrc_json_integer_array_new(3, im));

  mrc_json_t json_flds = mrc_json_object_new(2);
  mrc_json_object_push(json, "flds", json_flds);
  mrc_json_object_push_boolean(json_flds, "__field5d__", true);
  mrc_json_t json_flds_patches = mrc_json_array_new(n_patches);
  mrc_json_object_push(json_flds, "data", json_flds_patches);

  fields_single_t flds = get_host_fields();
  Fields3d<fields_single_t> F(flds);
  for (int p = 0; p < n_patches; p++) {
    copy_from_device(p, flds, 0, n_fields);

    mrc_json_t json_flds_comps = mrc_json_array_new(n_fields);
    mrc_json_array_push(json_flds_patches, json_flds_comps);
    for (int m = 0; m < n_fields; m++) {
      mrc_json_t json_fld_z = mrc_json_array_new(im[2]);
      mrc_json_array_push(json_flds_comps, json_fld_z);
      for (int k = ib[2]; k < ib[2] + im[2]; k++) {
	mrc_json_t json_fld_y = mrc_json_array_new(im[1]);
	mrc_json_array_push(json_fld_z, json_fld_y);
	for (int j = ib[1]; j < ib[1] + im[1]; j++) {
	  mrc_json_t json_fld_x = mrc_json_array_new(im[0]);
	  mrc_json_array_push(json_fld_y, json_fld_x);
	  for (int i = ib[0]; i < ib[0] + im[0]; i++) {
	    mrc_json_array_push_double(json_fld_x, F(m, i,j,k));
	  }
	}
      }
    }
  }
  flds.dtor();

  return json;
}

// ----------------------------------------------------------------------
// dump

void cuda_mfields::dump(const char *filename)
{
  mrc_json_t json = to_json();

  const char *buf = mrc_json_to_string(json);
  if (filename) {
    FILE *file = fopen(filename, "w");
    assert(file);
    fwrite(buf, 1, strlen(buf), file);
    fclose(file);
  } else {
    printf("cuda_mfields (json):\n%s\n", buf);
  }
  free((void *) buf);

  // FIXME free json
}

// ----------------------------------------------------------------------
// cast to DMFields

cuda_mfields::operator DMFields()
{
  return DMFields(d_flds_.data().get(), n_cells_per_patch * n_fields, im, ib);
}

// ----------------------------------------------------------------------
// operator[]

DFields cuda_mfields::operator[](int p)
{
  return static_cast<DMFields>(*this)[p];
}

// ----------------------------------------------------------------------
// get_host_fields

fields_single_t cuda_mfields::get_host_fields()
{
  return fields_single_t(grid(), ib, im, n_fields);
}

// ----------------------------------------------------------------------
// copy_to_device

void cuda_mfields::copy_to_device(int p, fields_single_t h_flds, int mb, int me)
{
  cudaError_t ierr;
  
  if (mb == me) {
    return;
  }
  assert(mb < me);

  uint size = n_cells_per_patch;
  ierr = cudaMemcpy((*this)[p].data() + mb * size,
		    h_flds.data() + mb * size,
		    (me - mb) * size * sizeof(float),
		    cudaMemcpyHostToDevice); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// copy_from_device

void cuda_mfields::copy_from_device(int p, fields_single_t h_flds, int mb, int me)
{
  cudaError_t ierr;

  if (mb == me) {
    return;
  }
  assert(mb < me);

  uint size = n_cells_per_patch;
  ierr = cudaMemcpy(h_flds.data() + mb * size,
		    (*this)[p].data() + mb * size,
		    (me - mb) * size * sizeof(float),
		    cudaMemcpyDeviceToHost); cudaCheck(ierr);
}

#define BND (2) // FIXME
#define X3_DEV_OFF_YZ(fldnr, jy,jz)					\
  ((((fldnr)								\
     *mz + ((jz)+2))							\
    *my + ((jy)+2))							\
   *1 + (0))

#define F3_DDEV(d_flds, fldnr,ix,jy,jz)		\
  (d_flds)[X3_DEV_OFF_YZ(fldnr, jy,jz)]

#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

// ----------------------------------------------------------------------
// axpy_comp_yz

__global__ static void
k_axpy_comp_yz(float *y_flds, int ym, float a, float *x_flds, int xm,
	     int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (iy >= my || iz >= mz) {
    return;
  }

  iy -= BND;
  iz -= BND;

  F3_DDEV(y_flds, ym, 0,iy,iz) += a * F3_DDEV(x_flds, xm, 0,iy,iz);
}

void cuda_mfields::axpy_comp_yz(int ym, float a, cuda_mfields *cmflds_x, int xm)
{
  int my = im[1];
  int mz = im[2];
  assert(ib[1] == -BND);
  assert(ib[2] == -BND);

  dim3 dimGrid((my + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
	       (mz + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z);
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);

  for (int p = 0; p < n_patches; p++) {
    k_axpy_comp_yz<<<dimGrid, dimBlock>>>((*this)[p].data(), ym, a,
					  (*cmflds_x)[p].data(), xm, my, mz);
  }
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// zero_comp_yz

// FIXME, this should be more easily doable by just a cudaMemset()

__global__ static void
k_zero_comp_yz(float *x_flds, int xm, int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (iy >= my || iz >= mz) {
    return;
  }

  iy -= BND;
  iz -= BND;

  F3_DDEV(x_flds, xm, 0,iy,iz) = 0.f;
}

__global__ static void
k_zero_comp_xyz(float *data, uint n, uint stride)
{
  uint i = blockIdx.x * blockDim.x + threadIdx.x;
  uint p = blockIdx.y;

  if (i < n) {
    data[i + p * stride] = 0.f;
  }
}

void cuda_mfields::zero_comp(int m, dim_yz tag)
{
  int my = im[1];
  int mz = im[2];
  assert(ib[1] == -BND);
  assert(ib[2] == -BND);

  dim3 dimGrid((my + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
	       (mz + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z);
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);

  // OPT, should be done in a single kernel
  for (int p = 0; p < n_patches; p++) {
    k_zero_comp_yz<<<dimGrid, dimBlock>>>((*this)[p].data(), m, my, mz);
  }
  cuda_sync_if_enabled();
}

void cuda_mfields::zero_comp(int m, dim_xyz tag)
{
  int n = n_cells_per_patch;
  int stride = n * n_fields;

  const int THREADS_PER_BLOCK = 512;
  dim3 dimBlock(THREADS_PER_BLOCK);
  dim3 dimGrid((n + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, n_patches);

  k_zero_comp_xyz<<<dimGrid, dimBlock>>>(data() + m * n, n, stride);
  cuda_sync_if_enabled();
}

