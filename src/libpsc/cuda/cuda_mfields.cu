
#include "cuda_mfields.h"
#include "cuda_mfields_const.h"
#include "cuda_bits.h"

#include "fields.hxx"

#include <cstdio>
#include <cassert>

// ======================================================================
// cuda_mfields

// ----------------------------------------------------------------------
// ctor

cuda_mfields::cuda_mfields(Grid_t& grid, int _n_fields, const Int3& ibn)
  : n_patches(grid.patches.size())
{
  n_fields = _n_fields;
  ldims = grid.ldims;
  for (int d = 0; d < 3; d++) {
    ib[d] = -ibn[d];
    im[d] = ldims[d] + 2 * ibn[d];
    dx[d] = grid.dx[d];
  }

  n_cells_per_patch = im[0] * im[1] * im[2];
  n_cells = n_patches * n_cells_per_patch;

  d_flds_.resize(n_fields * n_cells);
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
  mrc_json_object_push(json, "ldims", mrc_json_integer_array_new(3, ldims));
  double dx[3] = { this->dx[0], this->dx[1], this->dx[2] };
  mrc_json_object_push(json, "dx", mrc_json_double_array_new(3, dx));

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
// get_host_fields

fields_single_t cuda_mfields::get_host_fields()
{
  return fields_single_t(ib, im, n_fields);
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
  ierr = cudaMemcpy(DMFields(this)[p].d_flds() + mb * size,
		    h_flds.data + mb * size,
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
  ierr = cudaMemcpy(h_flds.data + mb * size,
		    DMFields(this)[p].d_flds() + mb * size,
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

  dim3 dimGrid((ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
	       (ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z);
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);

  for (int p = 0; p < n_patches; p++) {
    k_axpy_comp_yz<<<dimGrid, dimBlock>>>(DMFields(this)[p].d_flds(), ym, a,
					  DMFields(cmflds_x)[p].d_flds(), xm, my, mz);
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

void cuda_mfields::zero_comp_yz(int xm)
{
  int my = im[1];
  int mz = im[2];

  dim3 dimGrid((ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
	       (ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z);
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);

  // OPT, should be done in a single kernel
  for (int p = 0; p < n_patches; p++) {
    k_zero_comp_yz<<<dimGrid, dimBlock>>>(DMFields(this)[p].d_flds(), xm, my, mz);
  }
  cuda_sync_if_enabled();
}

