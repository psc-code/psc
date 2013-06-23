
#include <mrc_params.h>
#include <mrc_io.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>

static const float a1 = -1.586134342;
static const float a2 = -0.05298011854;
static const float a3 = 0.8829110762;
static const float a4 = 0.4435068522;
    
static const float k1 = 0.81289306611596146; // 1/1.230174104914
static const float k2 = 0.61508705245700002; // 1.230174104914/2

static const float k1i = 1.230174104914;
static const float k2i = 1.6257861322319229;

static void
wt97_x(struct mrc_fld *f, int *dims)
{
  struct mrc_fld *tmp = mrc_fld_duplicate(f);

  for (int iz = 0; iz < dims[2]; iz++) {
    for (int iy = 0; iy < dims[1]; iy++) {
      // Core 1D lifting process in this loop.
        
      // Predict 1. y1
      for (int ix = 1; ix < dims[0] - 1; ix += 2) {
	MRC_F3(f,0, ix,iy,iz) += a1 * (MRC_F3(f, 0, ix-1,iy,iz) + MRC_F3(f, 0, ix+1,iy,iz));
      }
      MRC_F3(f,0, dims[0] - 1,iy,iz) += 2 * a1 * MRC_F3(f, 0, dims[0] - 2,iy,iz); // Symmetric extension
      
      // Update 1. y0
      for (int ix = 2; ix < dims[0]; ix += 2) {
	MRC_F3(f,0, ix,iy,iz) += a2 * (MRC_F3(f, 0, ix-1,iy,iz) + MRC_F3(f, 0, ix+1,iy,iz));
      }
      MRC_F3(f,0, 0,iy,iz) += 2 * a2 * MRC_F3(f, 0, 1,iy,iz); // Symmetric extension
        
      // Predict 2.
      for (int ix = 1; ix < dims[0] - 1; ix += 2) {
	MRC_F3(f,0, ix,iy,iz) += a3 * (MRC_F3(f, 0, ix-1,iy,iz) + MRC_F3(f, 0, ix+1,iy,iz));
      }
      MRC_F3(f,0, dims[0] - 1,iy,iz) += 2 * a3 * MRC_F3(f, 0, dims[0] - 2,iy,iz); // Symmetric extension
        
      // Update 2.
      for (int ix = 2; ix < dims[0]; ix += 2) {
	MRC_F3(f,0, ix,iy,iz) += a4 * (MRC_F3(f, 0, ix-1,iy,iz) + MRC_F3(f, 0, ix+1,iy,iz));
      }
      MRC_F3(f,0, 0,iy,iz) += 2 * a4 * MRC_F3(f, 0, 1,iy,iz); // Symmetric extension

      // de-interleave
      for (int ix = 0; ix < dims[0] / 2; ix++) {
	// k1 and k2 scale the vals
	MRC_F3(tmp, 0, ix,iy,iz) = k1 * MRC_F3(f, 0, 2*ix,iy,iz);
	MRC_F3(tmp, 0, ix + dims[0]/2,iy,iz) = k2 * MRC_F3(f, 0, 2*ix+1,iy,iz);
      }
    }
  }
  for (int iz = 0; iz < dims[2]; iz++) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f, 0, ix,iy,iz) = MRC_F3(tmp, 0, ix,iy,iz);
      }
    }
  }
  mrc_fld_destroy(tmp);
}

static void
wt97_y(struct mrc_fld *f, int *dims)
{
  struct mrc_fld *tmp = mrc_fld_duplicate(f);

  for (int iz = 0; iz < dims[2]; iz++) {
    // Core 1D lifting process in this loop.
        
    // Predict 1. y1
    for (int iy = 1; iy < dims[1] - 1; iy += 2) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) += a1 * (MRC_F3(f, 0, ix,iy-1,iz) + MRC_F3(f, 0, ix,iy+1,iz));
      }
    }
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,dims[1] - 1,iz) += 2 * a1 * MRC_F3(f, 0, ix,dims[1] - 2,iz); // Symmetric extension
    }
      
      // Update 1. y0
    for (int iy = 2; iy < dims[1]; iy += 2) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) += a2 * (MRC_F3(f, 0, ix,iy-1,iz) + MRC_F3(f, 0, ix,iy+1,iz));
      }
    }
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,0,iz) += 2 * a2 * MRC_F3(f, 0, ix,1,iz); // Symmetric extension
    }
        
    // Predict 2.
    for (int iy = 1; iy < dims[1] - 1; iy += 2) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) += a3 * (MRC_F3(f, 0, ix,iy-1,iz) + MRC_F3(f, 0, ix,iy-1,iz));
      }
    }
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,dims[1] - 1,iz) += 2 * a3 * MRC_F3(f, 0, ix,dims[1] - 2,iz); // Symmetric extension
    }
        
    // Update 2.
    for (int iy = 2; iy < dims[1]; iy += 2) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) += a4 * (MRC_F3(f, 0, ix,iy-1,iz) + MRC_F3(f, 0, ix,iy+1,iz));
      }
    }
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,0,iz) += 2 * a4 * MRC_F3(f, 0, ix,1,iz); // Symmetric extension
    }

    // de-interleave
    for (int iy = 0; iy < dims[1] / 2; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	// k1 and k2 scale the vals
	MRC_F3(tmp, 0, ix,iy,iz) = k1 * MRC_F3(f, 0, ix,2*iy,iz);
	MRC_F3(tmp, 0, ix,iy + dims[1]/2,iz) = k2 * MRC_F3(f, 0, ix,2*iy+1,iz);
      }
    }
  }
  for (int iz = 0; iz < dims[2]; iz++) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f, 0, ix,iy,iz) = MRC_F3(tmp, 0, ix,iy,iz);
      }
    }
  }
  mrc_fld_destroy(tmp);
}

static void
wt97_z(struct mrc_fld *f, int *dims)
{
  struct mrc_fld *tmp = mrc_fld_duplicate(f);

  // Predict 1. y1
  for (int iz = 1; iz < dims[2] - 1; iz += 2) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) += a1 * (MRC_F3(f, 0, ix,iy,iz-1) + MRC_F3(f, 0, ix,iy,iz+1));
      }
    }
  }
  for (int iy = 0; iy < dims[1]; iy++) {
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,iy,dims[2] - 1) += 2 * a1 * MRC_F3(f, 0, ix,iy,dims[2] - 2); // Symmetric extension
    }
  }
      
  // Update 1. y0
  for (int iz = 2; iz < dims[2]; iz += 2) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) += a2 * (MRC_F3(f, 0, ix,iy,iz-1) + MRC_F3(f, 0, ix,iy,iz+1));
      }
    }
  }
  for (int iy = 0; iy < dims[1]; iy++) {
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,iy,0) += 2 * a2 * MRC_F3(f, 0, ix,iy,1); // Symmetric extension
    }
  }

  // Predict 2.
  for (int iz = 1; iz < dims[2] - 1; iz += 2) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) += a3 * (MRC_F3(f, 0, ix,iy,iz-1) + MRC_F3(f, 0, ix,iy,iz+1));
      }
    }
  }
  for (int iy = 0; iy < dims[1]; iy++) {
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,iy,dims[2] - 1) += 2 * a3 * MRC_F3(f, 0, ix,iy,dims[2] - 2); // Symmetric extension
    }
  }
    
  // Update 2.
  for (int iz = 2; iz < dims[1]; iz += 2) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) += a4 * (MRC_F3(f, 0, ix,iy,iz-1) + MRC_F3(f, 0, ix,iy,iz+1));
      }
    }
  }
  for (int iy = 0; iy < dims[1]; iy++) {
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,iy,0) += 2 * a4 * MRC_F3(f, 0, ix,iy,1); // Symmetric extension
    }
  }

  // de-interleave
  for (int iz = 0; iz < dims[2] / 2; iz++) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	// k1 and k2 scale the vals
	MRC_F3(tmp, 0, ix,iy,iz) = k1 * MRC_F3(f, 0, ix,iy,2*iz);
	MRC_F3(tmp, 0, ix,iy,iz + dims[2]/2) = k2 * MRC_F3(f, 0, ix,iy,2*iz+1);
      }
    }
  }
  for (int iz = 0; iz < dims[2]; iz++) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f, 0, ix,iy,iz) = MRC_F3(tmp, 0, ix,iy,iz);
      }
    }
  }

  mrc_fld_destroy(tmp);
}

static void
iwt97_x(struct mrc_fld *f, int *dims)
{
  struct mrc_fld *tmp = mrc_fld_duplicate(f);

  for (int iz = 0; iz < dims[2]; iz++) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(tmp, 0, ix,iy,iz) = MRC_F3(f, 0, ix,iy,iz);
      }
    }
  }
  for (int iz = 0; iz < dims[2]; iz++) {
    for (int iy = 0; iy < dims[1]; iy++) {
      // interleave
      for (int ix = 0; ix < dims[0] / 2; ix++) {
	// ik1 and ik2 scale the vals back
	MRC_F3(f, 0, 2*ix  ,iy,iz) = k1i * MRC_F3(tmp, 0, ix,iy,iz);
	MRC_F3(f, 0, 2*ix+1,iy,iz) = k2i * MRC_F3(tmp, 0, ix + dims[0]/2,iy,iz);
      }

      // inverse 1D lifting
        
      // inverse update 2.
      for (int ix = 2; ix < dims[0]; ix += 2) {
	MRC_F3(f,0, ix,iy,iz) -= a4 * (MRC_F3(f, 0, ix-1,iy,iz) + MRC_F3(f, 0, ix+1,iy,iz));
      }
      MRC_F3(f,0, 0,iy,iz) -= 2 * a4 * MRC_F3(f, 0, 1,iy,iz); // Symmetric extension
        
      // inverse predict 2.
      for (int ix = 1; ix < dims[0] - 1; ix += 2) {
	MRC_F3(f,0, ix,iy,iz) -= a3 * (MRC_F3(f, 0, ix-1,iy,iz) + MRC_F3(f, 0, ix+1,iy,iz));
      }
      MRC_F3(f,0, dims[0] - 1,iy,iz) -= 2 * a3 * MRC_F3(f, 0, dims[0] - 2,iy,iz); // Symmetric extension

      // inverse update 1. y0
      for (int ix = 2; ix < dims[0]; ix += 2) {
	MRC_F3(f,0, ix,iy,iz) -= a2 * (MRC_F3(f, 0, ix-1,iy,iz) + MRC_F3(f, 0, ix+1,iy,iz));
      }
      MRC_F3(f,0, 0,iy,iz) -= 2 * a2 * MRC_F3(f, 0, 1,iy,iz); // Symmetric extension
        
      // inverse predict 1. y1
      for (int ix = 1; ix < dims[0] - 1; ix += 2) {
	MRC_F3(f,0, ix,iy,iz) -= a1 * (MRC_F3(f, 0, ix-1,iy,iz) + MRC_F3(f, 0, ix+1,iy,iz));
      }
      MRC_F3(f,0, dims[0] - 1,iy,iz) -= 2 * a1 * MRC_F3(f, 0, dims[0] - 2,iy,iz); // Symmetric extension
    }
  }
  mrc_fld_destroy(tmp);
}

static void
iwt97_y(struct mrc_fld *f, int *dims)
{
  struct mrc_fld *tmp = mrc_fld_duplicate(f);

  for (int iz = 0; iz < dims[2]; iz++) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(tmp, 0, ix,iy,iz) = MRC_F3(f, 0, ix,iy,iz);
      }
    }
  }
  for (int iz = 0; iz < dims[2]; iz++) {
    // interleave
    for (int iy = 0; iy < dims[1] / 2; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	// k1 and k2 scale the vals
	MRC_F3(f, 0, ix,2*iy  ,iz) = k1i * MRC_F3(tmp, 0, ix,iy,iz);
	MRC_F3(f, 0, ix,2*iy+1,iz) = k2i * MRC_F3(tmp, 0, ix,iy + dims[1]/2,iz);
      }
    }
    // inverse 1D lifting
        
    // inverse update 2.
    for (int iy = 2; iy < dims[1]; iy += 2) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) -= a4 * (MRC_F3(f, 0, ix,iy-1,iz) + MRC_F3(f, 0, ix,iy+1,iz));
      }
    }
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,0,iz) -= 2 * a4 * MRC_F3(f, 0, ix,1,iz); // Symmetric extension
    }

    // inverse predict 2.
    for (int iy = 1; iy < dims[1] - 1; iy += 2) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) -= a3 * (MRC_F3(f, 0, ix,iy-1,iz) + MRC_F3(f, 0, ix,iy-1,iz));
      }
    }
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,dims[1] - 1,iz) -= 2 * a3 * MRC_F3(f, 0, ix,dims[1] - 2,iz); // Symmetric extension
    }
        
    // inverse update 1. y0
    for (int iy = 2; iy < dims[1]; iy += 2) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) -= a2 * (MRC_F3(f, 0, ix,iy-1,iz) + MRC_F3(f, 0, ix,iy+1,iz));
      }
    }
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,0,iz) -= 2 * a2 * MRC_F3(f, 0, ix,1,iz); // Symmetric extension
    }
        
    // inverse predict 1. y1
    for (int iy = 1; iy < dims[1] - 1; iy += 2) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) -= a1 * (MRC_F3(f, 0, ix,iy-1,iz) + MRC_F3(f, 0, ix,iy+1,iz));
      }
    }
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,dims[1] - 1,iz) -= 2 * a1 * MRC_F3(f, 0, ix,dims[1] - 2,iz); // Symmetric extension
    }
      
  }
  mrc_fld_destroy(tmp);
}

static void
iwt97_z(struct mrc_fld *f, int *dims)
{
  struct mrc_fld *tmp = mrc_fld_duplicate(f);

  for (int iz = 0; iz < dims[2]; iz++) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(tmp, 0, ix,iy,iz) = MRC_F3(f, 0, ix,iy,iz);
      }
    }
  }
  // interleave
  for (int iz = 0; iz < dims[2] / 2; iz++) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	// k1 and k2 scale the vals
	MRC_F3(f, 0, ix,iy,2*iz  ) = k1i * MRC_F3(tmp, 0, ix,iy,iz);
	MRC_F3(f, 0, ix,iy,2*iz+1) = k2i * MRC_F3(tmp, 0, ix,iy,iz + dims[2]/2);
      }
    }
  }

  // inverse update 2.
  for (int iz = 2; iz < dims[1]; iz += 2) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) -= a4 * (MRC_F3(f, 0, ix,iy,iz-1) + MRC_F3(f, 0, ix,iy,iz+1));
      }
    }
  }
  for (int iy = 0; iy < dims[1]; iy++) {
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,iy,0) -= 2 * a4 * MRC_F3(f, 0, ix,iy,1); // Symmetric extension
    }
  }

  // inverse predict 2.
  for (int iz = 1; iz < dims[2] - 1; iz += 2) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) -= a3 * (MRC_F3(f, 0, ix,iy,iz-1) + MRC_F3(f, 0, ix,iy,iz+1));
      }
    }
  }
  for (int iy = 0; iy < dims[1]; iy++) {
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,iy,dims[2] - 1) -= 2 * a3 * MRC_F3(f, 0, ix,iy,dims[2] - 2); // Symmetric extension
    }
  }
  // inverse update 1. y0
  for (int iz = 2; iz < dims[2]; iz += 2) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) -= a2 * (MRC_F3(f, 0, ix,iy,iz-1) + MRC_F3(f, 0, ix,iy,iz+1));
      }
    }
  }
  for (int iy = 0; iy < dims[1]; iy++) {
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,iy,0) -= 2 * a2 * MRC_F3(f, 0, ix,iy,1); // Symmetric extension
    }
  }
        
  // inverse predict 1. y1
  for (int iz = 1; iz < dims[2] - 1; iz += 2) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	MRC_F3(f,0, ix,iy,iz) -= a1 * (MRC_F3(f, 0, ix,iy,iz-1) + MRC_F3(f, 0, ix,iy,iz+1));
      }
    }
  }
  for (int iy = 0; iy < dims[1]; iy++) {
    for (int ix = 0; ix < dims[0]; ix++) {
      MRC_F3(f,0, ix,iy,dims[2] - 1) -= 2 * a1 * MRC_F3(f, 0, ix,iy,dims[2] - 2); // Symmetric extension
    }
  }
      
  mrc_fld_destroy(tmp);
}

static void
wt97(struct mrc_fld *f, int nr_levels)
{
  int dims[3];
  mrc_domain_get_global_dims(f->_domain, dims);

  for (int l = 0; l < nr_levels; l++) {
    wt97_x(f, dims);
    wt97_y(f, dims);
    wt97_z(f, dims);
    for (int d = 0; d < 3; d++) {
      dims[d] /= 2;
      assert(dims[d] % 2 == 0);
    }
  }
}


static void
iwt97(struct mrc_fld *f, int nr_levels)
{
  int dims[3];
  mrc_domain_get_global_dims(f->_domain, dims);

  for (int l = 0; l < nr_levels - 1; l++) {
    for (int d = 0; d < 3; d++) {
      dims[d] /= 2;
      assert(dims[d] % 2 == 0);
    }
  }
  for (int l = 0; l < nr_levels; l++) {
    iwt97_x(f, dims);
    iwt97_y(f, dims);
    iwt97_z(f, dims);
    for (int d = 0; d < 3; d++) {
      dims[d] *= 2;
    }
  }
}

static void
threshold(struct mrc_fld *f, float eps)
{
  int cnt = 0;
  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    if (fabs(MRC_F3(f, 0, ix,iy,iz)) < eps) {
      MRC_F3(f, 0, ix,iy,iz) = 0.;
    } else {
      cnt++;
    }
  } mrc_fld_foreach_end;
  
  int dims[3];
  mrc_domain_get_global_dims(f->_domain, dims);
  int len = dims[0] * dims[1] * dims[2];
  mprintf("cnt = %d / %d = %g %%\n", cnt, len, cnt * 100. / len);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  int nr_levels = 3;
  float eps = 1e-9;
  mrc_params_get_option_int("nr_levels", &nr_levels);
  mrc_params_get_option_float("eps", &eps);

  struct mrc_io *io = mrc_io_create(MPI_COMM_WORLD);
  mrc_io_set_param_string(io, "basename", "hc0002.3df");
  mrc_io_set_param_string(io, "outdir", "/Users/kai/ggcm/hc0002/target");
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "r", 300, 0.);
  struct mrc_fld *fld = mrc_fld_read(io, "mrc_fld_pp");
  mrc_fld_set_comp_name(fld, 0, "pp"); // FIXME
  mrc_io_close(io);
  mrc_io_destroy(io);

  struct mrc_fld *fw = mrc_fld_duplicate(fld);
  mrc_fld_set_name(fw, "fw");
  mrc_fld_set_comp_name(fw, 0, "pp");
  mrc_fld_copy(fw, fld);
  wt97(fw, nr_levels);
  threshold(fw, eps);

  struct mrc_fld *iw = mrc_fld_duplicate(fld);
  mrc_fld_set_name(iw, "iw");
  mrc_fld_set_comp_name(iw, 0, "pp");
  mrc_fld_copy(iw, fw);
  iwt97(iw, nr_levels);

  io = mrc_io_create(MPI_COMM_WORLD);
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_fld_write(fld, io);
  mrc_io_close(io);
  mrc_io_open(io, "w", 1, 1.);
  mrc_fld_write(fw, io);
  mrc_io_close(io);
  mrc_io_open(io, "w", 2, 2.);
  mrc_fld_write(iw, io);
  mrc_io_close(io);
  mrc_io_destroy(io);

  mrc_fld_destroy(fw);
  
  mrc_fld_destroy(fld);

  MPI_Finalize();
  return 0;
}
