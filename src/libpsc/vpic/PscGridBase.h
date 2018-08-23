
#ifndef PSC_GRID_BASE_H
#define PSC_GRID_BASE_H

#include "psc_vpic_bits.h"
#include "vec3.hxx"

#include <array>
#include <vector>
#include <cassert>
#include <mrc_common.h>

// ======================================================================
// PscMp

struct PscMp
{
  PscMp(int n_port) :
    recv_buf_(n_port),
    send_buf_(n_port),
    recv_req_(n_port, MPI_REQUEST_NULL),
    send_req_(n_port, MPI_REQUEST_NULL)
  {
  }

  void size_recv_buffer(int port, int size)
  {
    recv_buf_[port].resize(size);
  }

  void size_send_buffer(int port, int size)
  {
    send_buf_[port].resize(size);
  }

  char* send_buffer(int port)
  {
    return send_buf_[port].data();
  }

  char* recv_buffer(int port)
  {
    return recv_buf_[port].data();
  }

  void begin_send(int port, int size, int dst, int tag)
  {
    MPI_Isend(send_buf_[port].data(), size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, &send_req_[port]);
  }

  void end_send(int port)
  {
    MPI_Wait(&send_req_[port], MPI_STATUS_IGNORE);
  }

  void begin_recv(int port, int size, int src, int tag)
  {
    MPI_Irecv(recv_buf_[port].data(), size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &recv_req_[port]);
  }

  void end_recv(int port)
  {
    MPI_Wait(&recv_req_[port], MPI_STATUS_IGNORE);
  }

  std::vector<std::vector<char>> recv_buf_;
  std::vector<std::vector<char>> send_buf_;
  std::vector<MPI_Request> recv_req_;
  std::vector<MPI_Request> send_req_;
};

// ======================================================================
// PscGridBase

struct PscGridBase
{
  enum {
    anti_symmetric_fields = -1, // E_tang = 0
    pec_fields            = -1, // FIXME, matching vpic for now
    symmetric_fields      = -2, // B_tang = 0, B_norm = 0
    pmc_fields            = -3, // B_tang = 0, B_norm floats
    absorb_fields         = -4, 
  } Fbc;

  enum {
    reflect_particles = -1, // FIXME, matching vpic for now
    absorb_particles  = -2 
  } Pbc;

  PscGridBase()
  {
    memset(this, 0, sizeof(*this)); // FIXME
    neighbor = nullptr; // FIXME...
    for(int i = 0; i < 27; i++) {
      bc[i] = anti_symmetric_fields;
    }
    bc[BOUNDARY(0,0,0)] = psc_world_rank;
    mp = new PscMp(27);
  }
  
  ~PscGridBase()
  {
    delete[] neighbor;
    delete[] range;
    delete mp;
  }
  
  static PscGridBase* create()
  {
    return new PscGridBase;
  }

  void setup(double dx_[3], double dt_, double cvac_, double eps0_)
  {
    dx = dx_[0]; dy = dx_[1]; dz = dx_[2];
    dt = dt_;
    cvac = cvac_;
    eps0 = eps0_;
  }

  Int3 rank_to_idx(int rank, const Int3& np)
  {
    Int3 idx;
    int ix,iy,iz;
    ix = rank;
    iy = ix / np[0];
    ix -= iy * np[0];
    iz = iy / np[1];
    iy -= iz * np[1];
    idx[0] = ix;
    idx[1] = iy;
    idx[2] = iz;
    return idx;
  }

  unsigned int idx_to_rank(const int np[3], Int3 idx)
  {
    int ix = idx[0], iy = idx[1], iz = idx[2];
    while (ix >= np[0]) ix -= np[0];
    while (ix < 0) ix += np[0];
    while (iy >= np[1]) iy -= np[1];
    while (iy < 0) iy += np[1];
    while (iz >= np[2]) iz -= np[2];
    while (iz < 0) iz += np[2];
    return ix + np[0] * (iy + np[1] * iz);
  }
  
#define LOCAL_CELL_ID(x,y,z)  VOXEL(x,y,z, lnx,lny,lnz)
#define REMOTE_CELL_ID(x,y,z) VOXEL(x,y,z, rnx,rny,rnz)

  void size_grid(int lnx, int lny, int lnz)
  {
    int64_t x,y,z;
    int i, j, k;
    int64_t ii, jj, kk; 

    assert(lnx > 0 && lny > 0 && lnz > 0);

    // Setup phase 2 data structures
    sx = 1;
    sy = (lnx+2) * sx;
    sz = (lny+2) * sy;
    nv = (lnz+2) * sz;
    nx = lnx; ny = lny; nz = lnz;
    for (k = -1; k <= 1; k++) {
      for (j = -1; j <= 1; j++) {
	for (i = -1; i <= 1; i++) { 
	  bc[BOUNDARY(i,j,k)] = pec_fields;
	}
      }
    }
    bc[BOUNDARY(0,0,0)] = psc_world_rank;

    // Setup phase 3 data structures.  This is an ugly kludge to
    // interface phase 2 and phase 3 data structures
    delete[] range;
    range = new int64_t[psc_world_size + 1];
    ii = nv; // nv is not 64-bits
    MPI_Allgather(&ii, 1, MPI_LONG_LONG, range, 1, MPI_LONG_LONG, MPI_COMM_WORLD);
    jj = 0;
    range[psc_world_size] = 0;
    for(i = 0; i <= psc_world_size; i++) {
      kk = range[i];
      range[i] = jj;
      jj += kk;
    }
    rangel = range[psc_world_rank];
    rangeh = range[psc_world_rank+1]-1;

    delete[] neighbor;
    neighbor = new int64_t[6 * nv];

    for (z = 0; z <= lnz+1; z++) {
      for (y = 0; y <= lny+1; y++) {
	for (x = 0; x <= lnx+1; x++) {
	  i = 6*LOCAL_CELL_ID(x,y,z);
	  neighbor[i+0] = rangel + LOCAL_CELL_ID(x-1,y,z);
	  neighbor[i+1] = rangel + LOCAL_CELL_ID(x,y-1,z);
	  neighbor[i+2] = rangel + LOCAL_CELL_ID(x,y,z-1);
	  neighbor[i+3] = rangel + LOCAL_CELL_ID(x+1,y,z);
	  neighbor[i+4] = rangel + LOCAL_CELL_ID(x,y+1,z);
	  neighbor[i+5] = rangel + LOCAL_CELL_ID(x,y,z+1);
	  // Set boundary faces appropriately
	  if (x == 1  ) neighbor[i+0] = reflect_particles;
	  if (y == 1  ) neighbor[i+1] = reflect_particles;
	  if (z == 1  ) neighbor[i+2] = reflect_particles;
	  if (x == lnx) neighbor[i+3] = reflect_particles;
	  if (y == lny) neighbor[i+4] = reflect_particles;
	  if (z == lnz) neighbor[i+5] = reflect_particles;
	  // Set ghost cells appropriately
	  if (x==0 || x==lnx+1 || y==0 || y==lny+1 || z==0 || z==lnz+1) {
	    neighbor[i+0] = reflect_particles;
	    neighbor[i+1] = reflect_particles;
	    neighbor[i+2] = reflect_particles;
	    neighbor[i+3] = reflect_particles;
	    neighbor[i+4] = reflect_particles;
	    neighbor[i+5] = reflect_particles;
	  }
	}
      }
    }
  }

  void join_grid(int boundary, int rank)
  {
    int lx, ly, lz, lnx, lny, lnz, rx, ry, rz, rnx, rny, rnz, rnc;

    assert(boundary >= 0 && boundary < 27 && boundary != BOUNDARY(0,0,0));
    assert(rank >= 0 && rank < psc_world_size);

    // Join phase 2 data structures
    bc[boundary] = rank;

    // Join phase 3 data structures
    lnx = nx;
    lny = ny;
    lnz = nz;
    rnc = range[rank+1] - range[rank]; // Note: rnc <~ 2^31 / 6

# define GLUE_FACE(tag,i,j,k,X,Y,Z) BEGIN_PRIMITIVE {			\
      if (boundary == BOUNDARY(i,j,k)) {				\
	assert(rnc % ((ln##Y+2)*(ln##Z+2)) == 0);			\
	rn##X = (rnc / ((ln##Y+2) * (ln##Z+2))) - 2;			\
	rn##Y = ln##Y;							\
	rn##Z = ln##Z;							\
	l##X = (i+j+k)<0 ? 1     : ln##X;				\
	r##X = (i+j+k)<0 ? rn##X : 1;					\
	for (l##Z = 1; l##Z <= ln##Z; l##Z++) {				\
	  for (l##Y = 1; l##Y <= ln##Y; l##Y++) {			\
	    r##Y = l##Y;						\
	    r##Z = l##Z;						\
	    neighbor[6*LOCAL_CELL_ID(lx,ly,lz) + tag] =			\
	      range[rank] + REMOTE_CELL_ID(rx,ry,rz);			\
	  }								\
	}								\
	return;								\
      }									\
    } END_PRIMITIVE

    GLUE_FACE(0,-1, 0, 0,x,y,z);
    GLUE_FACE(1, 0,-1, 0,y,z,x);
    GLUE_FACE(2, 0, 0,-1,z,x,y);
    GLUE_FACE(3, 1, 0, 0,x,y,z);
    GLUE_FACE(4, 0, 1, 0,y,z,x);
    GLUE_FACE(5, 0, 0, 1,z,x,y);
# undef GLUE_FACE
  }

#undef LOCAL_CELL_ID
#undef REMOTE_CELL_ID

  void partition_periodic_box(const double xl[3], const double xh[3],
			      const int gdims[3], const Int3 np)
  {
    assert(np[0] > 0 && np[1] > 0 && np[2] > 0);
    assert(gdims[0] > 0 && gdims[1] > 0 && gdims[2] > 0);
    assert(gdims[0] % np[0] == 0 && gdims[1] % np[1] == 0 && gdims[2] % np[2] == 0);

    Int3 idx = rank_to_idx(psc_world_rank, np);

    double dxyz[3];
    for (int d = 0; d < 3; d++) {
      dxyz[d] = (xh[d] - xl[d]) / gdims[d];
    }
    dx = dxyz[0]; dy = dxyz[1]; dz = dxyz[2];
    dV = dxyz[0] * dxyz[1] * dxyz[2];

    rdx = 1./dxyz[0]; rdy = 1./dxyz[1]; rdz = 1./dxyz[2];
    r8V = 0.125 / (dxyz[0] * dxyz[1] * dxyz[2]);

    double xyz0[3], xyz1[3];
    for (int d = 0; d < 3; d++) {
      double f;
      f = (double) (idx[d]  ) / np[d]; xyz0[d] = xl[d] * (1-f) + xh[d] * f;
      f = (double) (idx[d]+1) / np[d]; xyz1[d] = xl[d] * (1-f) + xh[d] * f;
    }
    x0 = xyz0[0]; y0 = xyz0[1]; z0 = xyz0[2];
    x1 = xyz1[0]; y1 = xyz1[1]; z1 = xyz1[2];
    
    // Size the local grid
    size_grid(gdims[0] / np[0], gdims[1] / np[1], gdims[2] / np[2]);

    // Join the grid to neighbors
    int px = idx[0], py = idx[1], pz = idx[2];
    join_grid(BOUNDARY(-1, 0, 0), idx_to_rank(np, {px-1, py  , pz  }));
    join_grid(BOUNDARY( 0,-1, 0), idx_to_rank(np, {px  , py-1, pz  }));
    join_grid(BOUNDARY( 0, 0,-1), idx_to_rank(np, {px  , py  , pz-1}));
    join_grid(BOUNDARY( 1, 0, 0), idx_to_rank(np, {px+1, py  , pz  }));
    join_grid(BOUNDARY( 0, 1, 0), idx_to_rank(np, {px  , py+1, pz  }));
    join_grid(BOUNDARY( 0, 0, 1), idx_to_rank(np, {px  , py  , pz+1}));
  }

  void set_fbc(int boundary, int fbc)
  {
    assert(boundary >= 0 && boundary < 27 && boundary != BOUNDARY(0,0,0) &&
	   (fbc == anti_symmetric_fields || fbc == symmetric_fields ||
	    fbc == pmc_fields            || fbc == absorb_fields    ));
    
    bc[boundary] = fbc;
  }

  void set_pbc(int boundary, int pbc)
  {
    assert(boundary >= 0 && boundary < 27 && boundary != BOUNDARY(0,0,0) &&
	   pbc < 0);
    
#define LOCAL_CELL_ID(x,y,z)  VOXEL(x,y,z, nx,ny,nz)
#define SET_PBC(tag,i,j,k,X,Y,Z) BEGIN_PRIMITIVE {			\
      if (boundary == BOUNDARY(i,j,k)) {				\
	int l##X = (i+j+k)<0 ? 1 : n##X;				\
	for(int l##Z=1; l##Z <= n##Z; l##Z++ )				\
	  for(int l##Y=1; l##Y <= n##Y; l##Y++ )			\
	    neighbor[6*LOCAL_CELL_ID(lx,ly,lz) + tag] = pbc;		\
	return;								\
      }									\
    } END_PRIMITIVE
    
    SET_PBC(0,-1, 0, 0,x,y,z);
    SET_PBC(1, 0,-1, 0,y,z,x);
    SET_PBC(2, 0, 0,-1,z,x,y);
    SET_PBC(3, 1, 0, 0,x,y,z);
    SET_PBC(4, 0, 1, 0,y,z,x);
    SET_PBC(5, 0, 0, 1,z,x,y);
    
#undef SET_PBC
#undef LOCAL_CELL_ID
  }

  void mp_size_recv_buffer(int port, int size)
  {
    mp->size_recv_buffer(port, size);
  }

  void mp_size_send_buffer(int port, int size)
  {
    mp->size_send_buffer(port, size);
  }

  void* mp_recv_buffer(int port)
  {
    return mp->recv_buffer(port);
  }

  void* mp_send_buffer(int port)
  {
    return mp->send_buffer(port);
  }

  void mp_begin_recv(int port, int size, int src, int tag)
  {
    mp->begin_recv(port, size, src, tag);
  }

  void mp_begin_send(int port, int size, int dst, int tag)
  {
    mp->begin_send(port, size, dst, tag);
  }

  void mp_end_recv(int port)
  {
    mp->end_recv(port);
  }

  void mp_end_send(int port)
  {
    mp->end_send(port);
  }

  // ----------------------------------------------------------------------
  // for field communications
  
  void* size_send_port(int i, int j, int k, int size)
  {
    int port = BOUNDARY(i, j, k), dst = bc[port];
    if (dst < 0 || dst >= psc_world_size) {
      return nullptr;
    }
    mp_size_send_buffer(port, size);
    return mp_send_buffer(port);
  }

  void begin_send_port(int i, int j, int k, int size)
  {
    int port = BOUNDARY(i, j, k), dst = bc[port];
    if (dst < 0 || dst >= psc_world_size) {
      return;
    }
    mp_begin_send(port, size, dst, port);
  }

  void end_send_port(int i, int j, int k)
  {
    int port = BOUNDARY(i, j, k), dst = bc[port];
    if (dst < 0 || dst >= psc_world_size) {
      return;
    }
    mp_end_send(port);
  }

  void begin_recv_port(int i, int j, int k, int size)
  {
    int port = BOUNDARY(-i,-j,-k), src = bc[port];
    if (src < 0 || src >= psc_world_size) {
      return;
    }
    mp_size_recv_buffer(port, size);
    mp_begin_recv(port, size, src, BOUNDARY(i,j,k));
  }
  
  void* end_recv_port(int i, int j, int k)
  {
    int port = BOUNDARY(-i,-j,-k), src = bc[port];
    if (src < 0 || src >= psc_world_size) {
      return nullptr;
    }
    mp_end_recv(port);
    return mp_recv_buffer(port);
  }

  // System of units
  float dt, cvac, eps0;

  // Time stepper.  The simulation time is given by
  // t = g->t0 + (double)g->dt*(double)g->step
  int64_t step;             // Current timestep
  double t0;                // Simulation time corresponding to step 0

  // Phase 2 grid data structures 
  float x0, y0, z0;         // Min corner local domain (must be coherent)
  float x1, y1, z1;         // Max corner local domain (must be coherent)
  int   nx, ny, nz;         // Local voxel mesh resolution.  Voxels are
                            // indexed FORTRAN style 0:nx+1,0:ny+1,0:nz+1
                            // with voxels 1:nx,1:ny,1:nz being non-ghost
                            // voxels.
  float dx, dy, dz, dV;     // Cell dimensions and volume (CONVENIENCE ...
                            // USE x0,x1 WHEN DECIDING WHICH NODE TO USE!)
  float rdx, rdy, rdz, r8V; // Inverse voxel dimensions and one over
                            // eight times the voxel volume (CONVENIENCE)
  int   sx, sy, sz, nv;     // Voxel indexing x-, y-,z- strides and the
                            // number of local voxels (including ghosts,
                            // (nx+2)(ny+2)(nz+2)), (CONVENIENCE)
  int   bc[27];             // (-1:1,-1:1,-1:1) FORTRAN indexed array of
                            // boundary conditions to apply at domain edge
                            // 0 ... nproc-1 ... comm boundary condition
                            // <0 ... locally applied boundary condition

  // Phase 3 grid data structures
  // NOTE: VOXEL INDEXING LIMITS NUMBER OF VOXELS TO 2^31 (INCLUDING
  // GHOSTS) PER NODE.  NEIGHBOR INDEXING FURTHER LIMITS TO
  // (2^31)/6.  BOUNDARY CONDITION HANDLING LIMITS TO 2^28 PER NODE
  // EMITTER COMPONENT ID INDEXING FURTHER LIMITS TO 2^26 PER NODE.
  // THE LIMIT IS 2^63 OVER ALL NODES THOUGH.
  int64_t * ALIGNED(16) range;
                          // (0:nproc) indexed array giving range of
                          // global indexes of voxel owned by each
                          // processor.  Replicated on each processor.
                          // (range[rank]:range[rank+1]-1) are global
                          // voxels owned by processor "rank".  Note:
                          // range[rank+1]-range[rank] <~ 2^31 / 6

  int64_t * ALIGNED(128) neighbor;
                          // (0:5,0:local_num_voxel-1) FORTRAN indexed
                          // array neighbor(0:5,lidx) are the global
                          // indexes of neighboring voxels of the
                          // voxel with local index "lidx".  Negative
                          // if neighbor is a boundary condition.

  int64_t rangel, rangeh; // Redundant for move_p performance reasons:
                          //   rangel = range[rank]
                          //   rangeh = range[rank+1]-1.
                          // Note: rangeh-rangel <~ 2^26

  // Nearest neighbor communications ports
  PscMp* mp;
};

#endif
