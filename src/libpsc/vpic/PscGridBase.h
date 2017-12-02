
#ifndef PSC_GRID_BASE_H
#define PSC_GRID_BASE_H

#include <vpic.h>

#include <cassert>
#include <mrc_common.h>

// ======================================================================
// VpicMp

struct mp {
  int n_port;
  char * ALIGNED(128) * rbuf; char * ALIGNED(128) * sbuf;
  int * rbuf_sz;              int * sbuf_sz;
  int * rreq_sz;              int * sreq_sz;
  MPI_Request * rreq;         MPI_Request * sreq;
};

struct VpicMp : mp_t
{
  static VpicMp* create(int n_port)
  {
    return reinterpret_cast<VpicMp*>(::new_mp(n_port));
  }

  static void destroy(VpicMp* mp)
  {
    ::delete_mp(mp);
  }

  void size_recv_buffer(int port, int size)
  {
    ::mp_size_recv_buffer(this, port, size);
  }

  void size_send_buffer(int port, int size)
  {
    ::mp_size_send_buffer(this, port, size);
  }

  void* send_buffer(int port)
  {
    return ::mp_send_buffer(this, port);
  }

  void* recv_buffer(int port)
  {
    return ::mp_recv_buffer(this, port);
  }

  void begin_send(int port, int size,  int dst, int tag)
  {
    return ::mp_begin_send(this, port, size, dst, tag);
  }

  void end_send(int port)
  {
    return ::mp_end_send(this, port);
  }

  void begin_recv(int port, int size,  int src, int tag)
  {
    return ::mp_begin_recv(this, port, size, src, tag);
  }

  void end_recv(int port)
  {
    return ::mp_end_recv(this, port);
  }
};

// ======================================================================
// PscGridBase

struct PscGridBase : grid_t
{
  PscGridBase()
  {
    CLEAR(this, 1 );
    for(int i = 0; i < 27; i++) {
      bc[i] = anti_symmetric_fields;
    }
    bc[BOUNDARY(0,0,0)] = world_rank;
    mp = VpicMp::create(27);
  }
  
  ~PscGridBase()
  {
    FREE_ALIGNED(neighbor);
    FREE_ALIGNED(range);
    VpicMp::destroy(reinterpret_cast<VpicMp*>(mp));
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

  void partition_periodic_box(double xl[3], double xh[3], int gdims[3], int np[3])
  {
    ::partition_periodic_box(this, xl[0], xl[1], xl[2], xh[0], xh[1], xh[2],
			     gdims[0], gdims[1], gdims[2], np[0], np[1], np[2]);
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

  void mp_size_recv_buffer(int tag, int size)
  {
    VpicMp* mp = reinterpret_cast<VpicMp*>(this->mp);
    mp->size_recv_buffer(tag, size);
  }

  void mp_size_send_buffer(int tag, int size)
  {
    VpicMp* mp = reinterpret_cast<VpicMp*>(this->mp);
    mp->size_send_buffer(tag, size);
  }

  void* size_send_port(int i, int j, int k, int size)
  {
    VpicMp* mp = reinterpret_cast<VpicMp*>(this->mp);
    int port = BOUNDARY(i, j, k), dst = bc[port];
    if (dst < 0 || dst >= world_size ) {
      return nullptr;
    }
    mp->size_send_buffer(port, size);
    return mp->send_buffer(port);
  }

  void begin_send_port(int i, int j, int k, int size)
  {
    VpicMp* mp = reinterpret_cast<VpicMp*>(this->mp);
    int port = BOUNDARY(i, j, k), dst = bc[port];
    if (dst < 0 || dst >= world_size) {
      return;
    }
    mp->begin_send(port, size, dst, port);
  }

  void end_send_port(int i, int j, int k)
  {
    VpicMp* mp = reinterpret_cast<VpicMp*>(this->mp);
    int port = BOUNDARY(i, j, k), dst = bc[port];
    if (dst < 0 || dst >= world_size) {
      return;
    }
    mp->end_send(port);
  }

  void begin_recv_port(int i, int j, int k, int size)
  {
    VpicMp* mp = reinterpret_cast<VpicMp*>(this->mp);
    int port = BOUNDARY(-i,-j,-k), src = bc[port];
    if (src < 0 || src >= world_size) {
      return;
    }
    mp->size_recv_buffer(port, size);
    mp->begin_recv(port, size, src, BOUNDARY(i,j,k));
  }
  
  void* end_recv_port(int i, int j, int k)
  {
    VpicMp* mp = reinterpret_cast<VpicMp*>(this->mp);
    int port = BOUNDARY(-i,-j,-k), src = bc[port];
    if (src < 0 || src >= world_size) {
      return nullptr;
    }
    mp->end_recv(port);
    return mp->recv_buffer(port);
  }
};

#endif
