
#ifndef PSC_GRID_BASE_H
#define PSC_GRID_BASE_H

#include "psc_vpic_bits.h"

#include <cassert>
#include <mrc_common.h>

// ======================================================================
// PscMp

struct PscMp
{
  PscMp(int n_port_)
  {
    n_port = n_port_;
    MALLOC(rbuf,    n_port); MALLOC(sbuf,    n_port);
    MALLOC(rbuf_sz, n_port); MALLOC(sbuf_sz, n_port); 
    MALLOC(rreq_sz, n_port); MALLOC(sreq_sz, n_port); 
    MALLOC(rreq,    n_port); MALLOC(sreq,    n_port); 
    CLEAR( rbuf,    n_port); CLEAR( sbuf,    n_port);
    CLEAR( rbuf_sz, n_port); CLEAR( sbuf_sz, n_port); 
    CLEAR( rreq_sz, n_port); CLEAR( sreq_sz, n_port); 
    CLEAR( rreq,    n_port); CLEAR( sreq,    n_port); 
  }

  ~PscMp()
  {
    for(int port=0; port < n_port; port++) {
      FREE_ALIGNED(rbuf[port]); FREE_ALIGNED(sbuf[port]); 
    }
    FREE(rreq   ); FREE(sreq   ); 
    FREE(rreq_sz); FREE(sreq_sz); 
    FREE(rbuf_sz); FREE(sbuf_sz); 
    FREE(rbuf   ); FREE(sbuf   );
  }
  
  static PscMp* create(int n_port)
  {
    return new PscMp(n_port);
  }

  static void destroy(PscMp* mp)
  {
    delete mp;
  }

  constexpr static float RESIZE_FACTOR = 1.3125;
  
  void size_recv_buffer(int port, int size)
  {
    assert(port >= 0 && port < n_port);
    if (size <= rbuf_sz[port]) {
      return;
    }

    size *= RESIZE_FACTOR;
    FREE_ALIGNED(rbuf[port]);
    MALLOC_ALIGNED(rbuf[port], size, 128);
    rbuf_sz[port] = size;
  }

  void size_send_buffer(int port, int size)
  {
    assert(port >= 0 && port < n_port);
    if (size <= sbuf_sz[port]) {
      return;
    }

    size *= RESIZE_FACTOR;
    FREE_ALIGNED(sbuf[port]);
    MALLOC_ALIGNED(sbuf[port], size, 128);
    sbuf_sz[port] = size;
  }

  void* send_buffer(int port)
  {
    assert(port >= 0 && port < n_port);
    return sbuf[port];
  }

  void* recv_buffer(int port)
  {
    assert(port >= 0 && port < n_port);
    return rbuf[port];
  }

  void begin_send(int port, int size, int dst, int tag)
  {
    assert(port >= 0 && port < n_port);
    sreq_sz[port] = size;
    MPI_Isend(sbuf[port], size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, &sreq[port]);
  }

  void end_send(int port)
  {
    assert(port >= 0 && port < n_port);
    MPI_Wait(&sreq[port], MPI_STATUS_IGNORE);
  }

  void begin_recv(int port, int size,  int src, int tag)
  {
    assert(port >= 0 && port < n_port);
    rreq_sz[port] = size;
    MPI_Irecv(rbuf[port], size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &rreq[port]);
  }

  void end_recv(int port)
  {
    assert(port >= 0 && port < n_port);
    MPI_Wait(&rreq[port], MPI_STATUS_IGNORE);
  }

  int n_port;
  char * ALIGNED(128) * rbuf;
  char * ALIGNED(128) * sbuf;
  int * rbuf_sz;
  int * sbuf_sz;
  int * rreq_sz;
  int * sreq_sz;
  MPI_Request * rreq;
  MPI_Request * sreq;
};

typedef PscMp Mp;

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
    bc[BOUNDARY(0,0,0)] = psc_world_rank;
    mp = reinterpret_cast<mp_t*>(Mp::create(27));
  }
  
  ~PscGridBase()
  {
    FREE_ALIGNED(neighbor);
    FREE_ALIGNED(range);
    Mp::destroy(reinterpret_cast<Mp*>(mp));
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

  void mp_size_recv_buffer(int port, int size)
  {
    Mp* mp = reinterpret_cast<Mp*>(this->mp);
    mp->size_recv_buffer(port, size);
  }

  void mp_size_send_buffer(int port, int size)
  {
    Mp* mp = reinterpret_cast<Mp*>(this->mp);
    mp->size_send_buffer(port, size);
  }

  void* mp_recv_buffer(int port)
  {
    Mp* mp = reinterpret_cast<Mp*>(this->mp);
    return mp->recv_buffer(port);
  }

  void* mp_send_buffer(int port)
  {
    Mp* mp = reinterpret_cast<Mp*>(this->mp);
    return mp->send_buffer(port);
  }

  void mp_begin_recv(int port, int size, int src, int tag)
  {
    Mp* mp = reinterpret_cast<Mp*>(this->mp);
    mp->begin_recv(port, size, src, tag);
  }

  void mp_begin_send(int port, int size, int dst, int tag)
  {
    Mp* mp = reinterpret_cast<Mp*>(this->mp);
    mp->begin_send(port, size, dst, tag);
  }

  void mp_end_recv(int port)
  {
    Mp* mp = reinterpret_cast<Mp*>(this->mp);
    mp->end_recv(port);
  }

  void mp_end_send(int port)
  {
    Mp* mp = reinterpret_cast<Mp*>(this->mp);
    mp->end_send(port);
  }

  // ----------------------------------------------------------------------
  // for field communications
  
  void* size_send_port(int i, int j, int k, int size)
  {
    int port = BOUNDARY(i, j, k), dst = bc[port];
    if (dst < 0 || dst >= psc_world_size ) {
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
};

#endif
