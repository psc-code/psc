
#ifndef PSC_GRID_BASE_H
#define PSC_GRID_BASE_H

#include "psc_vpic_bits.h"

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

struct PscGridBase : grid_t
{
  PscGridBase()
  {
    CLEAR(this, 1);
    for(int i = 0; i < 27; i++) {
      bc[i] = anti_symmetric_fields;
    }
    bc[BOUNDARY(0,0,0)] = psc_world_rank;
    mp = reinterpret_cast<mp_t*>(new PscMp(27));
  }
  
  ~PscGridBase()
  {
    FREE_ALIGNED(neighbor);
    FREE_ALIGNED(range);
    delete reinterpret_cast<PscMp*>(mp);
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
    PscMp* mp = reinterpret_cast<PscMp*>(this->mp);
    mp->size_recv_buffer(port, size);
  }

  void mp_size_send_buffer(int port, int size)
  {
    PscMp* mp = reinterpret_cast<PscMp*>(this->mp);
    mp->size_send_buffer(port, size);
  }

  void* mp_recv_buffer(int port)
  {
    PscMp* mp = reinterpret_cast<PscMp*>(this->mp);
    return mp->recv_buffer(port);
  }

  void* mp_send_buffer(int port)
  {
    PscMp* mp = reinterpret_cast<PscMp*>(this->mp);
    return mp->send_buffer(port);
  }

  void mp_begin_recv(int port, int size, int src, int tag)
  {
    PscMp* mp = reinterpret_cast<PscMp*>(this->mp);
    mp->begin_recv(port, size, src, tag);
  }

  void mp_begin_send(int port, int size, int dst, int tag)
  {
    PscMp* mp = reinterpret_cast<PscMp*>(this->mp);
    mp->begin_send(port, size, dst, tag);
  }

  void mp_end_recv(int port)
  {
    PscMp* mp = reinterpret_cast<PscMp*>(this->mp);
    mp->end_recv(port);
  }

  void mp_end_send(int port)
  {
    PscMp* mp = reinterpret_cast<PscMp*>(this->mp);
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
