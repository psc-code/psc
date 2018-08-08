
#ifndef VPIC_GRID_BASE_H
#define VPIC_GRID_BASE_H

#include "vpic/vpic.h"

// ======================================================================
// VpicGridBase

struct VpicGridBase : grid_t
{
  enum {
    anti_symmetric_fields = ::anti_symmetric_fields,
    pec_fields            = ::pec_fields,
    symmetric_fields      = ::symmetric_fields,
    pmc_fields            = ::pmc_fields,
    absorb_fields         = ::absorb_fields,
  } Fbc;

  enum {
    reflect_particles = ::reflect_particles,
    absorb_particles  = ::absorb_particles,
  } Pbc;

  static VpicGridBase* create()
  {
    return static_cast<VpicGridBase*>(::new_grid());
  }

  void setup(double dx_[3], double dt_, double cvac_, double eps0_)
  {
    dx = dx_[0]; dy = dx_[1]; dz = dx_[2];
    dt = dt_;
    cvac = cvac_;
    eps0 = eps0_;
  }

  void partition_periodic_box(const double xl[3], const double xh[3], const int gdims[3], const int np[3])
  {
    ::partition_periodic_box(this, xl[0], xl[1], xl[2], xh[0], xh[1], xh[2],
			     gdims[0], gdims[1], gdims[2], np[0], np[1], np[2]);
  }

  void set_fbc(int boundary, int fbc) { ::set_fbc(this, boundary, fbc); }
  void set_pbc(int boundary, int pbc) { ::set_pbc(this, boundary, pbc); }

  void mp_size_recv_buffer(int port, int size) { ::mp_size_recv_buffer(mp, port, size); }
  void mp_size_send_buffer(int port, int size) { ::mp_size_send_buffer(mp, port, size); }
  void* mp_send_buffer(int port) { return ::mp_send_buffer(mp, port); };
  void* mp_recv_buffer(int port) { return ::mp_recv_buffer(mp, port); };
  void mp_begin_recv(int port, int size, int src, int tag) { ::mp_begin_recv(mp, port, size, src, tag); }
  void mp_end_recv(int port) { ::mp_end_recv(mp, port); }
  void mp_begin_send(int port, int size, int dst, int tag) { ::mp_begin_send(mp, port, size, dst, tag); }
  void mp_end_send(int port) { ::mp_end_send(mp, port); }
  
  void* size_send_port(int i, int j, int k, int sz) { return ::size_send_port(i, j, k, sz, this); }
  void begin_send_port(int i, int j, int k, int sz) { ::begin_send_port(i, j, k, sz, this); }
  void end_send_port(int i, int j, int k) { ::end_send_port(i, j, k, this); }
  void begin_recv_port(int i, int j, int k, int sz) { ::begin_recv_port(i, j, k, sz, this); }
  void* end_recv_port(int i, int j, int k) { return ::end_recv_port(i, j, k, this); }
};

#endif
