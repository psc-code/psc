
#ifndef PSC_FIELD_ARRAY_REMOTE_OPS_H
#define PSC_FIELD_ARRAY_REMOTE_OPS_H

// ======================================================================
// PscFieldArrayRemoteOps

template<class FA>
struct PscFieldArrayRemoteOps {
  typedef FA FieldArray;

  template<class F>
  static void foreach(F f, int ib, int ie, int jb, int je, int kb, int ke)
  {
    for (int k = kb; k <= ke; k++) {
      for (int j = jb; j <= je; j++) {
	for (int i = ib; i <= ie; i++) {
	  f(i,j,k);
	}
      }
    }
  };

  template<class F>
  static void foreach_edge(const grid_t *g, int Y, int Z, int face, F f)
  {
    if        (Y == 1 && Z == 2) { foreach(f, face, face, 1, g->ny  , 1, g->nz+1);
    } else if (Y == 2 && Z == 1) { foreach(f, face, face, 1, g->ny+1, 1, g->nz  );
    } else if (Y == 2 && Z == 0) { foreach(f, 1, g->nx+1, face, face, 1, g->nz  );
    } else if (Y == 0 && Z == 2) { foreach(f, 1, g->nx  , face, face, 1, g->nz+1);
    } else if (Y == 0 && Z == 1) { foreach(f, 1, g->nx  , 1, g->ny+1, face, face);
    } else if (Y == 1 && Z == 0) { foreach(f, 1, g->nx+1, 1, g->ny  , face, face);
    } else {
      assert(0);
    }
  }
  
  template<class F>
  static void foreach_node(const grid_t *g, int X, int face, F f)
  {
    if        (X == 0) { foreach(f, face, face, 1, g->ny+1, 1, g->nz+1);
    } else if (X == 1) { foreach(f, 1, g->nx+1, face, face, 1, g->nz+1);
    } else if (X == 2) { foreach(f, 1, g->nx+1, 1, g->ny+1, face, face);
    } else {
      assert(0);
    }
  }

  template<class F>
  static void foreach_center(const grid_t *g, int X, int face, F f)
  {
    if        (X == 0) { foreach(f, face, face, 1, g->ny  , 1, g->nz  );
    } else if (X == 1) { foreach(f, 1, g->nx  , face, face, 1, g->nz  );
    } else if (X == 2) { foreach(f, 1, g->nx  , 1, g->ny  , face, face);
    } else {
      assert(0);
    }
  }

  // ----------------------------------------------------------------------
  // Comm

  struct Comm
  {
    Comm(const grid_t* g) : g_(g)
    {
      nx_[0] = g_->nx;
      nx_[1] = g_->ny;
      nx_[2] = g_->nz;
    }

    const int to_ijk[3][2][3] = {
      { { -1, 0, 0 },
	{  1, 0, 0 }, },
      { {  0,-1, 0 },
	{  0, 1, 0 }, },
      { {  0, 0,-1 },
	{  0, 0, 1 }, },
    };
    // wrap what probably should be in Grid::

    float* get_send_buf(int i, int j, int k, int sz) const
    {
      return static_cast<float*>(::size_send_port(i, j, k, sz * sizeof(float), g_));
    }

    void begin_send_port(int i, int j, int k, int sz) const
    {
      ::begin_send_port(i, j, k, sz * sizeof(float), g_);
    }

    void begin_recv_port(int i, int j, int k, int sz) const
    {
      ::begin_recv_port(i, j, k, sz * sizeof(float), g_);
    }

    float* end_recv_port(int i, int j, int k) const
    {
      return static_cast<float*>(::end_recv_port(i, j, k, g_));
    }

    void end_send(int dir, int side) const
    {
      const int* ijk = to_ijk[dir][side];
      ::end_send_port(ijk[0], ijk[1], ijk[2], g_);
    }

    // ghost communication
    
    void nc_begin_recv(int dir, int side, int X, int Y, int Z) const
    {
      const int* ijk = to_ijk[dir][side];
      int sz = (nx_[Y] + 1) * (nx_[Z] + 1);
      begin_recv_port(ijk[0], ijk[1], ijk[2], sz);
    }

    void ec_begin_recv(int dir, int side, int X, int Y, int Z) const
    {
      const int* ijk = to_ijk[dir][side];
      int sz = nx_[Y] * (nx_[Z]+1) + nx_[Z] * (nx_[Y]+1);
      begin_recv_port(ijk[0], ijk[1], ijk[2], sz);
    }

    void cc_begin_recv(int dir, int side, int X, int Y, int Z) const
    {
      const int* ijk = to_ijk[dir][side];
      int sz = nx_[Y] * nx_[Z];
      begin_recv_port(ijk[0], ijk[1], ijk[2], sz);
    }

    template<class F3D>
    void nc_begin_send(int dir, int side, int X, int Y, int Z, F3D& F) const
    {
      const int* ijk = to_ijk[dir][side];
      int sz = (nx_[Y] + 1) * (nx_[Z] + 1);
      float *p = get_send_buf(ijk[0], ijk[1], ijk[2], sz);
      if (p) {
	int face = side ? nx_[X] : 1;
	foreach_node(g_, X, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).ex)[X]; });
	begin_send_port(ijk[0], ijk[1], ijk[2], sz);
      }
    }
    
    template<class F3D>
    void ec_begin_send(int dir, int side, int X, int Y, int Z, F3D& F) const
    {
      const int* ijk = to_ijk[dir][side];
      int sz = nx_[Y] * (nx_[Z]+1) + nx_[Z] * (nx_[Y]+1);
      float *p = get_send_buf(ijk[0], ijk[1], ijk[2], sz);
      if (p) {
	int face = side ? nx_[X] : 1;
	foreach_edge(g_, Z, Y, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).cbx)[Y]; });
	foreach_edge(g_, Y, Z, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).cbx)[Z]; });
	begin_send_port(ijk[0], ijk[1], ijk[2], sz);
      }
    }
    
    template<class F3D>
    void cc_begin_send(int dir, int side, int X, int Y, int Z, F3D& F) const
    {
      const int* ijk = to_ijk[dir][side];
      int sz = nx_[Y] * nx_[Z];
      float *p = get_send_buf(ijk[0], ijk[1], ijk[2], sz);
      if (p) {
	int face = side ? nx_[X] : 1;
	foreach_center(g_, X, face, [&](int x, int y, int z) { *p++ = F(x,y,z).div_b_err; });
	begin_send_port(ijk[0], ijk[1], ijk[2], sz);
      }
    }
    
    template<class F3D>
    void nc_end_recv(int dir, int side, int X, int Y, int Z, F3D& F) const
    {
      const int* ijk = to_ijk[dir][side];
      float* p = end_recv_port(ijk[0], ijk[1], ijk[2]);
      if (p) {
	int face = side ? 0 : nx_[X] + 1;
	foreach_node(g_, X, face, [&](int x, int y, int z) { (&F(x,y,z).ex)[X] = *p++; });
      }
    }

    template<class F3D>
    void ec_end_recv(int i, int j, int k, int X, int Y, int Z, F3D& F) const
    {
      float* p = end_recv_port(i,j,k);
      if (p) {
	int face = (i+j+k) < 0 ? nx_[X] + 1 : 0;
	foreach_edge(g_, Z, Y, face, [&](int x, int y, int z) { (&F(x,y,z).cbx)[Y] = *p++; });
	foreach_edge(g_, Y, Z, face, [&](int x, int y, int z) { (&F(x,y,z).cbx)[Z] = *p++; });
      }
    }

    template<class F3D>
    void cc_end_recv(int dir, int side, int X, int Y, int Z, F3D& F) const
    {
      const int* ijk = to_ijk[dir][side];
      float* p = end_recv_port(ijk[0], ijk[1], ijk[2]);
      if (p) {
	int face = side ? 0 : nx_[X] + 1;
	foreach_center(g_, X, face, [&](int x, int y, int z) { F(x,y,z).div_b_err = *p++; });
      }
    }

    //private:
    int nx_[3];
    const grid_t *g_;
  };

  void begin_remote_ghost_tang_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    Comm comm(fa.g);

    comm.ec_begin_recv(0, 0, 0,1,2);
    comm.ec_begin_recv(1, 0, 1,2,0);
    comm.ec_begin_recv(2, 0, 2,0,1);
    comm.ec_begin_recv(0, 1, 0,1,2);
    comm.ec_begin_recv(1, 1, 1,2,0);
    comm.ec_begin_recv(2, 1, 2,0,1);

    comm.ec_begin_send(0, 0, 0,1,2, F);
    comm.ec_begin_send(1, 0, 1,2,0, F);
    comm.ec_begin_send(2, 0, 2,0,1, F);
    comm.ec_begin_send(0, 1, 0,1,2, F);
    comm.ec_begin_send(1, 1, 1,2,0, F);
    comm.ec_begin_send(2, 1, 2,0,1, F);
  }

  void end_remote_ghost_tang_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    Comm comm(fa.g);
    
    comm.ec_end_recv(-1, 0, 0, 0,1,2, F);
    comm.ec_end_recv( 0,-1, 0, 1,2,0, F);
    comm.ec_end_recv( 0, 0,-1, 2,0,1, F);
    comm.ec_end_recv( 1, 0, 0, 0,1,2, F);
    comm.ec_end_recv( 0, 1, 0, 1,2,0, F);
    comm.ec_end_recv( 0, 0, 1, 2,0,1, F);

    comm.end_send(0, 0);
    comm.end_send(1, 0);
    comm.end_send(2, 0);
    comm.end_send(0, 1);
    comm.end_send(1, 1);
    comm.end_send(2, 1);
  }

  void begin_remote_ghost_norm_e(FieldArray &fa)
  {
    Field3D<FieldArray> F(fa);
    Comm comm(fa.g);

    comm.nc_begin_recv(0, 0, 0,1,2);
    comm.nc_begin_recv(1, 0, 1,2,0);
    comm.nc_begin_recv(2, 0, 2,0,1);
    comm.nc_begin_recv(0, 1, 0,1,2);
    comm.nc_begin_recv(1, 1, 1,2,0);
    comm.nc_begin_recv(2, 1, 2,0,1);

    comm.nc_begin_send(0, 0, 0,1,2, F);
    comm.nc_begin_send(1, 0, 1,2,0, F);
    comm.nc_begin_send(2, 0, 2,0,1, F);
    comm.nc_begin_send(0, 1, 0,1,2, F);
    comm.nc_begin_send(1, 1, 1,2,0, F);
    comm.nc_begin_send(2, 1, 2,0,1, F);
  }

  void end_remote_ghost_norm_e(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    Comm comm(fa.g);

    comm.nc_end_recv(0, 0, 0,1,2, F);
    comm.nc_end_recv(1, 0, 1,2,0, F);
    comm.nc_end_recv(2, 0, 2,0,1, F);
    comm.nc_end_recv(0, 1, 0,1,2, F);
    comm.nc_end_recv(1, 1, 1,2,0, F);
    comm.nc_end_recv(2, 1, 2,0,1, F);
    
    comm.end_send(0, 0);
    comm.end_send(1, 0);
    comm.end_send(2, 0);
    comm.end_send(0, 1);
    comm.end_send(1, 1);
    comm.end_send(2, 1);
  }

  void begin_remote_ghost_div_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    Comm comm(fa.g);

    comm.cc_begin_recv(0, 0, 0,1,2);
    comm.cc_begin_recv(1, 0, 1,2,0);
    comm.cc_begin_recv(2, 0, 2,0,1);
    comm.cc_begin_recv(0, 1, 0,1,2);
    comm.cc_begin_recv(1, 1, 1,2,0);
    comm.cc_begin_recv(2, 1, 2,0,1);

    comm.cc_begin_send(0, 0, 0,1,2, F);
    comm.cc_begin_send(1, 0, 1,2,0, F);
    comm.cc_begin_send(2, 0, 2,0,1, F);
    comm.cc_begin_send(0, 1, 0,1,2, F);
    comm.cc_begin_send(1, 1, 1,2,0, F);
    comm.cc_begin_send(2, 1, 2,0,1, F);
  }

  void end_remote_ghost_div_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    Comm comm(fa.g);

    comm.cc_end_recv(0, 0, 0,1,2, F);
    comm.cc_end_recv(1, 0, 1,2,0, F);
    comm.cc_end_recv(2, 0, 2,0,1, F);
    comm.cc_end_recv(0, 1, 0,1,2, F);
    comm.cc_end_recv(1, 1, 1,2,0, F);
    comm.cc_end_recv(2, 1, 2,0,1, F);
    
    comm.end_send(0, 0);
    comm.end_send(1, 0);
    comm.end_send(2, 0);
    comm.end_send(0, 1);
    comm.end_send(1, 1);
    comm.end_send(2, 1);
  }

};


#endif

