
#ifndef PSC_FIELD_ARRAY_REMOTE_OPS_H
#define PSC_FIELD_ARRAY_REMOTE_OPS_H

// ======================================================================
// foreach

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

// ======================================================================
// Comm

template<class F3D>
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

  float* get_send_buf(int dir, int side, int sz) const
  {
    const int* ijk = to_ijk[dir][side];
    return static_cast<float*>(::size_send_port(ijk[0], ijk[1], ijk[2], sz * sizeof(float), g_));
  }

  void begin_send_port(int dir, int side, int sz) const
  {
    const int* ijk = to_ijk[dir][side];
    ::begin_send_port(ijk[0], ijk[1], ijk[2], sz * sizeof(float), g_);
  }

  void begin_recv_port(int dir, int side, int sz) const
  {
    const int* ijk = to_ijk[dir][side];
    ::begin_recv_port(ijk[0], ijk[1], ijk[2], sz * sizeof(float), g_);
  }

  float* end_recv_port(int dir, int side) const
  {
    const int* ijk = to_ijk[dir][side];
    return static_cast<float*>(::end_recv_port(ijk[0], ijk[1], ijk[2], g_));
  }

  void end_send_port(int dir, int side) const
  {
    const int* ijk = to_ijk[dir][side];
    ::end_send_port(ijk[0], ijk[1], ijk[2], g_);
  }

  // ----------------------------------------------------------------------

  virtual void begin_send(int side, int dir, float* p, F3D& F) const = 0;
  virtual void end_recv(int side, int dir, float* p, F3D& F) const = 0;
  
  void begin(F3D& F) const
  {
    for (int side = 0; side <= 1; side++) {
      for (int dir = 0; dir < 3; dir++) {
	begin_recv_port(dir, side, buf_size_[dir]);
      }
    }

    for (int side = 0; side <= 1; side++) {
      for (int dir = 0; dir < 3; dir++) {
	float *p = get_send_buf(dir, side, buf_size_[dir]);
	if (p) {
	  begin_send(dir, side, p, F);
	  begin_send_port(dir, side, buf_size_[dir]);
	}
      }
    }
  }
  
  void end(F3D& F) const
  {
    for (int side = 0; side <= 1; side++) {
      for (int dir = 0; dir < 3; dir++) {
	float* p = end_recv_port(dir, side);
	if (p) {
	  end_recv(dir, side, p, F);
	}
      }
    }

    for (int side = 0; side <= 1; side++) {
      for (int dir = 0; dir < 3; dir++) {
	end_send_port(dir, side);
      }
    }
  }
  
protected:
  int nx_[3];
  int buf_size_[3];
  const grid_t *g_;
};

// ======================================================================
// PscFieldArrayRemoteOps

template<class FA>
struct PscFieldArrayRemoteOps {
  typedef FA FieldArray;

  // ----------------------------------------------------------------------
  // CommEC

  template<class F3D>
  struct CommEC : Comm<F3D>
  {
    typedef Comm<F3D> Base;
    using Base::begin;
    using Base::end;

    using Base::nx_;
    using Base::g_;
    using Base::buf_size_;

    CommEC(grid_t *g) : Base(g)
    {
      for (int X = 0; X < 3; X++) {
	int Y = (X + 1) % 3, Z = (X + 2) % 3;
	buf_size_[X] = nx_[Y] * (nx_[Z] + 1) + nx_[Z] * (nx_[Y] + 1);
      }
    }

    void begin_send(int X, int side, float* p, F3D& F) const
    {
      int Y = (X + 1) % 3, Z = (X + 2) % 3;
      int face = side ? nx_[X] : 1;
      foreach_edge(g_, Z, Y, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).cbx)[Y]; });
      foreach_edge(g_, Y, Z, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).cbx)[Z]; });
    }
    
    void end_recv(int X, int side, float* p, F3D& F) const
    {
      int Y = (X + 1) % 3, Z = (X + 2) % 3;
      int face = side ? 0 : nx_[X] + 1;
      foreach_edge(g_, Z, Y, face, [&](int x, int y, int z) { (&F(x,y,z).cbx)[Y] = *p++; });
      foreach_edge(g_, Y, Z, face, [&](int x, int y, int z) { (&F(x,y,z).cbx)[Z] = *p++; });
    }
  };

  void begin_remote_ghost_tang_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    CommEC<Field3D<FieldArray>> comm(fa.g);

    comm.begin(F);
  }

  void end_remote_ghost_tang_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    CommEC<Field3D<FieldArray>> comm(fa.g);
    
    comm.end(F);
  }

  // ----------------------------------------------------------------------
  // CommNC

  template<class F3D>
  struct CommNC : Comm<F3D>
  {
    typedef Comm<F3D> Base;
    using Base::begin;
    using Base::end;

    using Base::nx_;
    using Base::g_;
    using Base::buf_size_;

    CommNC(grid_t *g) : Base(g)
    {
      for (int X = 0; X < 3; X++) {
	int Y = (X + 1) % 3, Z = (X + 2) % 3;
	buf_size_[X] = (nx_[Y] + 1) * (nx_[Z] + 1);
      }
    }
    
    void begin_send(int X, int side, float *p, F3D& F) const
    {
      int face = side ? nx_[X] : 1;
      foreach_node(g_, X, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).ex)[X]; });
    }
    
    void end_recv(int X, int side, float *p, F3D& F) const
    {
      int face = side ? 0 : nx_[X] + 1;
      foreach_node(g_, X, face, [&](int x, int y, int z) { (&F(x,y,z).ex)[X] = *p++; });
    }
  };
  
  void begin_remote_ghost_norm_e(FieldArray &fa)
  {
    Field3D<FieldArray> F(fa);
    CommNC<Field3D<FieldArray>> comm(fa.g);

    comm.begin(F);
  }

  void end_remote_ghost_norm_e(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    CommNC<Field3D<FieldArray>> comm(fa.g);

    comm.end(F);
  }

  // ----------------------------------------------------------------------
  // CommCC

  template<class F3D>
  struct CommCC : Comm<F3D>
  {
    typedef Comm<F3D> Base;
    using Base::begin;
    using Base::end;
  
    using Base::nx_;
    using Base::g_;
    using Base::buf_size_;

    CommCC(grid_t *g) : Base(g)
    {
      for (int X = 0; X < 3; X++) {
	int Y = (X + 1) % 3, Z = (X + 2) % 3;
	buf_size_[X] = nx_[Y] * nx_[Z];
      }
    }
    
    void begin_send(int X, int side, float *p, F3D& F) const
    {
      int face = side ? nx_[X] : 1;
      foreach_center(g_, X, face, [&](int x, int y, int z) { *p++ = F(x,y,z).div_b_err; });
    }
    
    void end_recv(int X, int side, float *p, F3D& F) const
    {
      int face = side ? 0 : nx_[X] + 1;
      foreach_center(g_, X, face, [&](int x, int y, int z) { F(x,y,z).div_b_err = *p++; });
    }
  };
     
  void begin_remote_ghost_div_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    CommCC<Field3D<FieldArray>> comm(fa.g);

    comm.begin(F);
  }

  void end_remote_ghost_div_b(FieldArray& fa)
  {
    Field3D<FieldArray> F(fa);
    CommCC<Field3D<FieldArray>> comm(fa.g);

    comm.end(F);
  }

};


#endif

