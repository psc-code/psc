
#ifndef PSC_FIELD_ARRAY_REMOTE_OPS_H
#define PSC_FIELD_ARRAY_REMOTE_OPS_H

#include "GridLoop.h"

// ======================================================================
// PscFieldArrayRemoteOps

template<typename MfieldsState>
struct PscFieldArrayRemoteOps
{
  using Grid = typename MfieldsState::Grid;
  using F3D = Field3D<typename MfieldsState::Patch>;

  // ----------------------------------------------------------------------
  // CommEC

  template<class G, class F3D>
  struct CommEC : Comm<G, F3D>
  {
    typedef Comm<G, F3D> Base;
    typedef G Grid;
    using Base::begin;
    using Base::end;

    using Base::nx_;
    using Base::g_;
    using Base::buf_size_;

    CommEC(Grid *g) : Base(g)
    {
      for (int X = 0; X < 3; X++) {
	int Y = (X + 1) % 3, Z = (X + 2) % 3;
	buf_size_[X] = nx_[Y] * (nx_[Z] + 1) + nx_[Z] * (nx_[Y] + 1);
      }
    }

    void begin_send(int X, int side, float* p, F3D& F)
    {
      int Y = (X + 1) % 3, Z = (X + 2) % 3;
      int face = side ? nx_[X] : 1;
      foreach_edge(g_, Z, Y, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).cbx)[Y]; });
      foreach_edge(g_, Y, Z, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).cbx)[Z]; });
    }
    
    void end_recv(int X, int side, float* p, F3D& F)
    {
      int Y = (X + 1) % 3, Z = (X + 2) % 3;
      int face = side ? 0 : nx_[X] + 1;
      foreach_edge(g_, Z, Y, face, [&](int x, int y, int z) { (&F(x,y,z).cbx)[Y] = *p++; });
      foreach_edge(g_, Y, Z, face, [&](int x, int y, int z) { (&F(x,y,z).cbx)[Z] = *p++; });
    }
  };

  static void begin_remote_ghost_tang_b(MfieldsState &mflds)
  {
    auto& fa = mflds.getPatch(0);
    F3D F(fa);
    CommEC<Grid, F3D> comm(mflds.vgrid());

    comm.begin(F);
  }

  static void end_remote_ghost_tang_b(MfieldsState& mflds)
  {
    auto& fa = mflds.getPatch(0);
    F3D F(fa);
    CommEC<Grid, F3D> comm(mflds.vgrid());
    
    comm.end(F);
  }

  // ----------------------------------------------------------------------
  // CommNC

  template<class G, class F3D>
  struct CommNC : Comm<G, F3D>
  {
    typedef Comm<G, F3D> Base;
    typedef G Grid;
    using Base::begin;
    using Base::end;

    using Base::nx_;
    using Base::g_;
    using Base::buf_size_;

    CommNC(Grid *g) : Base(g)
    {
      for (int X = 0; X < 3; X++) {
	int Y = (X + 1) % 3, Z = (X + 2) % 3;
	buf_size_[X] = (nx_[Y] + 1) * (nx_[Z] + 1);
      }
    }
    
    void begin_send(int X, int side, float *p, F3D& F)
    {
      int face = side ? nx_[X] : 1;
      foreach_node(g_, X, face, [&](int x, int y, int z) { *p++ = (&F(x,y,z).ex)[X]; });
    }
    
    void end_recv(int X, int side, float *p, F3D& F)
    {
      int face = side ? 0 : nx_[X] + 1;
      foreach_node(g_, X, face, [&](int x, int y, int z) { (&F(x,y,z).ex)[X] = *p++; });
    }
  };
  
  static void begin_remote_ghost_norm_e(MfieldsState &mflds)
  {
    auto& fa = mflds.getPatch(0);
    F3D F(fa);
    CommNC<Grid, F3D> comm(mflds.vgrid());

    comm.begin(F);
  }

  static void end_remote_ghost_norm_e(MfieldsState& mflds)
  {
    auto& fa = mflds.getPatch(0);
    F3D F(fa);
    CommNC<Grid, F3D> comm(mflds.vgrid());

    comm.end(F);
  }

  // ----------------------------------------------------------------------
  // CommCC

  template<class G, class F3D>
  struct CommCC : Comm<G, F3D>
  {
    typedef Comm<G, F3D> Base;
    typedef G Grid;
    using Base::begin;
    using Base::end;
  
    using Base::nx_;
    using Base::g_;
    using Base::buf_size_;

    CommCC(Grid *g) : Base(g)
    {
      for (int X = 0; X < 3; X++) {
	int Y = (X + 1) % 3, Z = (X + 2) % 3;
	buf_size_[X] = nx_[Y] * nx_[Z];
      }
    }
    
    void begin_send(int X, int side, float *p, F3D& F)
    {
      int face = side ? nx_[X] : 1;
      foreach_face(g_, X, face, [&](int x, int y, int z) { *p++ = F(x,y,z).div_b_err; });
    }
    
    void end_recv(int X, int side, float *p, F3D& F)
    {
      int face = side ? 0 : nx_[X] + 1;
      foreach_face(g_, X, face, [&](int x, int y, int z) { F(x,y,z).div_b_err = *p++; });
    }
  };
     
  static void begin_remote_ghost_div_b(MfieldsState& mflds)
  {
    auto& fa = mflds.getPatch(0);
    F3D F(fa);
    CommCC<Grid, F3D> comm(mflds.vgrid());

    comm.begin(F);
  }

  static void end_remote_ghost_div_b(MfieldsState& mflds)
  {
    auto& fa = mflds.getPatch(0);
    F3D F(fa);
    CommCC<Grid, F3D> comm(mflds.vgrid());

    comm.end(F);
  }

};


#endif

