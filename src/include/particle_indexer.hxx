
#pragma once

// ======================================================================
// Particle positions in cell / patch

// Particle positions are stored as patch-relative positions, however,
// it is required to know the exact cell a particle is in in a number of places:
// - for performance sorting particles
// - for keeping particles sorted in, e.g., the CUDA particle pusher
// - for finding the appropriate EM fields to interpolate
// - for correctly depositing currents
// - for collisions
//
// The goal here is to establish rules / invariants of the position of
// particles to where (in which patch) they are stored and how to
// recover the cell they are in.
//
// To complicate things, there are currently two aspects to this: Cell
// position and block position, where the former refers to the
// computational mesh that E, B and J live on, whereas a block refers
// to a fixed-size super-cell (e.g., 4x4x4 cells), motivated by
// performance considerations. It is currently not necessarily clear
// that the calculated block indices and cell indices satisfy the
// anticipated relation (bpos[d] = cpos[d] / bs[d]) because of potential
// finite precision arithmetic
//
// Rules / invariants:
//
// for all particles in a given patch,
// (1) the cell position cpos
//     calculated in a given dimension d will satisfy 0 <= cpos < ldims[d]
// (2) the patch relative position xi satisfies
//     0 <= xi <= xm[d] = ldims[d] * dx[d]
//     with the equality at the upper limit only being allowed at a right/top
//     non-periodic boundary
//
// These invariants will be temporarily violated after the particle push, but
// will be restored by the bnd exchange.
//
// Tricky issues to deal with:
// - (1) and (2) should be closely related, but finite precision
//   arithmetic can cause surprises.
//
//   E.g.: dx = 1 ldims = 100. Particle ends up at position -1e-7. It
//   gets sent to the left, where it's new position will be -1e-7 +
//   100., which is actually = 100. (in single precision), meaning that
//   particle violates (1) and (2) in its new patch.
//
// - Calculating cpos correctly when legally xi == ldims[d] * xi[d] at a right
// boundary
//
// TODO:
// - have cell index be the primary quantity computed, always, and find
//   block index from that
// - boundary exchange should be based on cell, not block index
//   (though if the two indices are always consistent, it becomes a non-issue)

// ======================================================================
// ParticleIndexer

template <class R>
struct ParticleIndexer
{
  using real_t = R;
  using Real3 = Vec3<real_t>;

  ParticleIndexer(const Grid_t& grid)
    : dxi_(grid.domain.dx_inv), ldims_(grid.ldims)
  {
    n_cells_ = ldims_[0] * ldims_[1] * ldims_[2];
  }

  int cellPosition(real_t x, int d) const { return fint(x * dxi_[d]); }

  Int3 cellPosition(const real_t* pos) const
  {
    Int3 idx;
    for (int d = 0; d < 3; d++) {
      idx[d] = cellPosition(pos[d], d);
    }
    return idx;
  }

  int cellIndex(const Int3& cpos) const
  {
    if (uint(cpos[0]) >= ldims_[0] || uint(cpos[1]) >= ldims_[1] ||
        uint(cpos[2]) >= ldims_[2]) {
      return -1;
    }

    return (cpos[2] * ldims_[1] + cpos[1]) * ldims_[0] + cpos[0];
  }

  int cellIndex(const real_t* pos) const
  {
    Int3 cpos = cellPosition(pos);
    return cellIndex(cpos);
  }

  bool isValidCellPosition(const Int3& cpos) const
  {
    for (int d = 0; d < 3; d++) {
      // optimization trick: if cpos[d]<0, it becomes larger than any positive
      // int after casting both to uint
      // TODO (cursed?) store ldims as uints, so this happens implicitly
      if (uint(cpos[d]) >= ldims_[d]) {
        return false;
      }
    }
    return true;
  }

  int validCellIndex(const real_t* pos) const
  {
    Int3 cpos = cellPosition(pos);
    for (int d = 0; d < 3; d++) {
      // optimization trick: if cpos[d]<0, it becomes larger than any positive
      // int after casting both to uint
      if (uint(cpos[d]) >= ldims_[d]) {
        printf("validCellIndex: cpos[%d] = %d ldims_[%d] = %d // pos[%d] = %g "
               "pos[%d]*dxi_[%d] = %g\n",
               d, cpos[d], d, ldims_[d], d, pos[d], d, d, pos[d] * dxi_[d]);
        assert(0);
      }
    }
    int cidx = cellIndex(cpos);
    assert(cidx >= 0);
    return cidx;
  }

  void checkInPatchMod(real_t* xi) const
  {
    for (int d = 0; d < 3; d++) {
      int pos = cellPosition(xi[d], d);
      if (pos < 0 || pos >= ldims_[d]) {
        printf("checkInPatchMod xi %g %g %g\n", xi[0], xi[1], xi[2]);
        printf("checkInPatchMod d %d xi %g pos %d // %d\n", d, xi[d], pos,
               ldims_[d]);
        if (pos < 0) {
          xi[d] = 0.f;
        } else {
          xi[d] *= (1. - 1e-6);
        }
        pos = cellPosition(xi[d], d);
      }
      assert(pos >= 0 && pos < ldims_[d]);
    }
  }

  const Int3& ldims() const { return ldims_; }

  // private:
  Real3 dxi_;
  Int3 ldims_;
  uint n_cells_;
};
