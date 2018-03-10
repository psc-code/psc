
#include "fields_item.hxx"

#include <string>

#define DEPOSIT_TO_GRID_1ST_CC(part, flds, m, val) do {			\
    Fields F(flds);							\
    real_t *xi = &part->xi; /* don't shift back in time */		\
    real_t u = xi[0] * dxi - .5;					\
    real_t v = xi[1] * dyi - .5;					\
    real_t w = xi[2] * dzi - .5;					\
    int jx = fint(u);							\
    int jy = fint(v);							\
    int jz = fint(w);							\
    real_t h1 = u - jx;							\
    real_t h2 = v - jy;							\
    real_t h3 = w - jz;							\
    									\
    real_t g0x = 1.f - h1;						\
    real_t g0y = 1.f - h2;						\
    real_t g0z = 1.f - h3;						\
    real_t g1x = h1;							\
    real_t g1y = h2;							\
    real_t g1z = h3;							\
    									\
    int jxd = 1, jyd = 1, jzd = 1;					\
    if (ppsc->domain.gdims[0] == 1) {					\
      jx = 0; g0x = 1.; g1x = 0.; jxd = 0;				\
    }									\
    if (ppsc->domain.gdims[1] == 1) {					\
      jy = 0; g0y = 1.; g1y = 0.; jyd = 0;				\
    }									\
    if (ppsc->domain.gdims[2] == 1) {					\
      jz = 0; g0z = 1.; g1z = 0.; jzd = 0;				\
    }									\
    									\
    assert(jx >= -1 && jx < grid.ldims[0]);				\
    assert(jy >= -1 && jy < grid.ldims[1]);				\
    assert(jz >= -1 && jz < grid.ldims[2]);				\
    									\
    real_t fnq = prts.prt_wni(*part) * fnqs;				\
    									\
    F(m, jx    ,jy    ,jz    ) += fnq*g0x*g0y*g0z * (val);		\
    F(m, jx+jxd,jy    ,jz    ) += fnq*g1x*g0y*g0z * (val);		\
    F(m, jx    ,jy+jyd,jz    ) += fnq*g0x*g1y*g0z * (val);		\
    F(m, jx+jxd,jy+jyd,jz    ) += fnq*g1x*g1y*g0z * (val);		\
    F(m, jx    ,jy    ,jz+jzd) += fnq*g0x*g0y*g1z * (val);		\
    F(m, jx+jxd,jy    ,jz+jzd) += fnq*g1x*g0y*g1z * (val);		\
    F(m, jx    ,jy+jyd,jz+jzd) += fnq*g0x*g1y*g1z * (val);		\
    F(m, jx+jxd,jy+jyd,jz+jzd) += fnq*g1x*g1y*g1z * (val);		\
  } while (0)

#define DEPOSIT_TO_GRID_1ST_NC(part, flds, m, val) do {			\
    Fields F(flds);							\
    real_t *xi = &part->xi; /* don't shift back in time */		\
    real_t u = xi[0] * dxi;						\
    real_t v = xi[1] * dyi;						\
    real_t w = xi[2] * dzi;						\
    int jx = fint(u);							\
    int jy = fint(v);							\
    int jz = fint(w);							\
    real_t h1 = u - jx;							\
    real_t h2 = v - jy;							\
    real_t h3 = w - jz;							\
    									\
    real_t g0x = 1.f - h1;						\
    real_t g0y = 1.f - h2;						\
    real_t g0z = 1.f - h3;						\
    real_t g1x = h1;							\
    real_t g1y = h2;							\
    real_t g1z = h3;							\
    									\
    int jxd = 1, jyd = 1, jzd = 1;					\
    if (ppsc->domain.gdims[0] == 1) {					\
      jx = 0; g0x = 1.; g1x = 0.; jxd = 0;				\
    }									\
    if (ppsc->domain.gdims[1] == 1) {					\
      jy = 0; g0y = 1.; g1y = 0.; jyd = 0;				\
    }									\
    if (ppsc->domain.gdims[2] == 1) {					\
      jz = 0; g0z = 1.; g1z = 0.; jzd = 0;				\
    }									\
    									\
    assert(jx >= -1 && jx < grid.ldims[0]);				\
    assert(jy >= -1 && jy < grid.ldims[1]);				\
    assert(jz >= -1 && jz < grid.ldims[2]);				\
    									\
    real_t fnq = prts.prt_wni(*part) * fnqs;				\
									\
    F(m, jx    ,jy    ,jz    ) += fnq*g0x*g0y*g0z * (val);		\
    F(m, jx+jxd,jy    ,jz    ) += fnq*g1x*g0y*g0z * (val);		\
    F(m, jx    ,jy+jyd,jz    ) += fnq*g0x*g1y*g0z * (val);		\
    F(m, jx+jxd,jy+jyd,jz    ) += fnq*g1x*g1y*g0z * (val);		\
    F(m, jx    ,jy    ,jz+jzd) += fnq*g0x*g0y*g1z * (val);		\
    F(m, jx+jxd,jy    ,jz+jzd) += fnq*g1x*g0y*g1z * (val);		\
    F(m, jx    ,jy+jyd,jz+jzd) += fnq*g0x*g1y*g1z * (val);		\
    F(m, jx+jxd,jy+jyd,jz+jzd) += fnq*g1x*g1y*g1z * (val);		\
  } while (0)

#define DEPOSIT_TO_GRID_2ND_NC(part, flds, m, val) do {			\
    Fields F(flds);							\
    real_t *xi = &part->xi; /* don't shift back in time */		\
    real_t u = xi[0] * dxi;						\
    real_t v = xi[1] * dyi;						\
    real_t w = xi[2] * dzi;						\
    int jx = nint(u);							\
    int jy = nint(v);							\
    int jz = nint(w);							\
    real_t h1 = jx - u;							\
    real_t h2 = jy - v;							\
    real_t h3 = jz - w;							\
    									\
    real_t gmx = .5f*(.5f+h1)*(.5f+h1);					\
    real_t gmy = .5f*(.5f+h2)*(.5f+h2);					\
    real_t gmz = .5f*(.5f+h3)*(.5f+h3);					\
    real_t g0x = .75f-h1*h1;						\
    real_t g0y = .75f-h2*h2;						\
    real_t g0z = .75f-h3*h3;						\
    real_t g1x = .5f*(.5f-h1)*(.5f-h1);					\
    real_t g1y = .5f*(.5f-h2)*(.5f-h2);					\
    real_t g1z = .5f*(.5f-h3)*(.5f-h3);					\
    									\
    int jxd = 1, jyd = 1, jzd = 1;					\
    if (ppsc->domain.gdims[0] == 1) {					\
      jx = 0; g0x = 1.; g1x = 0.; gmx = 0.; jxd = 0;			\
    }									\
    if (ppsc->domain.gdims[1] == 1) {					\
      jy = 0; g0y = 1.; g1y = 0.; gmy = 0.; jyd = 0;			\
    }									\
    if (ppsc->domain.gdims[2] == 1) {					\
      jz = 0; g0z = 1.; g1z = 0.; gmz = 0.; jzd = 0;			\
    }									\
    									\
    assert(jx >= 0 && jx <= grid.ldims[0]);				\
    assert(jy >= 0 && jy <= grid.ldims[1]);				\
    assert(jz >= 0 && jz <= grid.ldims[2]);				\
    									\
    real_t fnq = prts.prt_wni(*part) * fnqs;				\
					     				\
    F(m, jx-jxd,jy-jyd,jz-jzd) += fnq*gmx*gmy*gmz * (val);		\
    F(m, jx    ,jy-jyd,jz-jzd) += fnq*g0x*gmy*gmz * (val);		\
    F(m, jx+jxd,jy-jyd,jz-jzd) += fnq*g1x*gmy*gmz * (val);		\
    F(m, jx-jxd,jy    ,jz-jzd) += fnq*gmx*g0y*gmz * (val);		\
    F(m, jx    ,jy    ,jz-jzd) += fnq*g0x*g0y*gmz * (val);		\
    F(m, jx+jxd,jy    ,jz-jzd) += fnq*g1x*g0y*gmz * (val);		\
    F(m, jx-jxd,jy+jyd,jz-jzd) += fnq*gmx*g1y*gmz * (val);		\
    F(m, jx    ,jy+jyd,jz-jzd) += fnq*g0x*g1y*gmz * (val);		\
    F(m, jx+jxd,jy+jyd,jz-jzd) += fnq*g1x*g1y*gmz * (val);		\
    F(m, jx-jxd,jy-jyd,jz    ) += fnq*gmx*gmy*g0z * (val);		\
    F(m, jx    ,jy-jyd,jz    ) += fnq*g0x*gmy*g0z * (val);		\
    F(m, jx+jxd,jy-jyd,jz    ) += fnq*g1x*gmy*g0z * (val);		\
    F(m, jx-jxd,jy    ,jz    ) += fnq*gmx*g0y*g0z * (val);		\
    F(m, jx    ,jy    ,jz    ) += fnq*g0x*g0y*g0z * (val);		\
    F(m, jx+jxd,jy    ,jz    ) += fnq*g1x*g0y*g0z * (val);		\
    F(m, jx-jxd,jy+jyd,jz    ) += fnq*gmx*g1y*g0z * (val);		\
    F(m, jx    ,jy+jyd,jz    ) += fnq*g0x*g1y*g0z * (val);		\
    F(m, jx+jxd,jy+jyd,jz    ) += fnq*g1x*g1y*g0z * (val);		\
    F(m, jx-jxd,jy-jyd,jz+jzd) += fnq*gmx*gmy*g1z * (val);		\
    F(m, jx    ,jy-jyd,jz+jzd) += fnq*g0x*gmy*g1z * (val);		\
    F(m, jx+jxd,jy-jyd,jz+jzd) += fnq*g1x*gmy*g1z * (val);		\
    F(m, jx-jxd,jy    ,jz+jzd) += fnq*gmx*g0y*g1z * (val);		\
    F(m, jx    ,jy    ,jz+jzd) += fnq*g0x*g0y*g1z * (val);		\
    F(m, jx+jxd,jy    ,jz+jzd) += fnq*g1x*g0y*g1z * (val);		\
    F(m, jx-jxd,jy+jyd,jz+jzd) += fnq*gmx*g1y*g1z * (val);		\
    F(m, jx    ,jy+jyd,jz+jzd) += fnq*g0x*g1y*g1z * (val);		\
    F(m, jx+jxd,jy+jyd,jz+jzd) += fnq*g1x*g1y*g1z * (val);		\
  } while (0)

// FIXME, this function exists about 100x all over the place, should
// be consolidated

static inline void
particle_calc_vxi(particle_t *part, real_t vxi[3])
{
  real_t root =
    1.f / std::sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

template<typename Moment_t>
struct ItemMomentWrap
{
  static void run(mparticles_t mprts, mfields_t mflds_res)
  {
    for (int p = 0; p < mprts->n_patches(); p++) {
      mflds_res[p].zero();
      Moment_t::run(mflds_res[p], mprts[p]);
      add_ghosts_boundary(mflds_res[p], p, 0, mflds_res->n_comps());
    }
  }

  // ----------------------------------------------------------------------
  // boundary stuff FIXME, should go elsewhere...

  static void add_ghosts_reflecting_lo(fields_t flds, int p, int d, int mb, int me)
  {
    Fields F(flds);
    const int *ldims = ppsc->grid().ldims;

    int bx = ldims[0] == 1 ? 0 : 1;
    if (d == 1) {
      for (int iz = -1; iz < ldims[2] + 1; iz++) {
	for (int ix = -bx; ix < ldims[0] + bx; ix++) {
	  int iy = 0; {
	    for (int m = mb; m < me; m++) {
	      F(m, ix,iy,iz) += F(m, ix,iy-1,iz);
	    }
	  }
	}
      }
    } else if (d == 2) {
      for (int iy = 0*-1; iy < ldims[1] + 0*1; iy++) {
	for (int ix = -bx; ix < ldims[0] + bx; ix++) {
	  int iz = 0; {
	    for (int m = mb; m < me; m++) {
	      F(m, ix,iy,iz) += F(m, ix,iy,iz-1);
	    }
	  }
	}
      }
    } else {
      assert(0);
    }
  }

  static void add_ghosts_reflecting_hi(fields_t flds, int p, int d, int mb, int me)
  {
    Fields F(flds);
    const int *ldims = ppsc->grid().ldims;

    int bx = ldims[0] == 1 ? 0 : 1;
    if (d == 1) {
      for (int iz = -1; iz < ldims[2] + 1; iz++) {
	for (int ix = -bx; ix < ldims[0] + bx; ix++) {
	  int iy = ldims[1] - 1; {
	    for (int m = mb; m < me; m++) {
	      F(m, ix,iy,iz) += F(m, ix,iy+1,iz);
	    }
	  }
	}
      }
    } else if (d == 2) {
      for (int iy = 0*-1; iy < ldims[1] + 0*1; iy++) {
	for (int ix = -bx; ix < ldims[0] + bx; ix++) {
	  int iz = ldims[2] - 1; {
	    for (int m = mb; m < me; m++) {
	      F(m, ix,iy,iz) += F(m, ix,iy,iz+1);
	    }
	  }
	}
      }
    } else {
      assert(0);
    }
  }

  static void add_ghosts_boundary(fields_t res, int p, int mb, int me)
  {
    // lo
    for (int d = 0; d < 3; d++) {
      if (psc_at_boundary_lo(ppsc, p, d)) {
	if (ppsc->domain.bnd_part_lo[d] == BND_PART_REFLECTING ||
	    ppsc->domain.bnd_part_lo[d] == BND_PART_OPEN) {
	  add_ghosts_reflecting_lo(res, p, d, mb, me);
	}
      }
    }
    // hi
    for (int d = 0; d < 3; d++) {
      if (psc_at_boundary_hi(ppsc, p, d)) {
	if (ppsc->domain.bnd_part_hi[d] == BND_PART_REFLECTING ||
	    ppsc->domain.bnd_part_hi[d] == BND_PART_OPEN) {
	  add_ghosts_reflecting_hi(res, p, d, mb, me);
	}
      }
    }
  }
};

template<typename Moment_t, typename mparticles_t>
struct ItemMoment : FieldsItemCRTP<ItemMoment<Moment_t, mparticles_t>>
{
  using Base = FieldsItemCRTP<ItemMoment<Moment_t, mparticles_t>>;
  using Base::Base;
  
  static const char* name()
  {
    return strdup((std::string(Moment_t::name) + "_" +
		   mparticles_traits<mparticles_t>::name).c_str());
  }
  constexpr static int n_comps = Moment_t::n_comps; 
  constexpr static fld_names_t fld_names() { return Moment_t::fld_names(); }
  constexpr static int flags = Moment_t::flags;

  void run(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) override
  {
    mparticles_t mprts = mprts_base.get_as<mparticles_t>();
    mfields_t mres = this->mres_base_->template get_as<mfields_t>(0, 0);

    ItemMomentWrap<Moment_t>::run(mprts, mres);
    
    mres.put_as(this->mres_base_, 0, this->mres_base_->nr_fields);
    mprts.put_as(mprts_base, MP_DONT_COPY);
  }
};
  
template<typename Item_t, typename mparticles_t>
using FieldsItemMomentOps = FieldsItemOps<ItemMoment<Item_t, mparticles_t>>;

