
#include "psc_marder_private.h"
#include "psc_bnd.h"
#include "psc_output_fields_item.h"
#include "fields.hxx"

#include <mrc_io.h>

template<typename MP, typename MF>
struct Marder_
{
  using Mparticles_t = MP;
  using Mfields = MF;
  using fields_t = typename Mfields::fields_t;
  using real_t = typename Mfields::real_t;
  using Fields = Fields3d<fields_t>;

  // FIXME: checkpointing won't properly restore state
  // FIXME: if the subclass creates objects, it'd be cleaner to have them
  // be part of the subclass

  // ----------------------------------------------------------------------
  // fld_create
  //
  // FIXME, should be consolidated with psc_checks.c, and probably other places

  static struct psc_mfields *
  fld_create(struct psc *psc, const char *name)
  {
    auto mflds = PscMfields<Mfields>::create(psc_comm(psc), psc->grid(), 1, psc->ibn);
    psc_mfields_set_comp_name(mflds.mflds(), 0, name);

    return mflds.mflds();
  }

  // ----------------------------------------------------------------------
  // setup

  static void
  setup(struct psc_marder *marder)
  {
    marder->bnd = psc_bnd_create(psc_marder_comm(marder));
    psc_bnd_set_name(marder->bnd, "marder_bnd");
    psc_bnd_set_type(marder->bnd, Mfields_traits<Mfields>::name);
    psc_bnd_set_psc(marder->bnd, ppsc);
    psc_bnd_setup(marder->bnd);

    // FIXME, output_fields should be taking care of their own psc_bnd?
    marder->item_div_e = psc_output_fields_item_create(psc_comm(ppsc));
    char name[20];
    snprintf(name, 20, "dive_%s", Mfields_traits<Mfields>::name);
    psc_output_fields_item_set_type(marder->item_div_e, name);
    psc_output_fields_item_set_psc_bnd(marder->item_div_e, marder->bnd);
    psc_output_fields_item_setup(marder->item_div_e);

    marder->item_rho = psc_output_fields_item_create(psc_comm(ppsc));
    auto s = std::string("rho_1st_nc_") + Mparticles_traits<Mparticles_t>::name;
    psc_output_fields_item_set_type(marder->item_rho, s.c_str());
    psc_output_fields_item_set_psc_bnd(marder->item_rho, marder->bnd);
    psc_output_fields_item_setup(marder->item_rho);

    if (marder->dump) {
      struct mrc_io *io = mrc_io_create(psc_comm(ppsc));
      mrc_io_set_type(io, "xdmf_collective");
      mrc_io_set_name(io, "mrc_io_marder");
      mrc_io_set_param_string(io, "basename", "marder");
      mrc_io_set_from_options(io);
      mrc_io_setup(io);

      marder->io = io;
    }
  }

  // ----------------------------------------------------------------------
  // destroy

  static void
  destroy(struct psc_marder *marder)
  {
    psc_output_fields_item_destroy(marder->item_div_e);
    psc_output_fields_item_destroy(marder->item_rho);

    psc_bnd_destroy(marder->bnd);

    if (marder->dump) {
      mrc_io_destroy(marder->io);
    }
  }

  // ----------------------------------------------------------------------
  // correct_patch
  //
  // Do the modified marder correction (See eq.(5, 7, 9, 10) in Mardahl and Verboncoeur, CPC, 1997)

#define define_dxdydz(dx, dy, dz)				\
  int dx _mrc_unused = (ppsc->grid().isInvar(0)) ? 0 : 1;	\
  int dy _mrc_unused = (ppsc->grid().isInvar(1)) ? 0 : 1;	\
  int dz _mrc_unused = (ppsc->grid().isInvar(2)) ? 0 : 1

#define psc_foreach_3d_more(psc, p, ix, iy, iz, l, r) {	\
  int __ilo[3] = { -l[0], -l[1], -l[2] };		\
  int __ihi[3] = { psc->grid().ldims[0] + r[0],		\
		   psc->grid().ldims[1] + r[1],		\
		   psc->grid().ldims[2] + r[2] };	\
  for (int iz = __ilo[2]; iz < __ihi[2]; iz++) {	\
  for (int iy = __ilo[1]; iy < __ihi[1]; iy++) {	\
  for (int ix = __ilo[0]; ix < __ihi[0]; ix++)

#define psc_foreach_3d_more_end			\
  } } }

  static void
  correct_patch(struct psc_marder *marder, fields_t flds, fields_t f, int p, real_t& max_err)
  {
    Fields F(flds), FF(f);
    define_dxdydz(dx, dy, dz);

    // FIXME: how to choose diffusion parameter properly?
    //double deltax = ppsc->patch[p].dx[0];
    double deltay = ppsc->grid().domain.dx[1]; // FIXME double/float
    double deltaz = ppsc->grid().domain.dx[2];
    double inv_sum = 0.;
    int nr_levels;
    const auto& grid = ppsc->grid();
    for (int d = 0; d < 3; d++) {
      if (!grid.isInvar(d)) {
	inv_sum += 1. / sqr(grid.domain.dx[d]);
      }
    }
    double diffusion_max = 1. / 2. / (.5 * ppsc->dt) / inv_sum;
    double diffusion     = diffusion_max * marder->diffusion;

    int l_cc[3] = {0, 0, 0}, r_cc[3] = {0, 0, 0};
    int l_nc[3] = {0, 0, 0}, r_nc[3] = {0, 0, 0};
    for (int d = 0; d < 3; d++) {
      if (grid.bc.fld_lo[d] == BND_FLD_CONDUCTING_WALL &&
	  psc_at_boundary_lo(ppsc, p, d)) {
	l_cc[d] = -1;
	l_nc[d] = -1;
      }
      if (grid.bc.fld_hi[d] == BND_FLD_CONDUCTING_WALL &&
	  psc_at_boundary_hi(ppsc, p, d)) {
	r_cc[d] = -1;
	r_nc[d] = 0;
      }
    }

#if 0
    psc_foreach_3d_more(ppsc, p, ix, iy, iz, l, r) {
      // FIXME: F3 correct?
      F(EX, ix,iy,iz) += 
	(FF(DIVE_MARDER, ix+dx,iy,iz) - FF(DIVE_MARDER, ix,iy,iz))
	* .5 * ppsc->dt * diffusion / deltax;
      F(EY, ix,iy,iz) += 
	(FF(DIVE_MARDER, ix,iy+dy,iz) - FF(DIVE_MARDER, ix,iy,iz))
	* .5 * ppsc->dt * diffusion / deltay;
      F(EZ, ix,iy,iz) += 
	(FF(DIVE_MARDER, ix,iy,iz+dz) - FF(DIVE_MARDER, ix,iy,iz))
	* .5 * ppsc->dt * diffusion / deltaz;
    } psc_foreach_3d_more_end;
#endif

    assert(ppsc->grid().isInvar(0));

    {
      int l[3] = { l_nc[0], l_cc[1], l_nc[2] };
      int r[3] = { r_nc[0], r_cc[1], r_nc[2] };
      psc_foreach_3d_more(ppsc, p, ix, iy, iz, l, r) {
	max_err = std::max(max_err, std::abs(FF(0, ix,iy,iz)));
	F(EY, ix,iy,iz) += 
	  (FF(0, ix,iy+dy,iz) - FF(0, ix,iy,iz))
	  * .5 * ppsc->dt * diffusion / deltay;
      } psc_foreach_3d_more_end;
    }

    {
      int l[3] = { l_nc[0], l_nc[1], l_cc[2] };
      int r[3] = { r_nc[0], r_nc[1], r_cc[2] };
      psc_foreach_3d_more(ppsc, p, ix, iy, iz, l, r) {
	F(EZ, ix,iy,iz) += 
	  (FF(0, ix,iy,iz+dz) - FF(0, ix,iy,iz))
	  * .5 * ppsc->dt * diffusion / deltaz;
      } psc_foreach_3d_more_end;
    }
  }

#undef psc_foreach_3d_more
#undef psc_foreach_3d_more_end

  // ----------------------------------------------------------------------
  // correct

  static void
  correct(struct psc_marder *marder, struct psc_mfields *_mflds_base,
	  struct psc_mfields *div_e)
  {
    auto mflds_base = PscMfieldsBase{_mflds_base};
    auto& mf = mflds_base->get_as<Mfields>(EX, EX + 3);
    auto& mf_div_e = *PscMfields<Mfields>{div_e}.sub();

    real_t max_err = 0.;
    for (int p = 0; p < mf_div_e.n_patches(); p++) {
      correct_patch(marder, mf[p], mf_div_e[p], p, max_err);
    }

    MPI_Allreduce(MPI_IN_PLACE, &max_err, 1, Mfields_traits<Mfields>::mpi_dtype(), MPI_MAX, psc_marder_comm(marder));
    mpi_printf(psc_marder_comm(marder), "marder: err %g\n", max_err);

    mflds_base->put_as(mf, EX, EX + 3);
  }

  static void run(struct psc_marder *marder, PscMfieldsBase mflds_base,
		  PscMparticlesBase mprts_base)
  {
    PscFieldsItemBase item_rho(marder->item_rho);
    PscFieldsItemBase item_div_e(marder->item_div_e);
    item_rho(mflds_base, mprts_base);

    // need to fill ghost cells first (should be unnecessary with only variant 1) FIXME
    auto bnd = PscBndBase(ppsc->bnd);
    bnd.fill_ghosts(mflds_base, EX, EX+3);

    for (int i = 0; i < marder->loop; i++) {
      marder_calc_aid_fields(marder, mflds_base, mprts_base);
      correct(marder, mflds_base.mflds(), item_div_e->mres().mflds());
      auto bnd = PscBndBase(ppsc->bnd);
      bnd.fill_ghosts(mflds_base, EX, EX+3);
    }
  }

  // ----------------------------------------------------------------------
  // marder_calc_aid_fields

  static void marder_calc_aid_fields(struct psc_marder *marder, 
				     PscMfieldsBase mflds_base, PscMparticlesBase mprts_base)
  {
    PscFieldsItemBase item_div_e(marder->item_div_e);
    PscFieldsItemBase item_rho(marder->item_rho);
    item_div_e(mflds_base, mprts_base); // FIXME, should accept NULL for particles
    
    if (marder->dump) {
      static int cnt;
      mrc_io_open(marder->io, "w", cnt, cnt);//ppsc->timestep, ppsc->timestep * ppsc->dt);
      cnt++;
      psc_mfields_write_as_mrc_fld(item_rho->mres().mflds(), marder->io);
      psc_mfields_write_as_mrc_fld(item_div_e->mres().mflds(), marder->io);
      mrc_io_close(marder->io);
    }
    
    item_div_e->mres()->axpy_comp(0, -1., *item_rho->mres().sub(), 0);
    // FIXME, why is this necessary?
    auto bnd = PscBndBase(marder->bnd);
    bnd.fill_ghosts(item_div_e->mres(), 0, 1);
  }
};

