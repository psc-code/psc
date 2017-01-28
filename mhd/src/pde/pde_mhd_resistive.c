
// ----------------------------------------------------------------------
// mhd_add_resistive_flux_const

static void
mhd_add_resistive_flux_const(fld1d_state_t F, fld1d_state_t V, int ib, int ie)
{
  for (int i = ib; i < ie; i++) {
    assert(!s_opt_background);
    // FIXME, we do need Btotal here

    mrc_fld_data_t jy = F1V(s_aux.j, 1, i);
    mrc_fld_data_t jz = F1V(s_aux.j, 2, i);
    mrc_fld_data_t by = .5f * (F1S(V, BY, i-1) + F1S(V, BY, i));
    mrc_fld_data_t bz = .5f * (F1S(V, BZ, i-1) + F1S(V, BZ, i));
    
    mrc_fld_data_t flux_BY = - s_eta * jz;
    mrc_fld_data_t flux_BZ =   s_eta * jy;
    mrc_fld_data_t flux_EE = by * flux_BY + bz * flux_BZ;

    F1S(F, BY, i) += flux_BY;
    F1S(F, BZ, i) += flux_BZ;
    F1S(F, EE, i) += flux_EE;
  }
}

// ----------------------------------------------------------------------
// mhd_add_resistive_flux

static void
mhd_add_resistive_flux(fld1d_state_t F, fld1d_state_t V, int ib, int ie)
{
  if (s_opt_resistivity == OPT_RESISTIVITY_NONE) {
  } else if (s_opt_resistivity == OPT_RESISTIVITY_CONST) {
    mhd_add_resistive_flux_const(F, V, ib, ie);
  } else {
    assert(0);
  }
}
