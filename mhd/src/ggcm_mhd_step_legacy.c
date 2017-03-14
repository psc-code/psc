
// ----------------------------------------------------------------------
// ggcm_mhd_step_legacy_setup_flds

static void
ggcm_mhd_step_legacy_setup_flds(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd *mhd = step->mhd;

  mrc_fld_set_type(mhd->fld, FLD_TYPE);
  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 2);
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_SCONS_FC_GGCM);
  mrc_fld_set_param_int(mhd->fld, "nr_comps", _NR_FLDS);
}
