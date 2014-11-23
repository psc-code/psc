
static void
ggcm_mhd_step_legacy_setup_flds(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd *mhd = step->mhd;

  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 2);
}
