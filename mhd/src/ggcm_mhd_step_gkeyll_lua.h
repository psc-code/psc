
#ifndef GGCM_MHD_STEP_GKEYLL_LUA_H
#define GGCM_MHD_STEP_GKEYLL_LUA_H

void ggcm_mhd_step_gkeyll_setup_flds_lua(const char*script_domain, 
            int *nr_comps, int *nr_ghosts, int *nr_moments, int *nr_fluids);
void ggcm_mhd_step_gkeyll_lua_setup(const char *script, const char *script_common,
    struct ggcm_mhd *mhd, struct mrc_fld *fld);
void ggcm_mhd_step_gkeyll_lua_run(struct ggcm_mhd *mhd, struct mrc_fld *fld);

#endif

