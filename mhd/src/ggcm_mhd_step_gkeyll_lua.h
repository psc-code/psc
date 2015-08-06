
#ifndef GGCM_MHD_STEP_GKEYLL_LUA_H
#define GGCM_MHD_STEP_GKEYLL_LUA_H

void ggcm_mhd_step_gkeyll_setup_flds_lua(struct mrc_fld *fld, const char *script_common); 
void ggcm_mhd_step_gkeyll_lua_setup(const char *script, const char *script_common,
    struct ggcm_mhd *mhd, struct mrc_fld *fld);
void ggcm_mhd_step_gkeyll_lua_run(struct ggcm_mhd *mhd, struct mrc_fld *fld);
void ggcm_mhd_step_gkeyll_lua_destroy(struct ggcm_mhd *mhd);

#endif

