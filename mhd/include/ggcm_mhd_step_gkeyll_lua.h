
#ifndef GGCM_MHD_STEP_GKEYLL_LUA_H
#define GGCM_MHD_STEP_GKEYLL_LUA_H

void ggcm_mhd_step_gkeyll_setup_flds_lua(struct ggcm_mhd *mhd);
void ggcm_mhd_step_gkeyll_lua_setup(void **lua_state_ptr, const char *script,
    struct ggcm_mhd *mhd, struct mrc_fld *qFlds[]);
void ggcm_mhd_step_gkeyll_lua_run(void *lua_state, struct ggcm_mhd *mhd, struct mrc_fld *fld);
void ggcm_mhd_step_gkeyll_lua_destroy(void *lua_state, struct ggcm_mhd *mhd);

#endif

