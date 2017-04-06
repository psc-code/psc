
#include <mpi.h>

extern "C" {

#include "ggcm_mhd_step_gkeyll_lua.h"

#include <ggcm_mhd_private.h>
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_crds_gen.h>
#include <ggcm_mhd_gkeyll.h>
#include <string.h>

}

#include <LcLogStream.h>
#include <LcLogger.h>
#include <LcLuaState.h>
#include <LcRegisterModules.h>
#include <LcStreamHandler.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>

static int ggcm_mhd_reduce_double_min_lua(lua_State *L) {
  double var = lua_tonumber(L, -1);
  struct ggcm_mhd *mhd = (struct ggcm_mhd *) lua_touserdata(L, -2);
  double temp_var = var;
  MPI_Allreduce(&temp_var, &var, 1, MPI_DOUBLE, MPI_MIN, mrc_domain_comm(mhd->domain));
  lua_pushnumber(L, var);
  return 1; 
}

static int ggcm_mhd_reduce_boolean_lua(lua_State *L) {
  bool val_wanted = lua_toboolean(L, -1);
  bool var = lua_toboolean(L, -2);
  struct ggcm_mhd *mhd = (struct ggcm_mhd *) lua_touserdata(L, -3);
  bool temp_var = var;
 if (val_wanted)
  MPI_Allreduce(&temp_var, &var, 1, MPI_C_BOOL, MPI_LAND, mrc_domain_comm(mhd->domain));
 else
  MPI_Allreduce(&temp_var, &var, 1, MPI_C_BOOL, MPI_LOR, mrc_domain_comm(mhd->domain));
  lua_pushboolean(L, var);
  return 1; 
}

static int ggcm_mhd_get_3d_fld_lua(lua_State *L) {
  int nr_comps = lua_tointeger(L, -1);
  struct ggcm_mhd *mhd = (struct ggcm_mhd *) lua_touserdata(L, -2);
  struct mrc_fld *fld = ggcm_mhd_get_3d_fld(mhd, nr_comps);
  mrc_fld_dict_add_int(fld, "mhd_type", MT_GKEYLL);
  lua_pushlightuserdata(L, fld);
  return 1;
}

static int ggcm_mhd_put_3d_fld_lua(lua_State *L) {
  struct mrc_fld *fld = (struct mrc_fld *) lua_touserdata(L, -1);
  struct ggcm_mhd *mhd = (struct ggcm_mhd *) lua_touserdata(L, -2);
  ggcm_mhd_put_3d_fld(mhd, fld);
  return 0;
}

static int mrc_fld_get_arr_lua(lua_State *L) {
  struct mrc_fld *fld = (struct mrc_fld *) lua_touserdata(L, -1);
  lua_pushlightuserdata(L, fld->_nd->arr);
  return 1;
}

static int ggcm_mhd_fill_ghosts_lua (lua_State *L) {
  double bntim = lua_tonumber(L, -1);
  struct mrc_fld *fld = (struct mrc_fld *) lua_touserdata(L, -2);
  struct ggcm_mhd *mhd = (struct ggcm_mhd *) lua_touserdata(L, -3);
  ggcm_mhd_fill_ghosts(mhd, fld, bntim);
  return 0;
}

// fill arr[] with lua array named arr_name
static void
lua_getarray(lua_State *L_temp, const char *arr_name, int nr_fluids, float arr[])
{
  lua_getglobal(L_temp, arr_name);
  for (int s = 0; s < nr_fluids; s++) {
    lua_pushnumber(L_temp, s + 1);
    lua_gettable(L_temp, -2);
    arr[s] = (float)lua_tonumber(L_temp, -1);
    lua_pop(L_temp, 1);
  }
  lua_pop(L_temp, 1);
}

// interface for lua function getCArray to get content an C array
static int
lua_pusharray(lua_State *L) {
  int nr_vals = (int) lua_tonumber(L, -1);
  float *arr = (float *) lua_touserdata(L, -2);
  for (int n = 0; n < nr_vals; n++) {
    lua_pushnumber(L, arr[n]);
  }
  return nr_vals;
}

void
ggcm_mhd_step_gkeyll_setup_flds_lua(struct ggcm_mhd *mhd)
{
  int nr_fluids = mhd->par.gk_nr_fluids;
  int nr_moments = mhd->par.gk_nr_moments;
  int nr_comps = nr_fluids * nr_moments + 8;
  mrc_fld_set_param_int(mhd->fld, "nr_comps", nr_comps);
  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 2);
}

void
ggcm_mhd_step_gkeyll_lua_setup(void **lua_state_ptr, const char *script,
    struct ggcm_mhd *mhd, struct mrc_fld *qFlds[])
{
  *lua_state_ptr = (void *) new Lucee::LuaState;
  Lucee::LuaState L = *((Lucee::LuaState *)(*lua_state_ptr));

  // determine input file
  std::string inpFile = script;

  // create top-level logger
  Lucee::Logger& logger = Lucee::Logger::create("lucee");
  logger.setLevel("debug");

  // create console logger
  Lucee::Logger& conLogger = Lucee::Logger::create("lucee.console");
  conLogger.setLevel("info");
  // create console stream
  Lucee::StreamHandler conStrm(std::cout);
  conStrm.attachToLogger("lucee.console");

  // load lua library: this must be done before loading input file
  Lucee::registerModules(L);

  // add command line options to the top-level module
  static const luaL_Reg topFuncs[] = { {NULL, NULL} };
  luaL_register(L, "Lucee", topFuncs);

  lua_pop(L, 1); // done adding command line stuff

  lua_pushcfunction(L, lua_pusharray);
  lua_setglobal(L, "getCArray");

  lua_pushcfunction(L, ggcm_mhd_reduce_double_min_lua);
  lua_setglobal(L, "ggcm_mhd_reduce_double_min");

  lua_pushcfunction(L, ggcm_mhd_reduce_boolean_lua);
  lua_setglobal(L, "ggcm_mhd_reduce_boolean");

  lua_pushcfunction(L, ggcm_mhd_get_3d_fld_lua);
  lua_setglobal(L, "ggcm_get_3d_fld");

  lua_pushcfunction(L, ggcm_mhd_put_3d_fld_lua);
  lua_setglobal(L, "ggcm_put_3d_fld");

  lua_pushcfunction(L, mrc_fld_get_arr_lua);
  lua_setglobal(L, "mrc_fld_get_arr");

  lua_pushcfunction(L, ggcm_mhd_fill_ghosts_lua);
  lua_setglobal(L, "ggcm_fill_ghosts");

  lua_pushlightuserdata(L, mhd);
  lua_setglobal(L, "ggcm_mhd");

  int rank;
  MPI_Comm_rank(mrc_domain_comm(mhd->domain), &rank);
  lua_pushinteger(L, rank);
  lua_setglobal(L, "rank");

  lua_pushlightuserdata(L, mhd->ymask);
  lua_setglobal(L, "ymask");

  lua_pushlightuserdata(L, mhd->b0);
  lua_setglobal(L, "b0");

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(mhd->domain, &nr_patches);
  assert(nr_patches > 0);
  const int *ldims = patches[0].ldims;
  lua_pushinteger(L, ldims[0]);
  lua_setglobal(L, "mx");
  lua_pushinteger(L, ldims[1]);
  lua_setglobal(L, "my");
  lua_pushinteger(L, ldims[2]);
  lua_setglobal(L, "mz");

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  bool nonuniform = strcmp(mrc_crds_type(crds), "uniform");
  lua_pushboolean(L, nonuniform);
  lua_setglobal(L, "nonuniform");

  lua_pushlightuserdata(L, crds->dcrd_nc[0]->_nd->arr);
  lua_setglobal(L, "crdx");
  lua_pushlightuserdata(L, crds->dcrd_nc[1]->_nd->arr);
  lua_setglobal(L, "crdy");
  lua_pushlightuserdata(L, crds->dcrd_nc[2]->_nd->arr);
  lua_setglobal(L, "crdz");

  double ll[3], lh[3];
  if (nonuniform) {
    // ll, lh not needed
  } else {
    double dx[3];
    mrc_crds_get_dx_base(crds, dx);
    ll[0] = MRC_DMCRDX(crds, 0, 0) - .5 * dx[0];
    ll[1] = MRC_DMCRDY(crds, 0, 0) - .5 * dx[1];
    ll[2] = MRC_DMCRDZ(crds, 0, 0) - .5 * dx[2];
    lh[0] = MRC_DMCRDX(crds, ldims[0]-1, 0) + .5 * dx[0];
    lh[1] = MRC_DMCRDY(crds, ldims[1]-1, 0) + .5 * dx[1];
    lh[2] = MRC_DMCRDZ(crds, ldims[2]-1, 0) + .5 * dx[2];
  }
  lua_pushnumber(L, ll[0]);
  lua_setglobal(L, "lx");
  lua_pushnumber(L, ll[1]);
  lua_setglobal(L, "ly");
  lua_pushnumber(L, ll[2]);
  lua_setglobal(L, "lz");
  lua_pushnumber(L, lh[0]);
  lua_setglobal(L, "hx");
  lua_pushnumber(L, lh[1]);
  lua_setglobal(L, "hy");
  lua_pushnumber(L, lh[2]);
  lua_setglobal(L, "hz");

  double l[3], h[3];
  mrc_crds_get_param_double3(mrc_domain_get_crds(mhd->domain), "l", l);
  mrc_crds_get_param_double3(mrc_domain_get_crds(mhd->domain), "h", h);
  lua_pushnumber(L, l[0]);
  lua_setglobal(L, "lxg");
  lua_pushnumber(L, l[1]);
  lua_setglobal(L, "lyg");
  lua_pushnumber(L, l[2]);
  lua_setglobal(L, "lzg");
  lua_pushnumber(L, h[0]);
  lua_setglobal(L, "hxg");
  lua_pushnumber(L, h[1]);
  lua_setglobal(L, "hyg");
  lua_pushnumber(L, h[2]);
  lua_setglobal(L, "hzg");

  lua_pushnumber(L, mhd->par.gamm);
  lua_setglobal(L, "gasGamma");

  lua_pushnumber(L, mhd->par.gk_speed_of_light);
  lua_setglobal(L, "lightSpeed");
  lua_pushnumber(L, mhd->par.gk_nr_fluids);
  lua_setglobal(L, "nr_fluids");
  lua_pushnumber(L, mhd->par.gk_nr_moments);
  lua_setglobal(L, "nr_moments");
  lua_pushlightuserdata(L, mhd->par.gk_charge.vals);
  lua_setglobal(L, "charge_array");
  lua_pushlightuserdata(L, mhd->par.gk_mass.vals);
  lua_setglobal(L, "mass_array");

  lua_pushnumber(L, mhd->par.thx);
  lua_setglobal(L, "cfl");

  lua_pushlightuserdata(L, mhd->fld);
  lua_setglobal(L, "qFld");

  if (qFlds[0]) {
    lua_pushlightuserdata(L, qFlds[0]);
    lua_setglobal(L, "qFldX");
  }

  if (qFlds[1]) {
    lua_pushlightuserdata(L, qFlds[1]);
    lua_setglobal(L, "qFldY");
  }

  if (qFlds[2]) {
    lua_pushlightuserdata(L, qFlds[2]);
    lua_setglobal(L, "qFldZ");
  }

  if (luaL_loadfile(L, inpFile.c_str())) {
    std::cerr << "Error loading file: " << inpFile << std::endl;
    std::string err(lua_tostring(L, -1));
    lua_pop(L, 1);
    std::cerr << err << std::endl;
    exit(1);
  }
  
  int nargs = 0;
  int nrets = 0;

  if (lua_pcall(L, nargs, nrets, 0)) {
    std::cerr << "Error in input file: " << inpFile << std::endl;
    std::string err(lua_tostring(L, -1));
    lua_pop(L, 1);
    std::cerr << err << std::endl;
    exit(1);
  }
}

void
ggcm_mhd_step_gkeyll_lua_run(void *lua_state,
    struct ggcm_mhd *mhd, struct mrc_fld *fld)
{
  if (!lua_state) return;
  Lucee::LuaState L = *((Lucee::LuaState *)lua_state);
  lua_getglobal(L, "runTimeStep");
  lua_pushnumber(L, mhd->dt_code);
  lua_pushnumber(L, mhd->time_code);
  lua_pushinteger(L, mhd->istep);
  lua_pushlightuserdata(L, fld->_nd->arr);

  if (lua_pcall(L, 4, 1, 0)) {
    std::cerr << "LUA Error:" << std::endl;
    std::string err(lua_tostring(L, -1));
    lua_pop(L, 1);
    std::cerr << err << std::endl;
    exit(1);
  }
  mhd->dt_code = lua_tonumber(L,-1);

  float myDt = mhd->dt_code;
  MPI_Allreduce(&myDt, &mhd->dt_code, 1, MPI_FLOAT, MPI_MIN, mrc_domain_comm(mhd->domain));

  lua_pop(L,1);
}

void
ggcm_mhd_step_gkeyll_lua_destroy(void *lua_state, struct ggcm_mhd *mhd)
{
  if (!lua_state) return;
  Lucee::LuaState L = *((Lucee::LuaState *)lua_state);
  lua_getglobal(L, "finalize");
  if (lua_pcall(L, 0, 0, 0)) {
    std::cerr << "LUA Error:" << std::endl;
    std::string err(lua_tostring(L, -1));
    lua_pop(L, 1);
    std::cerr << err << std::endl;
    exit(1);
  }
}

