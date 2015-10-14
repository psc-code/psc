
#include <mpi.h>

extern "C" {

#include "ggcm_mhd_ic_gkeyll_lua.h"

#include <ggcm_mhd_private.h>
#include <mrc_fld.h>
#include <mrc_domain.h>

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

void
ggcm_mhd_ic_gkeyll_lua_run(const char *script, const char *script_common,
    struct ggcm_mhd *mhd, struct mrc_fld *fld)
{
  // determine input file
  std::string inpFile = script;

  // create output prefix
  std::string outPrefix;
  // use input file name sans the .lua extension
  std::string snm = inpFile;
  unsigned trunc = inpFile.find_last_of(".", snm.size());
  if (trunc > 0)
    snm.erase(trunc, snm.size());
  outPrefix = snm;

  bool isRestarting = false;
  int rFrame = 0;

  // load input file using Lua
  Lucee::LuaState L;

  // load lua library: this must be done before loading input file
  Lucee::registerModules(L);

  // add command line options to the top-level module
  static const luaL_Reg topFuncs[] = { {NULL, NULL} };
  luaL_register(L, "Lucee", topFuncs);

  lua_pop(L, 1); // done adding command line stuff

  if (luaL_loadfile(L, inpFile.c_str())) {
    std::cerr << "Error loading file: " << inpFile << std::endl;
    std::string err(lua_tostring(L, -1));
    lua_pop(L, 1);
    std::cerr << err << std::endl;
    exit(1);
  }

  int nargs = 0;

  int rank;
  MPI_Comm_rank(mrc_domain_comm(mhd->domain), &rank);
  lua_pushinteger(L, rank);
  nargs += 1;

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  const int *ldims = mrc_fld_spatial_dims(fld);
  lua_pushinteger(L, ldims[0]);
  lua_pushinteger(L, ldims[1]);
  lua_pushinteger(L, ldims[2]);
  nargs += 3;

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  double dx[3];
  mrc_crds_get_dx_base(crds, dx);
  double ll[3], lh[3];
  ll[0] = MRC_DMCRDX(crds, 0, 0) - .5 * dx[0];
  ll[1] = MRC_DMCRDY(crds, 0, 0) - .5 * dx[1];
  ll[2] = MRC_DMCRDZ(crds, 0, 0) - .5 * dx[2];
  lh[0] = MRC_DMCRDX(crds, ldims[0]-1, 0) + .5 * dx[0];
  lh[1] = MRC_DMCRDY(crds, ldims[1]-1, 0) + .5 * dx[1];
  lh[2] = MRC_DMCRDZ(crds, ldims[2]-1, 0) + .5 * dx[2];
  lua_pushnumber(L, ll[0]);
  lua_pushnumber(L, ll[1]);
  lua_pushnumber(L, ll[2]);
  nargs += 3;
  lua_pushnumber(L, lh[0]);
  lua_pushnumber(L, lh[1]);
  lua_pushnumber(L, lh[2]);
  nargs += 3;

  double l[3], h[3];
  mrc_crds_get_param_double3(mrc_domain_get_crds(mhd->domain), "l", l);
  mrc_crds_get_param_double3(mrc_domain_get_crds(mhd->domain), "h", h);
  lua_pushnumber(L, l[0]);
  lua_pushnumber(L, l[1]);
  lua_pushnumber(L, l[2]);
  nargs += 3;
  lua_pushnumber(L, h[0]);
  lua_pushnumber(L, h[1]);
  lua_pushnumber(L, h[2]);
  nargs += 3;

  lua_pushnumber(L, mhd->par.gamm);
  nargs += 1;

  lua_pushstring(L, script_common);
  nargs++;

  lua_pushlightuserdata(L, fld->_arr);
  nargs++;

  if (lua_pcall(L, nargs, 0, 0)) {
    std::cerr << "Error parsing input file: " << inpFile << std::endl;
    std::string err(lua_tostring(L, -1));
    lua_pop(L, 1);
    std::cerr << err << std::endl;
    exit(1);
  }
}
