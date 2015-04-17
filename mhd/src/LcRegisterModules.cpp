/**
 * @file	LcRegisterModules.cpp
 *
 * @brief	Function for registering all modules.
 */

// lucee includes
#include <LcDataStructIfc.h>
#include <LcDataStructRegistry.h>
#include <LcDecompRegionCalcIfc.h>
#include <LcDecompRegionCalcRegistry.h>
#include <LcGridIfc.h>
#include <LcGridRegistry.h>
#include <LcHyperEquation.h>
#include <LcHyperEquationRegistry.h>
#include <LcLibRegistry.h>
#include <LcLuaModuleRegistry.h>
#include <LcLuceeMod.h>
#include <LcPoissonBracketEquation.h>
#include <LcPoissonBracketEquationRegistry.h>
#include <LcProtoSolverRegistry.h>
#include <LcRegisterModules.h>
#include <LcRteRegistry.h>
#include <LcSolverIfc.h>
#include <LcSolverRegistry.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Register all modules in Lucee.
 *
 * @param L Lua state object in which modules should be registered.
 */
  void registerModules(Lucee::LuaState& L)
  {
// register objects
    Lucee::registerProtoSolverObjects(L);
    Lucee::registerSolverObjects(L);
    Lucee::registerHyperEquationsObjects(L);
    Lucee::registerPoissonBracketEquationObjects(L);
    Lucee::registerRteObjects(L);
    Lucee::registerGridObjects(L);
    Lucee::registerDataStructObjects(L);
    Lucee::registerLibObjects(L);
    Lucee::registerDecompRegionCalc(L);

// register modules into Lua
    Lucee::LuaModuleRegistry<Lucee::SolverIfc>::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::HyperEquation>::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::PoissonBracketEquation>::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::GridIfc>::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::DataStructIfc>::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::UpdaterIfc>::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::LuceeMod>::registerModule(L);

    Lucee::LuaModuleRegistry<Lucee::DecompRegionCalcIfc<1> >::registerModule(L);    
    Lucee::LuaModuleRegistry<Lucee::DecompRegionCalcIfc<2> >::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::DecompRegionCalcIfc<3> >::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::DecompRegionCalcIfc<4> >::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::DecompRegionCalcIfc<5> >::registerModule(L);
  }
}
