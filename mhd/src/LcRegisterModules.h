/**
 * @file	LcRegisterModules.h
 *
 * @brief	Function for registering all modules.
 */

#ifndef LC_REGISTER_MODULES_H
#define LC_REGISTER_MODULES_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuaState.h>

namespace Lucee
{
  void registerModules(Lucee::LuaState& L);
}

#endif // LC_REGISTER_MODULES_H
