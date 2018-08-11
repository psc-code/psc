
#include "psc.h"
#include "fields3d.hxx"

#include <mrc_params.h>
#include <mrc_profile.h>
#include <mrc_io.h>
#include <stdlib.h>
#include <string.h>
#include <list>

std::list<MfieldsBase*> MfieldsBase::instances;

void MfieldsBase::convert(MfieldsBase& mf_from, MfieldsBase& mf_to, int mb, int me)
{
  // FIXME, implementing == wouldn't hurt
  assert(&mf_from.grid() == &mf_to.grid());
  
  auto convert_to = mf_from.convert_to().find(std::type_index(typeid(mf_to)));
  if (convert_to != mf_from.convert_to().cend()) {
    convert_to->second(mf_from, mf_to, mb, me);
    return;
  }
  
  auto convert_from = mf_to.convert_from().find(std::type_index(typeid(mf_from)));
  if (convert_from != mf_to.convert_from().cend()) {
    convert_from->second(mf_to, mf_from, mb, me);
    return;
  }

  fprintf(stderr, "ERROR: no conversion known from %s to %s!\n",
	  typeid(mf_from).name(), typeid(mf_to).name());
  assert(0);
}

#if 1
std::list<MfieldsStateBase*> MfieldsStateBase::instances;

void MfieldsStateBase::convert(MfieldsStateBase& mf_from, MfieldsStateBase& mf_to, int mb, int me)
{
  // FIXME, implementing == wouldn't hurt
  assert(&mf_from.grid() == &mf_to.grid());
  
  auto convert_to = mf_from.convert_to().find(std::type_index(typeid(mf_to)));
  if (convert_to != mf_from.convert_to().cend()) {
    convert_to->second(mf_from, mf_to, mb, me);
    return;
  }
  
  auto convert_from = mf_to.convert_from().find(std::type_index(typeid(mf_from)));
  if (convert_from != mf_to.convert_from().cend()) {
    convert_from->second(mf_to, mf_from, mb, me);
    return;
  }

  fprintf(stderr, "ERROR: no conversion known from %s to %s!\n",
	  typeid(mf_from).name(), typeid(mf_to).name());
  assert(0);
}
#endif
