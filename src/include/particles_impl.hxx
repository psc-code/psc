
#pragma once

void MparticlesBase::convert(MparticlesBase& mp_from, MparticlesBase& mp_to)
{
  // FIXME, implementing == wouldn't hurt
  assert(&mp_from.grid() == &mp_to.grid());
  
  auto convert_to = mp_from.convert_to().find(std::type_index(typeid(mp_to)));
  if (convert_to != mp_from.convert_to().cend()) {
    convert_to->second(mp_from, mp_to);
    return;
  }
  
  auto convert_from = mp_to.convert_from().find(std::type_index(typeid(mp_from)));
  if (convert_from != mp_to.convert_from().cend()) {
    convert_from->second(mp_to, mp_from);
    return;
  }

  fprintf(stderr, "ERROR: no conversion known from %s to %s!\n",
	  typeid(mp_from).name(), typeid(mp_to).name());
  assert(0);
}

 

