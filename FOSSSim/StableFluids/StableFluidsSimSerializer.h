#ifndef __STABLE_FLUIDS_SERIALIZER_H__
#define __STABLE_FLUIDS_SERIALIZER_H__

#include <fstream>
#include <iostream>

#include "StableFluidsSim.h"
#include "FOSSSim/StringUtilities.h"

class StableFluidsSimSerializer
{
public:
  void serializeScene( StableFluidsSim& scene, std::ofstream& outputstream ) const;

  void loadScene( StableFluidsSim& scene, std::ifstream& inputstream ) const;
};

#endif
