#ifndef _GETFOM_H
#define _GETFOM_H

#include "datatypes/ub_BeamHeader.h"
#include "datatypes/ub_BeamData.h"
#include <string>

namespace bmd
{
  float getFOM(std::string beam, const  gov::fnal::uboone::datatypes::ub_BeamHeader& bh, const std::vector<gov::fnal::uboone::datatypes::ub_BeamData>& bd);
}

#endif
