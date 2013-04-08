// ======================================================================
//
// LArRawInputSourceUBooNE_source.cc
//
// ======================================================================

#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/Source.h"
#include "RawData/utils/LArRawInputDriverUBooNE.h"

namespace lris {
  typedef art::Source<LArRawInputDriverUBooNE> LArRawInputSourceUBooNE;
}

DEFINE_ART_INPUT_SOURCE(lris::LArRawInputSourceUBooNE);
