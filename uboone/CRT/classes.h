/**
  \defgroup crt All things Cosmic Ray Tagger related
**/

#include "canvas/Persistency/Common/Wrapper.h"
#include "uboone/CRT/CRTData.hh"
#include <vector>

template class std::vector<crt::CRTData>;
template class art::Wrapper< std::vector<crt::CRTData> >;
template class art::Wrapper< crt::CRTData >;
