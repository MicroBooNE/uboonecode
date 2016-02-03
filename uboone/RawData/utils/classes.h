#include "art/Persistency/Common/Wrapper.h"

#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"

template class std::vector< bool >;
template class art::Wrapper< raw::ubdaqSoftwareTriggerData >;
template class std::pair< std::string, bool >;
template class std::vector< std::pair< std::string, bool > >;
