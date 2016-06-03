#include "canvas/Persistency/Common/Wrapper.h"

#include "uboone/OpticalDetectorAna/OpticalSubEvents/subevent_algo/SubEvent.hh"
#include "uboone/OpticalDetectorAna/OpticalSubEvents/subevent_algo/Flash.hh"
#include "uboone/OpticalDetectorAna/OpticalSubEvents/subevent_algo/SubEventList.hh"
#include "uboone/OpticalDetectorAna/OpticalSubEvents/subevent_algo/FlashList.hh"

template class std::vector< subevent::SubEvent >;
template class art::Wrapper< std::vector< subevent::SubEvent > >;
template class std::vector< subevent::Flash >;
template class art::Wrapper< subevent::FlashList >;
