
#include "art/Persistency/Common/Wrapper.h"

#include "MuCSData.h"
#include "MuCSRecoData.h"
#include "MuCSDTOffset.h"

template class std::vector<MuCS::MuCSDTOffset>;

template class art::Wrapper<MuCS::MuCSDTOffset>;

template class art::Wrapper<std::vector<MuCS::MuCSDTOffset> >;

template class std::vector<MuCS::MuCSData>;

template class art::Wrapper< MuCS::MuCSData >; 

template class art::Wrapper< std::vector<MuCS::MuCSData>  >;

template class std::vector<MuCS::MuCSRecoData>;

template class art::Wrapper< MuCS::MuCSRecoData >; 

template class art::Wrapper< std::vector<MuCS::MuCSRecoData>  >;
