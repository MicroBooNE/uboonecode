
#include "art/Persistency/Common/Wrapper.h"

#include "MuCSData.h"
#include "MuCSRecoData.h"

template class std::vector<MuCS::MuCSData>;

template class art::Wrapper< MuCS::MuCSData >; 

template class art::Wrapper< std::vector<MuCS::MuCSData>  >;

template class std::vector<MuCS::MuCSRecoData>;

template class art::Wrapper< MuCS::MuCSRecoData >; 

template class art::Wrapper< std::vector<MuCS::MuCSRecoData>  >;
