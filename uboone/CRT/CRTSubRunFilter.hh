/**
  \class CRTSubRunFilter
  \ingroup CRT
  \brief For a given subrun, uses an input text file to reduce the CRT data.
  \author (last to touch it) $Author: bv $
  \version $Revision: 0.1 $
  \date $Date: 2016/12/03 14:16:20 $
  \todo : Improve detailed description.
  \todo : Change Author version and date when complete
  \todo : Fill in methods
**/

#ifndef CRTSubRunFilter_hh_
#define CRTSubRunFilter_hh_

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include <string>
#include <fstream>

namespace crt{
  class CRTSubRunFilter : public art::EDFilter {
    std:string fSubRunFileName;  /// Filename holder for fhicl parameter
    std::ifstream fSubRunFile;   /// Text file holder
  public:
    /// Default Constructor
    explicit CRTSubRunFilter(fhicl::ParameterSet const & p);
    /// Default Destructor
    ~CRTSubRunFilter();
    /// At beginning of subrun, opens the subrun text file
    void beginSubRun(art::SubRun& subrun);
    /// Per each event filter if not in time window
    bool filter(art::Event & e) override;
    /// Close the subrun file
    void endSubRun(art::SubRun& subrun);
  };
}

#endif