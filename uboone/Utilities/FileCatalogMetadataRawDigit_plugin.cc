////////////////////////////////////////////////////////////////////////
// Name:  FileCatalogMetadataRawDigit_service.cc.  
//
// Purpose:  Plugin to generate metadata for counting events containing
//           RawDigit data products in an output stream.
//
// FCL parameters:
//
//   RawDigitLabel - Specify module label of RawDigits to count.
//
//
// Created:  10-Oct-2015,  H. Greenlee
//
////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "lardataobj/RawData/RawDigit.h"
#include "art/Framework/Core/FileCatalogMetadataPlugin.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace util {

  // Class declaration.

  class FileCatalogMetadataRawDigit : public art::FileCatalogMetadataPlugin
  {
  public:

    // Constructor, destructor.

    FileCatalogMetadataRawDigit(fhicl::ParameterSet const & pset);
    virtual ~FileCatalogMetadataRawDigit() = default;

    // Overrides.

    void collectMetadata(const art::Event& evt);
    collection_type produceMetadata();

  private:

    // Data members.

    // Fcl parameters.

    std::string fRawDigitLabel;

    // Statistics.

    int fNevents;
    int fTPCevents;
  };
}

//--------------------------------------------------------------------
// Constructor.

util::FileCatalogMetadataRawDigit::FileCatalogMetadataRawDigit(fhicl::ParameterSet const & pset):
  FileCatalogMetadataPlugin(pset),
  fRawDigitLabel(pset.get<std::string>("RawDigitLabel")),
  fNevents(0),
  fTPCevents(0)
{
  mf::LogInfo("FileCatalogMetadataRawDigit") 
    << "FileCatalogMetadataRawDigit configured to count events with RawDigit label " 
    << fRawDigitLabel;
}

//--------------------------------------------------------------------
// Method to collect statistics from event.

void util::FileCatalogMetadataRawDigit::collectMetadata(const art::Event& evt)
{
  ++fNevents;

  // Count events that contain a RawDigit data product with the specified label.

  art::Handle<std::vector<raw::RawDigit> > h;
  evt.getByLabel(fRawDigitLabel, h);
  if(h.isValid())
    ++fTPCevents;
}

//--------------------------------------------------------------------
// Method to generate metadata.

art::FileCatalogMetadataPlugin::collection_type util::FileCatalogMetadataRawDigit::produceMetadata()
{
  // Return value (intially empty collection).

  art::FileCatalogMetadataPlugin::collection_type result;

  // Send an informational message using the message facility.

  float ratio = 0.;
  if(fNevents != 0)
    ratio = float(fTPCevents) / float(fNevents);

  mf::LogInfo("FileCatalogMetadataRawDigit")
    << "FileCatalogMetadataRawDigit produce metadata:\n"
    << fTPCevents << " events with RawDigits\n"
    << fNevents << " total events\n"
    << "Fraction of events with RawDigits = " << ratio;

  // Fill metadata collection.

  std::ostringstream ostr;
  ostr << fTPCevents;
  result.emplace_back("filter.tpc_event_count", ostr.str());

  // Reset.

  fNevents = 0;
  fTPCevents = 0;

  // Done.

  return result;
}

DEFINE_ART_FILECATALOGMETADATA_PLUGIN(util::FileCatalogMetadataRawDigit)
