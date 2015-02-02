#include "PedestalRetrievalAlg.h"

#include <exception>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "Geometry/Geometry.h"


namespace dtbse{
  //==============================================================
  //Constructor
  //==============================================================
  PedestalRetrievalAlg::PedestalRetrievalAlg(fhicl::ParameterSet const& p) {
    this->reconfigure(p);
  }

  void PedestalRetrievalAlg::reconfigure(fhicl::ParameterSet const& p) {

    fDefaultCollectionPedMean = p.get<float>("CollectionPedMean", 400);
    fDefaultCollectionPedRMS  = p.get<float>("CollectionPedRMS", 0.3);

    fDefaultInductionPedMean  = p.get<float>("InductionPedMean", 2048);
    fDefaultInductionPedRMS   = p.get<float>("InductionPedRMS", 0.3);
  }

  //==============================================================
  //Get pedestal information
  //==============================================================
  void PedestalRetrievalAlg::GetPedestal(unsigned int channel, float& pedmean, float& pedrms, unsigned long timestamp /*=0*/) const {

    //If the timestamp is default, use the default pedestal values
    if (timestamp==0) {
      art::ServiceHandle<geo::Geometry> geo;
      geo::SigType_t sigtype = geo->SignalType(channel);
      if ( sigtype == geo::kCollection ) {
	pedmean = fDefaultCollectionPedMean;
	pedrms  = fDefaultCollectionPedRMS;
	return;
      }
      else if ( sigtype == geo::kInduction ) {
	pedmean = fDefaultInductionPedMean;
	pedrms  = fDefaultInductionPedRMS;
	return;
      }
      else { //throw exception
	throw cet::exception("PedestalRetrievalAlg")
              <<"  Unknown geometry signal wire of type "<<sigtype<<std::endl;
      }   
      return;
    }
    else { //retrieve pedestal from database

      if (!RetrieveFromDB(channel, pedmean, pedrms, timestamp)) {
	throw cet::exception("PedestalRetrievalAlg")
              <<"  Channel: "<<channel<<", timestamp: "<<timestamp
	      <<":  Failed to retrieve pedestal information from database"<<std::endl;
      }


      return;
    }
  }//end GetPedestal

  void PedestalRetrievalAlg::GetPedestalMean(unsigned int channel, float& pedmean, unsigned long timestamp /*=0*/) const {

    float dummy_rms = 0.0;
    GetPedestal(channel, pedmean, dummy_rms, timestamp);
    return;
  }

  void PedestalRetrievalAlg::GetPedestalRMS(unsigned int channel, float& pedrms, unsigned long timestamp /*=0*/) const {

    float dummy_mean = 0.0;
    GetPedestal(channel, dummy_mean, pedrms, timestamp);
    return;
  }



  //==============================================================
  //Retrieve pedestal from database.  
  //Return value is true if successful.  Right now we don't have a 
  //database, so return false
  //==============================================================
  bool PedestalRetrievalAlg::RetrieveFromDB(unsigned int channel, float& pedmean, float& pedrms, unsigned long timestamp) const {

    pedmean = 0.0;
    pedrms  = 0.0;
    return false;
  }
}//end namespace dtbse
