////////////////////////////////////////////////////////////////////////
// Class:       AnodeCathodePMAlg
//
// Algorithm for quickly identifying anode-to-cathode tracks based on a
// simple pattern-matching approach.
//
// Input:  vector<recob::Hit>
// Output: true if anode-cathode track candidate in event (false if not)
// 
////////////////////////////////////////////////////////////////////////

#include "AnodeCathodePMAlg.h"

//#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/Hit.h"

#include <math.h>
#include <algorithm>

pm::AnodeCathodePMAlg::AnodeCathodePMAlg()
{
}

void pm::AnodeCathodePMAlg::Configure(//fhicl::ParameterSet         const& p,
				      geo::GeometryCore           const& geo,
				      detinfo::DetectorProperties const& detp)
{
  /*
  //grab the wire information
  fPlaneID   = geo::PlaneID(p.get<geo::CryostatID::CryostatID_t>("CryostatID",0),
			    p.get<geo::TPCID::TPCID_t>("TPCID",0),
			    p.get<geo::PlaneID::PlaneID_t>("PlaneID"));
  fStartWire = geo::WireID(fPlaneID,p.get<geo::WireID::WireID_t>("StartWire",0));
  fEndWire = geo::WireID(fPlaneID,p.get<geo::WireID::WireID_t>("EndWire",geo.Nwires(fPlaneID)-1));
  */
  fPlaneID = geo::PlaneID(0,0,2);
  fStartWire = geo::WireID(fPlaneID,0);
  fEndWire = geo::WireID(0,0,2,geo.Nwires(fPlaneID)-1);

  //check wire indices: end > start
  if(fEndWire.Wire < fStartWire.Wire){
    geo::WireID tmp(fEndWire);
    fEndWire = fStartWire;
    fStartWire = tmp;
  }
  //fill nwires to run over. By default, this is all wires in the plane
  fNwires = fEndWire.Wire - fStartWire.Wire + 1;

  //set wires per bin, and set it right
  //fWiresPerBin = p.get<unsigned int>("WiresPerBin");
  fWiresPerBin = 400;
  if(fNwires < fWiresPerBin)
    fWiresPerBin = fNwires;  
  fNbinsWires = (fNwires / fWiresPerBin) + 1;

  //grab time tick info
  fStartTimeTick = 0.0;
  fEndTimeTick   = (float)(detp.NumberTimeSamples());

  //fTimeTicksPerBin = p.get<float>("TimeTicksPerBin");  
  fTimeTicksPerBin = 100;
  fNbinsTimeTicks = (unsigned int)((fEndTimeTick - fStartTimeTick)/fTimeTicksPerBin) + 1;
  
  fPatterns.resize(fNbinsWires*fNbinsWires);
  fHitmap.resize(fNbinsWires,std::vector<float>(fNbinsTimeTicks));

  //fNtimeTicksPattern = (unsigned int)(std::round(p.get<float>("NtimeTicksPattern")/fTimeTicksPerBin));
  fNtimeTicksPattern = (unsigned int)(std::round(4500./fTimeTicksPerBin));

  //fIntegralThreshold = p.get<float>("IntegralThreshold");
  fIntegralThreshold = 50.0;
  //fMatchFractionThreshold = p.get<float>("MatchFractionThreshold",10.0);
  fMatchFractionThreshold = 10.0;
  
  MakePatterns();
}

void pm::AnodeCathodePMAlg::ClearHitmap(){
  for(auto & w : fHitmap)
    std::fill(w.begin(),w.end(),0.0);
}

unsigned int pm::AnodeCathodePMAlg::GetBinTime(float time){
  if(time<fStartTimeTick || time>fEndTimeTick) return -1;
  return (time - fStartTimeTick)/fTimeTicksPerBin;
}

unsigned int pm::AnodeCathodePMAlg::GetBinWire(float wire){
  if(wire<fStartWire.Wire || wire>fEndWire.Wire) return -1;
  return (wire - fStartWire.Wire)/fWiresPerBin;
}

void pm::AnodeCathodePMAlg::MakePatterns(){

  float wire_a,wire_c;
  float slope, intercept, time, end_time;

  for(size_t i_p=0; i_p<fPatterns.size(); ++i_p){

    wire_a = (float)(fStartWire.Wire + (i_p/fNbinsWires)*fWiresPerBin);
    wire_c = (float)(fStartWire.Wire + (i_p%fNbinsWires)*fWiresPerBin);
    slope = (wire_c - wire_a)/(float)(fNtimeTicksPattern);

    for(unsigned int i_s_t=0; i_s_t<(fNbinsTimeTicks-fNtimeTicksPattern); ++i_s_t){
      intercept = wire_a - (fStartTimeTick+i_s_t*fTimeTicksPerBin)*slope;
      
      for(time=fStartTimeTick+i_s_t*fTimeTicksPerBin, end_time = fStartTimeTick+(i_s_t+fNtimeTicksPattern)*fTimeTicksPerBin;
	  time<=end_time;
	  time += 0.1*fTimeTicksPerBin){
	fPatterns[i_p].emplace_back(GetBinWire(wire_a+time*slope-intercept),GetBinTime(time));
      }
    }
  }

}

void pm::AnodeCathodePMAlg::FillHitmap( std::vector<recob::Hit> const& hits)
{
  ClearHitmap();
  for(auto const& hit : hits){
    if(hit.WireID().planeID()!=fPlaneID) continue;
    fHitmap[GetBinWire(hit.WireID())][GetBinTime(hit.PeakTime())] += hit.Integral();
  }
}

float pm::AnodeCathodePMAlg::GetFractionMatched()
{
  unsigned int n_matched=0;
  float max_fraction_matched = -1.0;
  //std::vector< std::pair<size_t,float> > matched_patterns;
  for(size_t i_p=0, tot_pat=fPatterns.size(); i_p<tot_pat; ++i_p){
    n_matched=0;

    for(auto const& i_b : fPatterns[i_p])
      if(fHitmap[i_b.first][i_b.second]>fIntegralThreshold)
	++n_matched;

    if( (float)n_matched/(float)(fPatterns[i_p].size()) > max_fraction_matched)
      max_fraction_matched = (float)n_matched/(float)(fPatterns[i_p].size());

    if(max_fraction_matched>fMatchFractionThreshold)
      break;
    
    //if(n_matched>MATCH_FRACTION*(fPatterns[i_p].size()))
    //matched_patterns.emplace_back(i_p,(float)n_matched / (float)(fPatterns[i_p].size()));
  }

  return max_fraction_matched;
}

void pm::AnodeCathodePMAlg::RunPatternMatching(std::vector<recob::Hit> const& hitVector,
					       float & frac_matched)
{
  FillHitmap(hitVector);
  frac_matched = GetFractionMatched();
}
