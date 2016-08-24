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

#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/Hit.h"

#include <math.h>
#include <algorithm>
#include <iostream>

pm::AnodeCathodePMAlg::AnodeCathodePMAlg()
{
}

void pm::AnodeCathodePMAlg::Configure(fhicl::ParameterSet         const& p,
				      geo::GeometryCore           const& geo,
				      detinfo::DetectorProperties const& detp)
{

  //grab the wire information
  
  fPlaneID   = geo::PlaneID(p.get<unsigned int>("CryostatID",0),
			    p.get<unsigned int>("TPCID",0),
			    p.get<unsigned int>("PlaneID"));
  fStartWire = geo::WireID(fPlaneID,p.get<unsigned int>("StartWire",0));
  fEndWire = geo::WireID(fPlaneID,p.get<unsigned int>("EndWire",geo.Nwires(fPlaneID)-1));
  
  //check wire indices: end > start
  if(fEndWire.Wire < fStartWire.Wire){
    geo::WireID tmp(fEndWire);
    fEndWire = fStartWire;
    fStartWire = tmp;
  }

  if(fEndWire.Wire > geo.Nwires(fPlaneID)-1 )
    fEndWire = geo::WireID(fPlaneID,geo.Nwires(fPlaneID)-1);
  
  //fill nwires to run over. By default, this is all wires in the plane
  fNwires = fEndWire.Wire - fStartWire.Wire + 1;

  //set wires per bin, and set it right
  fWiresPerBin = p.get<unsigned int>("WiresPerBin");
  if(fNwires < fWiresPerBin)
    fWiresPerBin = fNwires;
  fNbinsWires = (fNwires / fWiresPerBin);
  
  //grab time tick info
  fStartTimeTick = p.get<float>("StartTime",0.0);
  fEndTimeTick   = p.get<float>("EndTime",(float)(detp.NumberTimeSamples()));

  //check wire indices: end > start
  if(fEndTimeTick < fStartTimeTick){
    float tmp = fEndTimeTick;
    fEndTimeTick = fStartTimeTick;
    fStartTimeTick = tmp;
  }
  
  if(fStartTimeTick < 0.0)
    fStartTimeTick = 0.0;
  if(fEndTimeTick > (float)(detp.NumberTimeSamples()) )
    fEndTimeTick = (float)(detp.NumberTimeSamples());
  
  fTimeTicksPerBin = p.get<float>("TimeTicksPerBin");  
  fNbinsTimeTicks = (unsigned int)((fEndTimeTick - fStartTimeTick)/fTimeTicksPerBin);

  
  fHitmap.resize(fNbinsWires,std::vector<float>(fNbinsTimeTicks));

  fNtimeTicksPattern = p.get<float>("NtimeTicksPattern");
  fNbinsPattern = (unsigned int)(std::round(fNtimeTicksPattern/fTimeTicksPerBin));
  
  fIntegralThreshold = p.get<float>("IntegralThreshold");
  fMatchFractionThreshold = p.get<float>("MatchFractionThreshold",10.0);
  
  fPatterns.reserve(fNbinsWires*fNbinsWires*fNbinsPattern);

  MakePatterns();
  
}

void pm::AnodeCathodePMAlg::ClearHitmap(){
  for(auto & w : fHitmap)
    std::fill(w.begin(),w.end(),0.0);
}

unsigned int pm::AnodeCathodePMAlg::GetBinTime(float time){
  if(time<fStartTimeTick || time>fEndTimeTick) {
    std::cout << "ERROR!!! " << time << " not inside bounds [ "
	      << fStartTimeTick << "," << fEndTimeTick << "]" << std::endl;
    sleep(1);
    return -1;
  }
  return (time - fStartTimeTick)/fTimeTicksPerBin;
}

unsigned int pm::AnodeCathodePMAlg::GetBinWire(float wire){
  if(wire<fStartWire.Wire || wire>fEndWire.Wire) {
    std::cout << "ERROR!!! " << wire << " not inside bounds [ "
	      << fStartWire.Wire << "," << fEndWire.Wire << "]" << std::endl;
    sleep(1);
    return -1;
  }
  return (wire - fStartWire.Wire)/fWiresPerBin;
}

void pm::AnodeCathodePMAlg::MakePatterns(){

  float wire_a,wire_c;
  float slope, intercept, time, time_start, time_end;
  for(size_t i_p=0; i_p<fNbinsWires*fNbinsWires; ++i_p){

    wire_a = (float)(fStartWire.Wire + (i_p/fNbinsWires)*fWiresPerBin + fWiresPerBin/2);
    wire_c = (float)(fStartWire.Wire + (i_p%fNbinsWires)*fWiresPerBin + fWiresPerBin/2);

    for(time_start=fStartTimeTick, time_end=fStartTimeTick+fNtimeTicksPattern;
	time_end<=fEndTimeTick;
	time_start+=fTimeTicksPerBin, time_end+= fTimeTicksPerBin){
      slope = (wire_c - wire_a)/(float)(time_end-time_start);
      intercept = (time_start)*slope;

      fPatterns.push_back(Pattern_t());
      for(time=time_start; time<=time_end; time += 0.1*fTimeTicksPerBin){
	
	fPatterns.back().emplace(GetBinWire(wire_a+time*slope-intercept),GetBinTime(time));
      }
    }
  }

}

void pm::AnodeCathodePMAlg::FillHitmap( std::vector<recob::Hit> const& hits)
{
  ClearHitmap();
  for(auto const& hit : hits){
    if(hit.WireID().planeID()!=fPlaneID) continue;
    if(hit.PeakTime() > fEndTimeTick || hit.PeakTime() < fStartTimeTick) continue;
    if(hit.WireID().Wire < fStartWire.Wire || hit.WireID().Wire > fEndWire.Wire) continue;
    fHitmap[GetBinWire(hit.WireID())][GetBinTime(hit.PeakTime())] += hit.Integral();
  }
}

float pm::AnodeCathodePMAlg::GetFractionMatched()
{
  unsigned int n_matched=0;
  float max_fraction_matched = -1.0;

  for(size_t i_p=0, tot_pat=fPatterns.size(); i_p<tot_pat; ++i_p){
    n_matched=0;

    for(auto const& i_b : fPatterns[i_p])
      if(fHitmap[i_b.first][i_b.second]>fIntegralThreshold)
	++n_matched;

    if( (float)n_matched/(float)(fPatterns[i_p].size()) > max_fraction_matched)
      max_fraction_matched = (float)n_matched/(float)(fPatterns[i_p].size());
    
    if(max_fraction_matched>fMatchFractionThreshold)
      break;    
  }

  return max_fraction_matched;
}

void pm::AnodeCathodePMAlg::RunPatternMatching(std::vector<recob::Hit> const& hitVector,
					       float & frac_matched)
{
  FillHitmap(hitVector);
  frac_matched = GetFractionMatched();
}

