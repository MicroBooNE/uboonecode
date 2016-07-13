#ifndef OVERLAY_DATAOVERLAY_RAWDIGITADDER_CXX
#define OVERLAY_DATAOVERLAY_RAWDIGITADDER_CXX

#include "RawDigitAdder.h"
#include <limits>
#include <stdexcept>

void mix::RawDigitAdder::AddRawDigits( std::vector<short> const& inputVec1,
				  std::vector<short> const& inputVec2,
				  std::vector<short>& outputVec)
{
  CheckVectorSize(inputVec1,inputVec2);
  outputVec.clear(); outputVec.resize(inputVec1.size());
  for(size_t i=0; i<inputVec1.size(); i++)
    AddRawDigit(inputVec1[i]-(short)_ped1,inputVec2[i]-(short)_ped2,outputVec[i]);
}

void mix::RawDigitAdder::AddRawDigits( std::vector<short> const& inputVec1,
				  std::vector<short>& inputVec2)
{
  CheckVectorSize(inputVec1,inputVec2);
  for(size_t i=0; i<inputVec1.size(); i++)
    AddRawDigit(inputVec1[i]-(short)_ped1,inputVec2[i]-(short)_ped2,inputVec2[i]);
}

void mix::RawDigitAdder::AddRawDigits( std::vector< std::vector<short> > const& inputVecList,
				  std::vector<short>& outputVec)
{
  if(inputVecList.size()==0){
    if(_throw)
      throw std::runtime_error("Error in RawDigitAdder::AddRawDigits : Input vector list is empty.");
    return;
  }
  for(auto const& vec : inputVecList)
    CheckVectorSize(inputVecList[0],vec);

  outputVec.clear(); outputVec.resize(inputVecList[0].size(),0.0);

  for(auto const& inputVec : inputVecList)
    AddRawDigits(inputVec,outputVec);

}

void mix::RawDigitAdder::CheckVectorSize(std::vector<short> const& vec1,
				    std::vector<short> const& vec2)
{
  if(vec1.size() != vec2.size()){
    if(_throw)
      throw std::runtime_error("Error in RawDigitAdder: Failed CheckVectorSize()");
    return;
  }
}

void mix::RawDigitAdder::AddRawDigit(short const& d1, short const& d2, short& d_out)
{
  d_out = d1 + d2;
  FixOverflow(d_out);
}

void mix::RawDigitAdder::AddRawDigit(short const& d1, short& d2)
{
  AddRawDigit(d1,d2,d2);
}

void mix::RawDigitAdder::FixOverflow(short& d)
{ 
  if(d<0)
    d = std::numeric_limits<short>::max();
}


void mix::RawDigitAdder::SetPedestalInput(float f, float& _ped)
{
  _ped = f;
}

#endif
