
////////////////////////////////////////////////////////////////////////
//
//    trivia : Class to hold the Muon Counter System (MuCS) data, 
//             a.k.a. the best class there's ever been, May 2015 
//    author : Odysseas Kanenas
//    e-mail : kalousis@vt.edu
//
////////////////////////////////////////////////////////////////////////

#include "MuCSData.h"

#include "cetlib/exception.h"

namespace MuCS {

MuCSData::MuCSData()
{}

MuCSData::~MuCSData()
{}
  
MuCSData::MuCSData( Float_t t0, Float_t adc1[24], Float_t adc2[24], Float_t adc3[24], Float_t adc7[24],
		    std::vector<Int_t> hits1, std::vector<Int_t> hits2, std::vector<Int_t> hits3, std::vector<Int_t> hits7 ) 
{
  ft0 = t0;
  
  for ( Int_t i=0; i<24; i++ ) 
    { 
      fadc1[i]=adc1[i]; fadc2[i]=adc2[i];
      fadc3[i]=adc3[i]; fadc7[i]=adc7[i];
      
    }
  
  Int_t s1 = hits1.size(); fhits1.clear();
  for ( Int_t i=0; i<s1; i++ ) fhits1.push_back( hits1.at(i) );
  
  Int_t s2 = hits2.size(); fhits2.clear();
  for ( Int_t i=0; i<s2; i++ ) fhits2.push_back( hits2.at(i) );
  
  Int_t s3 = hits3.size(); fhits3.clear();
  for ( Int_t i=0; i<s3; i++ ) fhits3.push_back( hits3.at(i) );
  
  Int_t s7 = hits7.size(); fhits7.clear();
  for ( Int_t i=0; i<s7; i++ ) fhits7.push_back( hits7.at(i) );
    
}
  
Float_t MuCSData::T0() const
{
  return ft0;
  
}

std::vector<Float_t> MuCSData::ADC1() const
{
  std::vector<Float_t> fadc; fadc.clear();
  for ( Int_t i=0; i<24; i++ ) fadc.push_back( fadc1[i] );
  return fadc;
  
}

std::vector<Float_t> MuCSData::ADC2() const
{
  std::vector<Float_t> fadc; fadc.clear();
  for ( Int_t i=0; i<24; i++ ) fadc.push_back( fadc2[i] );
  return fadc;
  
}

std::vector<Float_t> MuCSData::ADC3() const
{
  std::vector<Float_t> fadc; fadc.clear();
  for ( Int_t i=0; i<24; i++ ) fadc.push_back( fadc3[i] );
  return fadc;
  
}

std::vector<Float_t> MuCSData::ADC7() const
{
  std::vector<Float_t> fadc; fadc.clear();
  for ( Int_t i=0; i<24; i++ ) fadc.push_back( fadc7[i] );
  return fadc;
  
}

std::vector<Int_t> MuCSData::Hits1() const
{
  return fhits1;
  
}

std::vector<Int_t> MuCSData::Hits2() const
{
  return fhits2;
  
}

std::vector<Int_t> MuCSData::Hits3() const
{
  return fhits3;
  
}
  
std::vector<Int_t> MuCSData::Hits7() const
{
  return fhits7;
  
}

}

////////////////////////////////////////////////////////////////////////
//
//    The end !
//
////////////////////////////////////////////////////////////////////////
