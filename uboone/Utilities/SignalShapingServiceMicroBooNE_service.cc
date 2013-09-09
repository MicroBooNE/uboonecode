////////////////////////////////////////////////////////////////////////
/// \file   SignalShapingServiceMicroBooNE_service.cc
/// \author H. Greenlee 
////////////////////////////////////////////////////////////////////////

#include "Utilities/SignalShapingServiceMicroBooNE.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Utilities/LArFFT.h"
#include "TFile.h"

//----------------------------------------------------------------------
// Constructor.
util::SignalShapingServiceMicroBooNE::SignalShapingServiceMicroBooNE(const fhicl::ParameterSet& pset,
								    art::ActivityRegistry& reg) 
  : fInit(false)
{
  reconfigure(pset);
}


//----------------------------------------------------------------------
// Destructor.
util::SignalShapingServiceMicroBooNE::~SignalShapingServiceMicroBooNE()
{}


//----------------------------------------------------------------------
// Reconfigure method.
void util::SignalShapingServiceMicroBooNE::reconfigure(const fhicl::ParameterSet& pset)
{
  // Reset initialization flag.

  fInit = false;

  // Reset kernels.

  fColSignalShaping.Reset();
  fIndSignalShaping.Reset();

  // Fetch fcl parameters.

  fNFieldBins = pset.get<int>("FieldBins");
  fCol3DCorrection = pset.get<double>("Col3DCorrection");
  fInd3DCorrection = pset.get<double>("Ind3DCorrection");
  fColFieldRespAmp = pset.get<double>("ColFieldRespAmp");
  fIndFieldRespAmp = pset.get<double>("IndFieldRespAmp");
  fShapeTimeConst = pset.get<std::vector<double> >("ShapeTimeConst");

  fUseFunctionFieldShape= pset.get<bool>("UseFunctionFieldShape");
  fGetFilterFromHisto= pset.get<bool>("GetFilterFromHisto");
  
  // Construct parameterized collection filter function.
 if(!fGetFilterFromHisto)
 {
  mf::LogInfo("SignalShapingServiceMicroBooNE") << "Getting Filter from .fcl file" ;
  std::string colFilt = pset.get<std::string>("ColFilter");
  std::vector<double> colFiltParams =
  pset.get<std::vector<double> >("ColFilterParams");
  fColFilterFunc = new TF1("colFilter", colFilt.c_str());
  for(unsigned int i=0; i<colFiltParams.size(); ++i)
    fColFilterFunc->SetParameter(i, colFiltParams[i]);

  // Construct parameterized induction filter function.

  std::string indFilt = pset.get<std::string>("IndFilter");
  std::vector<double> indFiltParams =
  pset.get<std::vector<double> >("IndFilterParams");
  fIndFilterFunc = new TF1("indFilter", indFilt.c_str());
  for(unsigned int i=0; i<indFiltParams.size(); ++i)
    fIndFilterFunc->SetParameter(i, indFiltParams[i]);
 }
 else
 {
  
   std::string histoname = pset.get<std::string>("FilterHistoName");
   mf::LogInfo("SignalShapingServiceMicroBooNE") << " using filter from .root file " ;
   int fNPlanes=3;
   
  // constructor decides if initialized value is a path or an environment variable
  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(pset.get<std::string>("FilterFunctionFname"), fname);
    
  TFile * in=new TFile(fname.c_str(),"READ");
   for(int i=0;i<fNPlanes;i++){
     TH1D * temp=(TH1D *)in->Get(Form(histoname.c_str(),i));
     fFilterHist[i]=new TH1D(Form(histoname.c_str(),i),Form(histoname.c_str(),i),temp->GetNbinsX(),0,temp->GetNbinsX());
     temp->Copy(*fFilterHist[i]); 
    }
   
   in->Close();
   
 }
 
 /////////////////////////////////////
 if(fUseFunctionFieldShape)
 {
  std::string colField = pset.get<std::string>("ColFieldShape");
  std::vector<double> colFieldParams =
    pset.get<std::vector<double> >("ColFieldParams");
  fColFieldFunc = new TF1("colField", colField.c_str());
  for(unsigned int i=0; i<colFieldParams.size(); ++i)
    fColFieldFunc->SetParameter(i, colFieldParams[i]);

  // Construct parameterized induction filter function.

  std::string indField = pset.get<std::string>("IndFieldShape");
  std::vector<double> indFieldParams =
    pset.get<std::vector<double> >("IndFieldParams");
  fIndFieldFunc = new TF1("indField", indField.c_str());
  for(unsigned int i=0; i<indFieldParams.size(); ++i)
    fIndFieldFunc->SetParameter(i, indFieldParams[i]);
 }   // Warning, last parameter needs to be multiplied by the FFTSize, in current version of the code,
 
 
}


//----------------------------------------------------------------------
// Accessor for single-plane signal shaper.
const util::SignalShaping&
util::SignalShapingServiceMicroBooNE::SignalShaping(unsigned int channel) const
{
  if(!fInit)
    init();

  // Figure out plane type.

  art::ServiceHandle<geo::Geometry> geom;
  geo::SigType_t sigtype = geom->SignalType(channel);

  // Return appropriate shaper.

  if(sigtype == geo::kInduction)
    return fIndSignalShaping;
  else if(sigtype == geo::kCollection)
    return fColSignalShaping;
  else
    throw cet::exception("SignalShapingServiceMicroBooNE")<< "can't determine"
                                                          << " SignalType";
							  
return fColSignalShaping;
}


//----------------------------------------------------------------------
// Initialization method.
// Here we do initialization that can't be done in the constructor.
// All public methods should ensure that this method is called as necessary.
void util::SignalShapingServiceMicroBooNE::init()
{
  if(!fInit) {
    fInit = true;

    // Do microboone-specific configuration of SignalShaping by providing
    // microboone response and filter functions.

    // Calculate field and electronics response functions.

    SetFieldResponse();
    SetElectResponse();

    // Configure convolution kernels.

    fColSignalShaping.AddResponseFunction(fColFieldResponse);
    fColSignalShaping.AddResponseFunction(fElectResponse);
    fColSignalShaping.SetPeakResponseTime(0.);

    fIndSignalShaping.AddResponseFunction(fIndFieldResponse);
    fIndSignalShaping.AddResponseFunction(fElectResponse);
    fIndSignalShaping.SetPeakResponseTime(0.);

    // Calculate filter functions.

    SetFilters();

    // Configure deconvolution kernels.

    fColSignalShaping.AddFilterFunction(fColFilter);
    fColSignalShaping.CalculateDeconvKernel();

    fIndSignalShaping.AddFilterFunction(fIndFilter);
    fIndSignalShaping.CalculateDeconvKernel();
  }
}


//----------------------------------------------------------------------
// Calculate microboone field response.
void util::SignalShapingServiceMicroBooNE::SetFieldResponse()
{
  // Get services.

  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArProperties> larp;

  // Get plane pitch.
 
  double xyz1[3] = {0.};
  double xyz2[3] = {0.};
  double xyzl[3] = {0.};
  // should always have at least 2 planes
  geo->Plane(0).LocalToWorld(xyzl, xyz1);
  geo->Plane(1).LocalToWorld(xyzl, xyz2);

  // this assumes all planes are equidistant from each other,
  // probably not a bad assumption
  double pitch = xyz2[0] - xyz1[0]; ///in cm

  fColFieldResponse.resize(fNFieldBins, 0.);
  fIndFieldResponse.resize(fNFieldBins, 0.);

  // set the response for the collection plane first
  // the first entry is 0

  double driftvelocity=larp->DriftVelocity()/1000.;  
  int nbinc = TMath::Nint(fCol3DCorrection*(std::abs(pitch))/(driftvelocity*detprop->SamplingRate())); ///number of bins //KP
  
  ////////////////////////////////////////////////////
   if(fUseFunctionFieldShape)
  {
  art::ServiceHandle<util::LArFFT> fft;
  int signalSize = fft->FFTSize();
  std::vector<double> ramp(signalSize);
   // TComplex kernBin;
   // int size = signalSize/2;
   // int bin=0;
    //std::vector<TComplex> freqSig(size+1);
  std::vector<double> bipolar(signalSize);
    
    
  fColFieldResponse.resize(signalSize, 0.);
  fIndFieldResponse.resize(signalSize, 0.);
   
  // Hardcoding. Bad. Temporary hopefully.
  fIndFieldFunc->SetParameter(4,fIndFieldFunc->GetParameter(4)*signalSize);
  
  
  double integral = 0.;
    for(int i = 0; i < signalSize; i++) {
          ramp[i]=fColFieldFunc->Eval(i);
          fColFieldResponse[i]=ramp[i];
          integral += fColFieldResponse[i];
     // rampc->Fill(i,ramp[i]);
      bipolar[i]=fIndFieldFunc->Eval(i);
      fIndFieldResponse[i]=bipolar[i];
     // bipol->Fill(i,bipolar[i]);
    }
     
   for(int i = 0; i < signalSize; ++i){
      fColFieldResponse[i] *= fColFieldRespAmp/integral;
     }
      
    //this might be not necessary if the function definition is not defined in the middle of the signal range  
    fft->ShiftData(fIndFieldResponse,signalSize/2.0);
  }
  else
  {
  //////////////////////////////////////////////////
  mf::LogInfo("SignalShapingServiceMicroBooNE") << " using the old field shape " ;
  double integral = 0.;
  for(int i = 1; i < nbinc; ++i){
    fColFieldResponse[i] = fColFieldResponse[i-1] + 1.0;
    integral += fColFieldResponse[i];
  }

  for(int i = 0; i < nbinc; ++i){
    fColFieldResponse[i] *= fColFieldRespAmp/integral;
  }

  // now the induction plane
    
  int nbini = TMath::Nint(fInd3DCorrection*(std::abs(pitch))/(driftvelocity*detprop->SamplingRate()));//KP
  for(int i = 0; i < nbini; ++i){
    fIndFieldResponse[i] = fIndFieldRespAmp/(1.*nbini);
    fIndFieldResponse[nbini+i] = -fIndFieldRespAmp/(1.*nbini);
    }

  }
  
  return;
}


//----------------------------------------------------------------------
// Calculate microboone field response.
void util::SignalShapingServiceMicroBooNE::SetElectResponse()
{
  // Get services.

  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArFFT> fft;

  LOG_DEBUG("SignalShapingMicroBooNE") << "Setting MicroBooNE electronics response function...";

  int nticks = fft->FFTSize();
  fElectResponse.resize(nticks, 0.);
  std::vector<double> time(nticks,0.);

  //Gain and shaping time variables from fcl file:    
  double Ao = fShapeTimeConst[0];  //gain
  double To = fShapeTimeConst[1];  //peaking time
    
  // this is actually sampling time, in ns
  // mf::LogInfo("SignalShapingMicroBooNE") << "Check sampling intervals: " 
  //                                  << fSampleRate << " ns" 
  //                                  << "Check number of samples: " << fNTicks;

  // The following sets the microboone electronics response function in 
  // time-space. Function comes from BNL SPICE simulation of MicroBooNE 
  // electronics. SPICE gives the electronics transfer function in 
  // frequency-space. The inverse laplace transform of that function 
  // (in time-space) was calculated in Mathematica and is what is being 
  // used below. Parameters Ao and To are cumulative gain/timing parameters 
  // from the full (ASIC->Intermediate amp->Receiver->ADC) electronics chain. 
  // They have been adjusted to make the SPICE simulation to match the 
  // actual electronics response. Default params are Ao=1.4, To=0.5us. 
  double integral=0.;
  
  for(int i = 0; i < nticks; ++i){

    //convert time to microseconds, to match fElectResponse[i] definition
    time[i] = (1.*i)*detprop->SamplingRate()*1e-3; 
    fElectResponse[i] = 
      4.31054*exp(-2.94809*time[i]/To)*Ao - 2.6202*exp(-2.82833*time[i]/To)*cos(1.19361*time[i]/To)*Ao
      -2.6202*exp(-2.82833*time[i]/To)*cos(1.19361*time[i]/To)*cos(2.38722*time[i]/To)*Ao
      +0.464924*exp(-2.40318*time[i]/To)*cos(2.5928*time[i]/To)*Ao
      +0.464924*exp(-2.40318*time[i]/To)*cos(2.5928*time[i]/To)*cos(5.18561*time[i]/To)*Ao
      +0.762456*exp(-2.82833*time[i]/To)*sin(1.19361*time[i]/To)*Ao
      -0.762456*exp(-2.82833*time[i]/To)*cos(2.38722*time[i]/To)*sin(1.19361*time[i]/To)*Ao
      +0.762456*exp(-2.82833*time[i]/To)*cos(1.19361*time[i]/To)*sin(2.38722*time[i]/To)*Ao
      -2.6202*exp(-2.82833*time[i]/To)*sin(1.19361*time[i]/To)*sin(2.38722*time[i]/To)*Ao 
      -0.327684*exp(-2.40318*time[i]/To)*sin(2.5928*time[i]/To)*Ao + 
      +0.327684*exp(-2.40318*time[i]/To)*cos(5.18561*time[i]/To)*sin(2.5928*time[i]/To)*Ao
      -0.327684*exp(-2.40318*time[i]/To)*cos(2.5928*time[i]/To)*sin(5.18561*time[i]/To)*Ao
      +0.464924*exp(-2.40318*time[i]/To)*sin(2.5928*time[i]/To)*sin(5.18561*time[i]/To)*Ao;

      integral+=fElectResponse[i];
  }// end loop over time buckets
    
  LOG_DEBUG("SignalShapingMicroBooNE") << " Done.";

 //normalize fElectResponse[i], before the convolution   
   for(int i = 0; i < nticks; ++i){
     fElectResponse[i]/=integral;
   }
  
  
  return;

}


//----------------------------------------------------------------------
// Calculate microboone filter functions.
void util::SignalShapingServiceMicroBooNE::SetFilters()
{
  // Get services.

  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArFFT> fft;

  double ts = detprop->SamplingRate();
  int n = fft->FFTSize() / 2;

  // Calculate collection filter.

  fColFilter.resize(n+1);
  fIndFilter.resize(n+1);
  
  if(!fGetFilterFromHisto)
  {
  fColFilterFunc->SetRange(0, double(n));

  for(int i=0; i<=n; ++i) {
    double freq = 500. * i / (ts * n);      // Cycles / microsecond.
    double f = fColFilterFunc->Eval(freq);
    fColFilter[i] = TComplex(f, 0.);
  }
  

  // Calculate induction filter.

 
  fIndFilterFunc->SetRange(0, double(n));

  for(int i=0; i<=n; ++i) {
    double freq = 500. * i / (ts * n);      // Cycles / microsecond.
    double f = fIndFilterFunc->Eval(freq);
    fIndFilter[i] = TComplex(f, 0.);
    }
  
  }
  else
  {
    
    for(int i=0; i<=n; ++i) {
      double f = fFilterHist[2]->GetBinContent(i);  // hardcoded plane numbers. Bad. To change later.
      fColFilter[i] = TComplex(f, 0.);
      double g = fFilterHist[0]->GetBinContent(i);
      fIndFilter[i] = TComplex(g, 0.);
      
    }
  }
  
  fIndSignalShaping.AddFilterFunction(fIndFilter);
  fColSignalShaping.AddFilterFunction(fColFilter);
  
}



namespace util {

  DEFINE_ART_SERVICE(SignalShapingServiceMicroBooNE)

}
