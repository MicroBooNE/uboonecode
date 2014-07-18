
////////////////////////////////////////////////////////////////////////
/// \file   SignalShapingServiceMicroBooNE_service.cc
/// \author H. Greenlee 
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
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
								    art::ActivityRegistry& /* reg */) 
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
  fIndUSignalShaping.Reset();
  fIndVSignalShaping.Reset();

  // Fetch fcl parameters.

  fADCTicksPerPCAtLowestASICGainSetting = pset.get<double>("ADCTicksPerPCAtLowestASICGainSetting");
  fASICGainInMVPerFC = pset.get<double>("ASICGainInMVPerFC");
  fDefaultDriftVelocity = pset.get<double>("DefaultDriftVelocity");
  fNFieldBins = pset.get<int>("FieldBins");
  fInputFieldRespSamplingRate = pset.get<double>("InputFieldRespSamplingRate");
  fCol3DCorrection = pset.get<double>("Col3DCorrection");
  fInd3DCorrection = pset.get<double>("Ind3DCorrection");
  fColFieldRespAmp = pset.get<double>("ColFieldRespAmp");
  fIndUFieldRespAmp = pset.get<double>("IndUFieldRespAmp");
  fIndVFieldRespAmp = pset.get<double>("IndVFieldRespAmp");
  fShapeTimeConst = pset.get<std::vector<double> >("ShapeTimeConst");

  // Currently we have three options of field response shapes:
  // 1. UseFunctionFieldShape
  // 2. UseHistogramFieldShape, i.e. the waveforms from Leon Rochester
  // 3. Otherwise, a square function for induction planes, and
  //    a ramp function for collection planes
  fUseFunctionFieldShape= pset.get<bool>("UseFunctionFieldShape");
  fUseHistogramFieldShape = pset.get<bool>("UseHistogramFieldShape");

  fGetFilterFromHisto= pset.get<bool>("GetFilterFromHisto");

  // Construct parameterized collection filter function.
  if(!fGetFilterFromHisto) {

    mf::LogInfo("SignalShapingServiceMicroBooNE") << "Getting Filter from .fcl file" ;
    std::string colFilt = pset.get<std::string>("ColFilter");
    std::vector<double> colFiltParams = pset.get<std::vector<double> >("ColFilterParams");
    fColFilterFunc = new TF1("colFilter", colFilt.c_str());
    for(unsigned int i=0; i<colFiltParams.size(); ++i)
      fColFilterFunc->SetParameter(i, colFiltParams[i]);
    
    // Construct parameterized induction filter function.

    std::string indUFilt = pset.get<std::string>("IndUFilter");
    std::vector<double> indUFiltParams = pset.get<std::vector<double> >("IndUFilterParams");
    fIndUFilterFunc = new TF1("indUFilter", indUFilt.c_str());
    for(unsigned int i=0; i<indUFiltParams.size(); ++i)
      fIndUFilterFunc->SetParameter(i, indUFiltParams[i]);

    std::string indVFilt = pset.get<std::string>("IndVFilter");
    std::vector<double> indVFiltParams = pset.get<std::vector<double> >("IndVFilterParams");
    fIndVFilterFunc = new TF1("indVFilter", indVFilt.c_str());
    for(unsigned int i=0; i<indVFiltParams.size(); ++i)
      fIndVFilterFunc->SetParameter(i, indVFiltParams[i]);
  } else {
  
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
    delete in; 
  }
 
  /////////////////////////////////////
  if(fUseFunctionFieldShape) {

    std::string colField = pset.get<std::string>("ColFieldShape");
    std::vector<double> colFieldParams = pset.get<std::vector<double> >("ColFieldParams");
    fColFieldFunc = new TF1("colField", colField.c_str());
    for(unsigned int i=0; i<colFieldParams.size(); ++i)
      fColFieldFunc->SetParameter(i, colFieldParams[i]);

    // Construct parameterized induction filter function.
    
    std::string indUField = pset.get<std::string>("IndUFieldShape");
    std::vector<double> indUFieldParams = pset.get<std::vector<double> >("IndUFieldParams");
    fIndUFieldFunc = new TF1("indUField", indUField.c_str());
    for(unsigned int i=0; i<indUFieldParams.size(); ++i)
      fIndUFieldFunc->SetParameter(i, indUFieldParams[i]);
    // Warning, last parameter needs to be multiplied by the FFTSize, in current version of the code,
    
    std::string indVField = pset.get<std::string>("IndVFieldShape");
    std::vector<double> indVFieldParams = pset.get<std::vector<double> >("IndVFieldParams");
    fIndVFieldFunc = new TF1("indVField", indVField.c_str());
    for(unsigned int i=0; i<indVFieldParams.size(); ++i)
      fIndVFieldFunc->SetParameter(i, indVFieldParams[i]);
    // Warning, last parameter needs to be multiplied by the FFTSize, in current version of the code,

  } else if ( fUseHistogramFieldShape ) {
    mf::LogInfo("SignalShapingServiceMicroBooNE") << " using the field response provided from a .root file " ;
    int fNPlanes = 3;

    // constructor decides if initialized value is a path or an environment variable
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file( pset.get<std::string>("FieldResponseFname"), fname );
    std::string histoname = pset.get<std::string>("FieldResponseHistoName");

    std::unique_ptr<TFile> fin(new TFile(fname.c_str(), "READ"));
    if ( !fin->IsOpen() ) throw art::Exception( art::errors::NotFound ) << "Could not find the field response file " << fname << "!" << std::endl;

    std::string iPlane[3] = { "U", "V", "Y" };

    for ( int i = 0; i < fNPlanes; i++ ) {
      TString iHistoName = Form( "%s_%s", histoname.c_str(), iPlane[i].c_str());
      TH1F *temp = (TH1F*) fin->Get( iHistoName );
      if ( !temp ) throw art::Exception( art::errors::NotFound ) << "Could not find the field response histogram " << iHistoName << std::endl;

      if ( temp->GetNbinsX() > fNFieldBins ) throw art::Exception( art::errors::InvalidNumber ) << "FieldBins should always be larger than or equal to the number of the bins in the input histogram!" << std::endl;

      fFieldResponseHist[i] = new TH1F( iHistoName, iHistoName, temp->GetNbinsX(), temp->GetBinLowEdge(1), temp->GetBinLowEdge( temp->GetNbinsX() + 1) );
      temp->Copy(*fFieldResponseHist[i]);
    }

    fin->Close();
  } else {
    fNFieldBins = 300;
  }  
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
  //geo::SigType_t sigtype = geom->SignalType(channel);

  // we need to distiguish the U and V planes
  geo::View_t view = geom->View(channel);

  // Return appropriate shaper.

  if(view == geo::kU)
    return fIndUSignalShaping;
  else if(view == geo::kV)
    return fIndVSignalShaping;
  else if(view == geo::kZ)
    return fColSignalShaping;
  else
    throw cet::exception("SignalShapingServiceMicroBooNE")<< "can't determine"
                                                          << " View\n";
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
    // fColSignalShaping.SetPeakResponseTime(0.);

    fIndUSignalShaping.AddResponseFunction(fIndUFieldResponse);
    fIndUSignalShaping.AddResponseFunction(fElectResponse);
    // fIndUSignalShaping.SetPeakResponseTime(0.);

    fIndVSignalShaping.AddResponseFunction(fIndVFieldResponse);
    fIndVSignalShaping.AddResponseFunction(fElectResponse);
    // fIndVSignalShaping.SetPeakResponseTime(0.);

    // Currently we only have fine binning "fInputFieldRespSamplingRate"
    // for the field and electronic responses.
    // Now we are sampling the convoluted field-electronic response
    // with the nominal sampling rate.
    // We may consider to do the same for the filters as well.
    SetResponseSampling();

    // Calculate filter functions.

    SetFilters();

    // Configure deconvolution kernels.

    fColSignalShaping.AddFilterFunction(fColFilter);
    fColSignalShaping.CalculateDeconvKernel();

    fIndUSignalShaping.AddFilterFunction(fIndUFilter);
    fIndUSignalShaping.CalculateDeconvKernel();

    fIndVSignalShaping.AddFilterFunction(fIndVFilter);
    fIndVSignalShaping.CalculateDeconvKernel();
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
  fIndUFieldResponse.resize(fNFieldBins, 0.);
  fIndVFieldResponse.resize(fNFieldBins, 0.);

  // set the response for the collection plane first
  // the first entry is 0

  double driftvelocity=larp->DriftVelocity()/1000.; // in cm/nsec 
  double integral = 0.;
  ////////////////////////////////////////////////////
  if(fUseFunctionFieldShape) {

    art::ServiceHandle<util::LArFFT> fft;
    int signalSize = fft->FFTSize();
    std::vector<double> ramp(signalSize);
    // TComplex kernBin;
    // int size = signalSize/2;
    // int bin=0;
    //std::vector<TComplex> freqSig(size+1);
    std::vector<double> bipolar(signalSize);    
    
    fColFieldResponse.resize(signalSize, 0.);
    fIndUFieldResponse.resize(signalSize, 0.);
    fIndVFieldResponse.resize(signalSize, 0.);
   
    // Hardcoding. Bad. Temporary hopefully.
    fIndUFieldFunc->SetParameter(4,fIndUFieldFunc->GetParameter(4)*signalSize);
    fIndVFieldFunc->SetParameter(4,fIndVFieldFunc->GetParameter(4)*signalSize);

    for(int i = 0; i < signalSize; i++) {
      ramp[i]=fColFieldFunc->Eval(i);
      fColFieldResponse[i]=ramp[i];
      integral += fColFieldResponse[i];
      // rampc->Fill(i,ramp[i]);
      bipolar[i]=fIndUFieldFunc->Eval(i);
      fIndUFieldResponse[i]=bipolar[i];
      bipolar[i]=fIndVFieldFunc->Eval(i);
      fIndVFieldResponse[i]=bipolar[i];
      // bipol->Fill(i,bipolar[i]);
    }
     
    for(int i = 0; i < signalSize; ++i){
      fColFieldResponse[i] *= fColFieldRespAmp/integral;
    }
      
    //this might be not necessary if the function definition is not defined in the middle of the signal range  
    fft->ShiftData(fIndUFieldResponse,signalSize/2.0);
    fft->ShiftData(fIndVFieldResponse,signalSize/2.0);

  } else if ( fUseHistogramFieldShape ) {
    // Ticks in nanosecond
    // Calculate the normalization of the collection plane
    for ( int ibin = 1; ibin <= fFieldResponseHist[2]->GetNbinsX(); ibin++ )
      integral += fFieldResponseHist[2]->GetBinContent( ibin ) * fInputFieldRespSamplingRate;

    // Induction plane
    for ( int ibin = 1; ibin <= fFieldResponseHist[0]->GetNbinsX(); ibin++ )
      fIndUFieldResponse[ibin-1] = fIndUFieldRespAmp*fFieldResponseHist[0]->GetBinContent( ibin )/integral;

    for ( int ibin = 1; ibin <= fFieldResponseHist[1]->GetNbinsX(); ibin++ )
      fIndVFieldResponse[ibin-1] = fIndVFieldRespAmp*fFieldResponseHist[1]->GetBinContent( ibin )/integral;

    for ( int ibin = 1; ibin <= fFieldResponseHist[2]->GetNbinsX(); ibin++ )
      fColFieldResponse[ibin-1] = fColFieldRespAmp*fFieldResponseHist[2]->GetBinContent( ibin )/integral;

  } else {

    //////////////////////////////////////////////////
    mf::LogInfo("SignalShapingServiceMicroBooNE") << " using the old field shape " ;
    int nbinc = TMath::Nint(fCol3DCorrection*(std::abs(pitch))/(driftvelocity* fInputFieldRespSamplingRate)); ///number of bins //KP
    
    double integral = 0.;
    for(int i = 1; i < nbinc; ++i){
      fColFieldResponse[i] = fColFieldResponse[i-1] + 1.0;
      integral += fColFieldResponse[i];
    }

    for(int i = 0; i < nbinc; ++i){
      fColFieldResponse[i] *= fColFieldRespAmp/integral;
    }

    // now the induction plane
    
    int nbini = TMath::Nint(fInd3DCorrection*(std::abs(pitch))/(driftvelocity*fInputFieldRespSamplingRate));//KP
    for(int i = 0; i < nbini; ++i){
      fIndUFieldResponse[i] = fIndUFieldRespAmp/(1.*nbini);
      fIndUFieldResponse[nbini+i] = -fIndUFieldRespAmp/(1.*nbini);
    }

    for(int i = 0; i < nbini; ++i){
      fIndVFieldResponse[i] = fIndVFieldRespAmp/(1.*nbini);
      fIndVFieldResponse[nbini+i] = -fIndVFieldRespAmp/(1.*nbini);
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
    time[i] = (1.*i) * fInputFieldRespSamplingRate * 1e-3; 
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

  // normalize fElectResponse[i], before the convolution
  // Put in overall normalization in a pedantic way:
  // first put in the pulse area per eleectron at the lowest gain setting,
  // then normalize by the actual ASIC gain setting used.
  // This code is executed only during initialization of service,
  // so don't worry about code inefficiencies here.
   for(int i = 0; i < nticks; ++i){
     // The sampling rate was not taken into account in the older version
     fElectResponse[i] /= integral * fInputFieldRespSamplingRate;
     fElectResponse[i] *= fADCTicksPerPCAtLowestASICGainSetting*1.6e-7;
     fElectResponse[i] *= fASICGainInMVPerFC/4.7;
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
  fIndUFilter.resize(n+1);
  fIndVFilter.resize(n+1);
  
  if(!fGetFilterFromHisto)
  {
  fColFilterFunc->SetRange(0, double(n));

  for(int i=0; i<=n; ++i) {
    double freq = 500. * i / (ts * n);      // Cycles / microsecond.
    double f = fColFilterFunc->Eval(freq);
    fColFilter[i] = TComplex(f, 0.);
  }
  

  // Calculate induction filters.
 
  fIndUFilterFunc->SetRange(0, double(n));

  for(int i=0; i<=n; ++i) {
    double freq = 500. * i / (ts * n);      // Cycles / microsecond.
    double f = fIndUFilterFunc->Eval(freq);
    fIndUFilter[i] = TComplex(f, 0.);
    }

   fIndVFilterFunc->SetRange(0, double(n));

  for(int i=0; i<=n; ++i) {
    double freq = 500. * i / (ts * n);      // Cycles / microsecond.
    double f = fIndVFilterFunc->Eval(freq);
    fIndVFilter[i] = TComplex(f, 0.);
    }
  
  }
  else
  {
    
    for(int i=0; i<=n; ++i) {
      double f = fFilterHist[2]->GetBinContent(i);  // hardcoded plane numbers. Bad. To change later.
      fColFilter[i] = TComplex(f, 0.);
      double g = fFilterHist[1]->GetBinContent(i);
      fIndVFilter[i] = TComplex(g, 0.);
      double h = fFilterHist[0]->GetBinContent(i);
      fIndUFilter[i] = TComplex(h, 0.);
      
    }
  }
  
  fIndUSignalShaping.AddFilterFunction(fIndUFilter);
  fIndVSignalShaping.AddFilterFunction(fIndVFilter);
  fColSignalShaping.AddFilterFunction(fColFilter);
  
}

//----------------------------------------------------------------------
// Sample microboone response (the convoluted field and electronic
// response), will probably add the filter later
void util::SignalShapingServiceMicroBooNE::SetResponseSampling()
{
  // Get services
  art::ServiceHandle<util::LArProperties> larp;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArFFT> fft;

  // We allow different drift velocities
  double kDVel = larp->DriftVelocity() / fDefaultDriftVelocity;
  // Rounding
  if ( std::abs( 1. - kDVel ) < 10e-5 ) kDVel = 1.;
  if ( kDVel > 1. ) mf::LogInfo("SignalShapingServiceMicroBooNE") << "The drift velocity "
    << larp->DriftVelocity() << "cm/usec is faster than the default " << fDefaultDriftVelocity
    << "cm/usec!" << std::endl;

  int nticks = fft->FFTSize();
  std::vector<double> InputTime( nticks, 0. );
  std::vector<double> SamplingTime( nticks, 0. );

  for ( int itime = 0; itime < nticks; itime++ ) {
    InputTime[itime] = (1.*itime) * fInputFieldRespSamplingRate;
    SamplingTime[itime] = (1.*itime) * detprop->SamplingRate() / kDVel;
  }

  // Sampling
  int fNPlanes = 3;
  for ( int iplane = 0; iplane < fNPlanes; iplane++ ) {
    const std::vector<double>* pResp;
    std::vector<double> SamplingResp( nticks, 0. );
    int SamplingCount = 0;

    switch ( iplane ) {
      case 0: pResp = &(fIndUSignalShaping.Response()); break;
      case 1: pResp = &(fIndVSignalShaping.Response()); break;
      default: pResp = &(fColSignalShaping.Response()); break;
    }

    /* Check the convoluted response, will remove this piece
    std::string iPlaneName[3] = { "u", "v", "y" };
    std::cout << "The convoluted response, before sampling: " << std::endl;
    std::cout << iPlaneName[iplane] << " Field Response: " << std::endl;
    std::cout << "   const int nbin" << iPlaneName[iplane] << "Plane = " << 1200 << ";" << std::endl;
    std::cout << "   const double " << iPlaneName[iplane] << "PlaneResponse[nbin" << iPlaneName[iplane] << "Plane] = {";
    for ( unsigned j = 0; j < 1200; j++ ) {
      if ( j%4 == 0 )  std::cout << std::endl;
      std::cout << "   " << (*pResp)[j] << ", ";
    }
    std::cout << "};" << std::endl;
    */

    for ( int itime = 0; itime < nticks; itime++ ) {
      int low = -1, up = -1;
      for ( int jtime = 0; jtime < nticks; jtime++ ) {
        if ( InputTime[jtime] == SamplingTime[itime] ) {
          SamplingResp[itime] = kDVel * (*pResp)[jtime];
          SamplingCount++;
          break;
        } else if ( InputTime[jtime] > SamplingTime[itime] ) {
          low = jtime - 1;
          up = jtime;
          SamplingResp[itime] = (*pResp)[low] + ( SamplingTime[itime] - InputTime[low] ) * ( (*pResp)[up] - (*pResp)[low] ) / ( InputTime[up] - InputTime[low] );
          SamplingResp[itime] *= kDVel;
          SamplingCount++;
          break;
        } else {
          SamplingResp[itime] = 0.;
        }
      } // for ( int jtime = 0; jtime < nticks; jtime++ )
    } // for ( int itime = 0; itime < nticks; itime++ )

    SamplingResp.resize( SamplingCount, 0.);

    /* Check the convoluted, sampled response, will remove this piece
    std::cout << "The convoluted response, after sampling: " << std::endl;
    std::cout << iPlaneName[iplane] << " Field Response: " << std::endl;
    std::cout << "   const int nbin" << iPlaneName[iplane] << "Plane = " << SamplingResp.size() << ";" << std::endl;
    std::cout << "   const double " << iPlaneName[iplane] << "PlaneResponse[nbin" << iPlaneName[iplane] << "Plane] = {";
    for ( unsigned j = 0; j < SamplingResp.size(); j++ ) {
      if ( j%4 == 0 )  std::cout << std::endl;
      std::cout << "   " << SamplingResp[j] << ", ";
    }
    std::cout << "};" << std::endl;
    */

    switch ( iplane ) {
      case 0: fIndUSignalShaping.AddResponseFunction( SamplingResp, true ); break;
      case 1: fIndVSignalShaping.AddResponseFunction( SamplingResp, true ); break;
      default: fColSignalShaping.AddResponseFunction( SamplingResp, true ); break;
    }


    /* Check the convoluted response, will remove this piece
    const std::vector<double>* pTemp;
    switch ( iplane ) {
      case 0: pTemp = &(fIndUSignalShaping.Response()); break;
      case 1: pTemp = &(fIndVSignalShaping.Response()); break;
      default: pTemp = &(fColSignalShaping.Response()); break;
    }

    std::cout << "The convoluted response, after sampling and pushing back to the SignalShaping object: " << std::endl;
    std::cout << iPlaneName[iplane] << " Field Response: " << std::endl;
    std::cout << "   const int nbin" << iPlaneName[iplane] << "Plane = 1200;" << std::endl;
    std::cout << "   const double " << iPlaneName[iplane] << "PlaneResponse[nbin" << iPlaneName[iplane] << "Plane] = {";
    for ( unsigned j = 0; j < 1200; j++ ) {
      if ( j%4 == 0 )  std::cout << std::endl;
      std::cout << "   " << (*pTemp)[j] << ", ";
    }
    std::cout << "};" << std::endl;
    */

  } // for ( int iplane = 0; iplane < fNPlanes; iplane++ )

  return;
}

namespace util {

  DEFINE_ART_SERVICE(SignalShapingServiceMicroBooNE)

}
