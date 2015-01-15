
////////////////////////////////////////////////////////////////////////
/// \file   SignalShapingServiceMicroBooNE_service.cc
/// \author H. Greenlee
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Utilities/LArFFT.h"
#include "TFile.h"

namespace {

  // loop indices for plane, view, wire, bin
  size_t _pl = 0;
  size_t _vw = 0;
  size_t _wr = 0;
  size_t _bn = 0;
  size_t _ind = 0;
}

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
  // add a comment here
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::LArProperties> larp;
  // Reset initialization flag.

  fInit = false;

  _pl = 0; // just to use it

  fNPlanes = geo->Nplanes();
  fNViews  = pset.get<size_t>("NViews");
  fViewIndex = pset.get<std::vector<size_t> >("ViewIndex");
  for(_vw=0; _vw<fNViews; ++_vw) { fViewMap[_vw] = fViewIndex[_vw]; }
  fViewForNormalization = pset.get<size_t>("ViewForNormalization");
  fNResponses       = pset.get<std::vector<std::vector<size_t> > >("NResponses");
  fNActiveResponses = pset.get<std::vector<std::vector<size_t> > >("NActiveResponses");
  fPrintResponses   = pset.get<bool>("PrintResponses");

  for(size_t ktype=0; ktype<2; ++ktype) {
    // here are some checks
    bool badN = false;

    for(size_t ktype=0; ktype<2; ++ktype) {
      for(_vw=0; _vw<fNViews; ++_vw) {
        if(fNResponses[ktype][_vw] < fNActiveResponses[ktype][_vw])  {
          std::cout
          << "NActiveResponses[" << _vw << "] = " << fNActiveResponses[ktype][_vw] <<
          " > fNResponses[" << _vw << "] = " << fNResponses[ktype][_vw] << std::endl;

          badN = true;
        }
      }
    }
    if(badN) {
      throw art::Exception( art::errors::InvalidNumber )
      << "check NResponses/NActiveRespones" << std::endl;
    }
  }


  // Reset kernels.

  size_t ktype = 0;
  fSignalShapingVec.resize(2);
  std::cout  << fSignalShapingVec.size() << " sets of kernels" << std::endl;
  for(auto& kset : fSignalShapingVec) {
    kset.resize(fNViews);
    _vw = 0;
    for(auto& plane : kset) {
      size_t nWires = fNResponses[ktype][_vw];
      plane.resize(nWires);
      _wr = 0;
      for (auto& ss : plane) {
        ss.Reset();
        _wr++;
      }
      _vw++;
    }
    ktype++;
  }

  fFieldResponseVec.resize(2);
  for(size_t ktype=0; ktype<2; ++ktype) {
    fFieldResponseVec[ktype].resize(fNViews);
    for(_vw=0; _vw<fNViews; ++_vw) {
      size_t nWires = fNResponses[ktype][_vw];
      fFieldResponseVec[ktype][_vw].resize(nWires);
    }
  }

  // Fetch fcl parameters.

  fADCPerPCAtLowestASICGain = pset.get<double>("ADCPerPCAtLowestASICGain");
  fASICGainInMVPerFC = pset.get<double>("ASICGainInMVPerFC");
  fDefaultDriftVelocity = pset.get< DoubleVec >("DefaultDriftVelocity");
  fFieldResponseTOffset = pset.get< std::vector<DoubleVec> >("FieldResponseTOffset");

  for(size_t ktype=0;ktype<2;++ktype) {
    if(fDefaultDriftVelocity.size() != geo->Nplanes() ||
       fFieldResponseTOffset[ktype].size() != geo->Nplanes() )
      throw cet::exception(__FUNCTION__)
      << "\033[93m"
      << "Drift velocity vector and Field response time offset fcl parameter must have length = Nplanes!"
      << "\033[00m" << std::endl;
  }

  fNoiseFactVec =  pset.get<std::vector<DoubleVec> >("NoiseFactVec");
  fCol3DCorrection = pset.get<double>("Col3DCorrection");
  fInd3DCorrection = pset.get<double>("Ind3DCorrection");
  fFieldRespAmpVec = pset.get<DoubleVec>("FieldRespAmpVec");

  fShapeTimeConst = pset.get<std::vector<double> >("ShapeTimeConst");
  fDeconvPol = pset.get<std::vector<int> >("DeconvPol");

  fGetFilterFromHisto= pset.get<bool>("GetFilterFromHisto");

  // Construct parameterized collection filter function.
  if(!fGetFilterFromHisto) {

    fFilterFuncVec.resize(fNViews);
    mf::LogInfo("SignalShapingServiceMicroBooNE") << "Getting Filters from .fcl file" ;

    std::vector<DoubleVec> params = pset.get<std::vector<DoubleVec> >("FilterParamsVec");
    fFilterFuncVec = pset.get<std::vector<std::string> > ("FilterFuncVec");

    fFilterTF1Vec.resize(fNViews);
    for(_vw=0;_vw<fNViews; ++_vw) {
      std::string name = Form("Filter_vw%02i_wr%02i", (int)_vw, (int)_wr);
      fFilterTF1Vec[_vw] = new TF1(name.c_str(), fFilterFuncVec[_vw].c_str() );
      for(_ind=0; _ind<params[_vw].size(); ++_ind) {
        fFilterTF1Vec[_vw]->SetParameter(_ind, params[_vw][_ind]);
      }
    }
  } else {

    std::string histoname = pset.get<std::string>("FilterHistoName");
    mf::LogInfo("SignalShapingServiceMicroBooNE") << " using filter from .root file " ;

    // constructor decides if initialized value is a path or an environment variable
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(pset.get<std::string>("FilterFunctionFname"), fname);

    TFile * in=new TFile(fname.c_str(),"READ");
    for(_vw=0;_vw<fNViews;_vw++){
      std::string name = Form("%s_vw%02i", histoname.c_str(), (int)_vw);
      fFilterHistVec[_vw] = (TH1D *)in->Get(name.c_str());
    }

    in->Close();
    delete in;
  }

  /*
   We allow different drift velocities.
   kDVel is ratio of what was used in LArG4 to field response simulation.
   If drift velocity used for field response is set to <0, then we assume
   the same drift velocity as used in LArG4.
   */
  for(size_t plane = 0; plane < geo->Nplanes(); ++plane) {

    double larg4_velocity = larp->DriftVelocity( larp->Efield(plane), larp->Temperature() );

    if(fDefaultDriftVelocity.at(plane) < 0) fDefaultDriftVelocity.at(plane) = larg4_velocity;

  }

  mf::LogInfo("SignalShapingServiceMicroBooNE") << " using the field response provided from a .root file " ;

  // constructor decides if initialized value is a path or an environment variable
  std::string fileNameBase = pset.get<std::string>("FieldResponseFNameBase");
  std::vector<std::string> version      = pset.get<std::vector<std::string> >("FieldResponseFVersion");
  std::string histNameBase = pset.get<std::string>("FieldResponseHNameBase");
  cet::search_path sp("FW_SEARCH_PATH");

  fFieldResponseHistVec.resize(2);
  for(size_t ktype=0;ktype<2;++ktype) {
    fFieldResponseHistVec[ktype].resize(fNViews);

    _vw = 0;
    std::cout << std::endl;
    for(auto& plane : fFieldResponseHistVec[ktype]) {
      std::string fname0 = Form("%s_vw%02i_%s.root", fileNameBase.c_str(), (int)_vw, version[ktype].c_str());
      std::string fname;
      sp.find_file(fname0, fname);
      //std::cout << "name " << fname0 << std::endl;
      ////std::cout << "path " << fname << std::endl;
      //std::cout << plane.size() << std::endl;
      //std::cout << fNResponses[ktype].size() << std::endl;
      plane.resize(fNResponses[ktype][_vw]);
      std::unique_ptr<TFile> fin(new TFile(fname.c_str(), "READ"));
      _wr = 0;
      std::cout << " fFieldResponseHistVec size, view " << _vw << ": "<< plane.size() << std::endl;
      for(auto& resp : plane) {
        TString histName = Form("%s_vw%02i_wr%02i", histNameBase.c_str(), (int)_vw, (int)_wr);
        //std::cout << histName << std::endl;
        resp = (TH1F*)fin->Get(histName);
        //std::cout << resp->GetNbinsX() << std::endl;

        auto Xaxis = resp->GetXaxis();
        fNFieldBins[ktype] = Xaxis->GetNbins();
        // internal time is in nsec
        fFieldBinWidth[ktype] = resp->GetBinWidth(1)*1000.;
        _wr++;
      }
      fin->Close();
      _vw++;
    }
  }
}

//----------------------------------------------------------------------
// Accessor for single-plane signal shaper.
const util::SignalShaping&
util::SignalShapingServiceMicroBooNE::SignalShaping(unsigned int channel, unsigned int wire, size_t ktype) const
{
  if(!fInit)
    init();

  // Figure out plane type.

  art::ServiceHandle<geo::Geometry> geom;
  //geo::SigType_t sigtype = geom->SignalType(channel);

  auto view = (size_t)geom->View(channel);

  // Return appropriate shaper.

  if(view<fViewIndex[0]||view>=fViewIndex[fNViews]) {
    throw cet::exception("SignalShapingServiceMicroBooNE")<< "can't determine"
    << " View\n";
  }

  // something very wrong here... why doesn't map behave the way it's supposed to???
  //thisView = fViewMap[view];
  return fSignalShapingVec[ktype][view][wire];
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

    SetFilters();

    // Calculate field and electronics response functions.

    std::string kset[2] = { "Convolution ", "Deconvolution "};

    for(size_t ktype=0;ktype<2;++ktype) {

      std::cout << std::endl << kset[ktype] << "functions:" << std::endl;

      // call this first, so that the binning will be known to SetElectResponse
      SetFieldResponse(ktype);
      SetElectResponse(ktype);

      // Configure convolution kernels.

      //Electronic response
      std::cout << "Electonic response " << fElectResponse[ktype].size() << " bins" << std::endl;

      if(fPrintResponses) {
        for(size_t i = 0; i<100; ++i) {
          std::cout << fElectResponse[ktype][i] << " " ;
          if((i+1)%10==0) std::cout << std::endl;
        }
        std::cout << std::endl;
      }

      std::cout << "Input field responses" << std::endl;

      for(_vw=0;_vw<fNViews; ++_vw) {
        for(_wr=0; _wr<fNResponses[ktype][_vw]; ++_wr) {

          if(fPrintResponses) {          std::cout << "Input field response for view " << _vw << " wire " << _wr
            << ", " << (fFieldResponseVec[ktype][_vw][_wr]).size() << " bins" << std::endl;
            for(size_t i = 0; i<(fFieldResponseVec[ktype][_vw][_wr]).size(); ++i) {
              std::cout << fFieldResponseVec[ktype][_vw][_wr][i] << " " ;
              if((i+1)%10==0) std::cout << std::endl;
            }
            std::cout << std::endl;
          }

          (fSignalShapingVec[ktype][_vw][_wr]).AddResponseFunction(fFieldResponseVec[ktype][_vw][_wr]);
          (fSignalShapingVec[ktype][_vw][_wr]).AddResponseFunction(fElectResponse[ktype]);
        }
      }

      // Currently we only have fine binning "fFieldBinWidth"
      // for the field and electronic responses.
      // Now we are sampling the convoluted field-electronic response
      // with the nominal sampling.
      // We may consider to do the same for the filters as well.

      SetResponseSampling(ktype);

            
      // Calculate filter functions.

      // Configure deconvolution kernels.

      for(_vw=0;_vw<fNViews; ++_vw) {
        //std::cout << "filtervec size" << fFilterVec[_vw].size() << std::endl;
        for(_wr=0; _wr<fNResponses[ktype][_vw]; ++_wr) {
          (fSignalShapingVec[ktype][_vw][_wr]).AddFilterFunction(fFilterVec[_vw]);
          (fSignalShapingVec[ktype][_vw][_wr]).SetDeconvKernelPolarity( fDeconvPol.at(_vw));
          (fSignalShapingVec[ktype][_vw][_wr]).CalculateDeconvKernel();
        }
      }
    }
  }
}


//----------------------------------------------------------------------
// Calculate microboone field response.
void util::SignalShapingServiceMicroBooNE::SetFieldResponse(size_t ktype)
{
  // Get services.

  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArProperties> larp;

  //  std::cout << fNResponses.size() << std::endl;
  //  for (auto& ktype : fNResponses) {
  //    std::cout << ktype.size() << std::endl;
  //    for(auto& view : ktype) {
  //      std::cout << view << std::endl;
  //    }
  //  }


  //  //fFieldResponseVec.resize(2);
  //  std::cout << fFieldResponseVec.size() << std::endl;
  //  //fFieldResponseVec[ktype].resize(fNViews);
  //  std::cout << fFieldResponseVec[ktype].size() << std::endl;
  //  _vw = 0;
  //  for(auto& plane : fFieldResponseVec[kt]) {
  //    std::cout << fNResponses[ktype][_wr] << std::endl;
  //    //plane.resize(fNResponses[ktype][_vw]);
  //    std::cout << plane.size() << std::endl;
  //    _vw++;
  //  }



  //double driftvelocity=larp->DriftVelocity()/1000.; // in cm/nsec
  ////////////////////////////////////////////////////

  // Ticks in nanosecond
  // Calculate the normalization of the collection plane
  double integral = 0.;
  integral = fFieldResponseHistVec[ktype][fViewForNormalization][0]->Integral();
  double weight = 1.0/integral;
  //
  //  std::cout << "nResp " << fFieldResponseVec.size() << std::endl;
  //  for(size_t ktype=0;ktype<2;++ktype) {
  //    std::cout << "ktype " << fFieldResponseVec[ktype].size() << std::endl;
  //    for(_vw=0;_vw<fNViews;++_vw) {
  //      std::cout << "view " << (fFieldResponseVec[ktype][_vw]).size() << std::endl;
  //    }
  //  }

  for(_vw=0; _vw<fNViews; ++_vw) {
    for(_wr=0; _wr<fNResponses[ktype][_vw]; ++_wr) {
      size_t nBins = fFieldResponseHistVec[ktype][_vw][_wr]->GetNbinsX();
      (fFieldResponseVec[ktype][_vw][_wr]).resize(nBins);
      for(_bn=1; _bn<=nBins; ++_bn) {
        fFieldResponseVec[ktype][_vw][_wr][_bn-1] =
        fFieldRespAmpVec[_vw]*fFieldResponseHistVec[ktype][_vw][_wr]->GetBinContent(_bn)*weight;
      }
    }
  }

  return;
}


//----------------------------------------------------------------------
// Calculate microboone electronic response.
void util::SignalShapingServiceMicroBooNE::SetElectResponse(size_t ktype)
{
  // Get services.

  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArFFT> fft;

  LOG_DEBUG("SignalShapingMicroBooNE") << "Setting MicroBooNE electronics response function...";

  int nticks = fft->FFTSize();
  std::vector<double> time(nticks,0.);

  fElectResponse.resize(2);
  for(auto& resp : fElectResponse) {
    resp.resize(nticks, 0.);
  }

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

  double max = 0;

  for(int i=0; i<nticks;++i) {
    time[i] = (1.*i) * fFieldBinWidth[ktype]*1.e-3;
  }
  int i = 0;
  for(auto& element :fElectResponse[ktype]) {
    //convert time to microseconds, to match fElectResponse[i] definition
    element =
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

    if(element > max) max = element;
    i++;
  }// end loop over time buckets

  LOG_DEBUG("SignalShapingMicroBooNE") << " Done.";

  // normalize fElectResponse[i], before the convolution
  // Put in overall normalization in a pedantic way:
  // first put in the pulse area per eleectron at the lowest gain setting,
  // then normalize by the actual ASIC gain setting used.
  // This code is executed only during initialization of service,
  // so don't worry about code inefficiencies here.
  double last_integral=0;
  double last_max=0;
  for(auto& element : fElectResponse[ktype]){
    element /= max;
    element *= fADCPerPCAtLowestASICGain * 1.60217657e-7;
    //fElectResponse.at(i) *= fFieldBinWidth / detprop->SamplingRate();
    element *= fASICGainInMVPerFC / 4.7;

    if(element > last_max) last_max = element;
    last_integral += element * fFieldBinWidth[ktype] / detprop->SamplingRate();
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
  size_t nFFT2 = fft->FFTSize() / 2;

  // Calculate collection filter.

  fFilterVec.resize(fNViews);
  for(auto& filter : fFilterVec) {
    filter.resize(nFFT2+1);
  }

  if(!fGetFilterFromHisto)
  {
    _vw = 0;
    for(auto& func : fFilterTF1Vec) {
      func->SetRange(0, double(nFFT2));
      size_t count = 0;
      for(_bn=0; _bn<=nFFT2; ++_bn) {
        //std::cout << "checking TF1 generation " << _bn << " " <<nFFT2 << std::endl;
        double freq = 500.*_bn/(ts*nFFT2);
        double f = func->Eval(freq);
        if(f!=0.0) count++;
        fFilterVec[_vw][_bn] = TComplex(f, 0.);
      }
      //std::cout << count << " non-zero bins out of " << nFFT2 << std::endl;
      _vw++;
    }
  } else{

    _vw = 0;
    for(auto hist : fFilterHistVec) {
      for(_bn=1; _bn<=nFFT2+1; ++_bn) {
        double f = hist->GetBinContent(_bn);
        fFilterVec[_wr][_bn-1] = TComplex(f, 0.);
      }
      _vw++;
    }
  }

}

//----------------------------------------------------------------------
// Sample microboone response (the convoluted field and electronic
// response), will probably add the filter later
void util::SignalShapingServiceMicroBooNE::SetResponseSampling(size_t ktype)
{
  // Get services
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::LArProperties> larp;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArFFT> fft;

  int samplingRate = detprop->SamplingRate();

  /* This could be a warning, but in principle, there's no reason to restrict the binning

   // Operation permitted only if output of rebinning has a larger bin size
   if( fFieldBinWidth > samplingRate )
   throw cet::exception(__FUNCTION__) << "\033[93m"
   << "Invalid operation: cannot rebin to a more finely binned vector!"
   << "\033[00m" << std::endl;

   */

  int nticks = fft->FFTSize();
  std::vector<double> InputTime( nticks, 0. );
  std::vector<double> SamplingTime( nticks, 0. );
  double deltaInputTime = fFieldBinWidth[ktype];
  for ( int itime = 0; itime < nticks; itime++ ) {
    InputTime[itime] = (1.*itime) * deltaInputTime;
    SamplingTime[itime] = (1.*itime) * samplingRate;
    /// VELOCITY-OUT ... comment out kDVel usage here
    //SamplingTime[itime] = (1.*itime) * detprop->SamplingRate() / kDVel;
  }

  // Sampling
  //int fNPlanes = geo->Nplanes();
  //
  std::cout << "Calculating sampled field responses\n";
  for(_vw=0; _vw<fNViews; ++_vw) {
    for(_wr=0; _wr<fNResponses[ktype][_vw]; ++_wr) {
      const std::vector<double>* pResp = &((fSignalShapingVec[ktype][_vw][_wr]).Response());


      std::vector<double> SamplingResp( pResp->size(), 0. );

      /*
       We allow different drift velocities.
       kDVel is ratio of what was used in LArG4 to field response simulation.
       If drift velocity used for field response is set to <0, then we assume
       the same drift velocity as used in LArG4.
       */
      double larg4_velocity = larp->DriftVelocity( larp->Efield(_vw), larp->Temperature() );
      double kDVel = larg4_velocity / fDefaultDriftVelocity.at(_vw);

      // Warning
      if ( kDVel < 1. )
        mf::LogInfo("SignalShapingServiceMicroBooNE") << "The drift velocity "
        << larg4_velocity << "cm/usec is less than the default "
        << fDefaultDriftVelocity.at(_vw) << "cm/usec!"
        << " ... (view=" << _vw << ")"
        << std::endl;

      /*
       Linear interpolation...
       Note from LSR: this is okay, but by restarting jtime at zero for every itime,
       the inner loop gets run way too many times. Since we've already stepped through some
       jtimes, we don't have to do them again.

       Also, just taste, but the "==" case is just a special case of the general
       formula... I prefer to leave it general and simplify the code, assuming that
       it won't be a big effect.
       */

      int SamplingCount = 0;
      //      if(fManualInterpolation) {
      int startJ = 0;
      for ( int itime = 0; itime < nticks; itime++ ) {
        int low = -1, up = -1;
        for ( int jtime = startJ; jtime < nticks; jtime++ ) {
          if ( InputTime[jtime] >= SamplingTime[itime] ) {
            low = jtime - 1;
            up = jtime;
            double interpolationFactor =
            ( (*pResp)[up] - (*pResp)[low] ) /  deltaInputTime;
            SamplingResp[itime] = (*pResp)[low] + ( SamplingTime[itime] - InputTime[low] ) * interpolationFactor;
            /// VELOCITY-OUT ... comment out kDVel usage here
            //SamplingResp[itime] *= kDVel;
            SamplingCount++;
            startJ = jtime;
            break;
          }
        } // for ( int jtime = 0; jtime < nticks; jtime++ )
      } // for ( int itime = 0; itime < nticks; itime++ )
      SamplingResp.resize( SamplingCount, 0.);
      //      } else {
      // the idea here would  be to turn the responses into histograms, and let ROOT do the interpolation.
      // The only advantage is that it's standard code.
      //      }

      if(fPrintResponses) {
        std::cout << "Sampled response (ticks) for view " << _vw << " wire " << _wr << std::endl;
        for(int i = 0; i< 100; ++i) {
          std::cout << SamplingResp[i] << " " ;
          if((i+1)%10==0) std::cout << std::endl;
        }
        std::cout << std::endl;
      }

      (fSignalShapingVec[ktype][_vw][_wr]).AddResponseFunction( SamplingResp, true);
      //std::cout << "Finished with wire " << _wr << ", view " << _vw << std::endl;
      
    }  //  loop over wires
  } // loop over views
  
  std::cout << "Done with field responses" << std::endl;
  return;
}

namespace util {
  
  DEFINE_ART_SERVICE(SignalShapingServiceMicroBooNE)
  
}
