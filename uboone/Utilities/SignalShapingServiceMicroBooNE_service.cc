
////////////////////////////////////////////////////////////////////////
/// \file   SignalShapingServiceMicroBooNE_service.cc
/// \author H. Greenlee
/// Modified by X. Qian 1/6/2015
/// if histogram is used, inialize
/// Response_Offset, Response_Sampling, FieldBins from histogram
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

#include <fstream>

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
: fInit(false), fInitConfigMap(false)
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
  fInitConfigMap = false;

  _pl = 0; // just to use it


  // test
//  fTestParams = pset.get< DoubleVec3 > ("TestParams");
//  // and print it out
//  for(auto& config : fTestParams) {
//    for(auto& plane : config) {
//      for(auto& item : plane) {
//        std::cout << item << " ";
//      }
//      std::cout << " / " ;
//    }
//    std::cout << std::endl;
//  }

  fNConfigs = pset.get<size_t>("NConfigs");
  std::cout << fNConfigs << " TPC ASIC configs are activated" << std::endl;
  //std::cout << "init flag " << fInitConfigMap << std::endl;

  if(fNConfigs>1 && !fInitConfigMap) {

    fConfigMap.clear();
    ifstream configList;

    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file("quietWires.txt", fname);

    configList.open(fname, std::ios::in);
//    if(!configList.isOpen() {
//      std::cout << "file quietWires.txt not found" << std::endl;
//    }

    while(!configList.eof()) {
      size_t item = 10000;
      size_t config;
      configList >> item >> config;
      if (item==10000) break;
      fConfigMap[item] = config;
      //std::cout << item << " " << config << std::endl;
    }
    fInitConfigMap = true;
    std::cout << fConfigMap.size() << " channels read in" << std::endl;

    // now find first and last to speed up search
    if(fConfigMap.size()) {
      fConfigMapFirstChannel = fConfigMap.begin()->first;
      fConfigMapLastChannel  = fConfigMap.rbegin()->first;
    }
    std::cout << "Config map first/last channels: " << fConfigMapFirstChannel << " " << fConfigMapLastChannel << std::endl;
  }

  fASICGainInMVPerFC    = pset.get< DoubleVec2 >("ASICGainInMVPerFC");

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
          /*
           std::cout
           << "NActiveResponses[" << _vw << "] = " << fNActiveResponses[ktype][_vw] <<
           " > fNResponses[" << _vw << "] = " << fNResponses[ktype][_vw] << std::endl;
           */
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

  fSignalShapingVec.resize(fNConfigs);
  for(auto& config : fSignalShapingVec) {
    config.resize(2);
    ktype = 0;
    for(auto& kset : config) {
      kset.resize(fNViews);
      _vw = 0;
      for(auto& plane : kset) {
        size_t nWires = fNResponses[ktype][_vw];
        plane.resize(nWires);
        for (auto& ss : plane) {
          ss.Reset();
        }
        _vw++;
      }
      ktype++;
    }
  }

  fFieldResponseVec.resize(fNConfigs);
  for(auto& config : fFieldResponseVec) {
    config.resize(2);
    int ktype = 0;
    for(auto& kset : config) {
      kset.resize(fNViews);
      _vw = 0;
      for(auto& plane : kset) {
        size_t nWires = fNResponses[ktype][_vw];
        plane.resize(nWires);
      }
      _vw++;
    }
    ktype++;
  }

  // Fetch fcl parameters.
  fDeconNorm = pset.get<double>("DeconNorm");

  fADCPerPCAtLowestASICGain = pset.get<double>("ADCPerPCAtLowestASICGain");
  fDefaultDriftVelocity = pset.get< DoubleVec >("DefaultDriftVelocity");
  fFieldResponseTOffset.resize(2);
  for(size_t ktype=0;ktype<2;++ktype) {
    fFieldResponseTOffset[ktype].resize(fNViews);
  }
  fCalibResponseTOffset = pset.get< DoubleVec >("CalibResponseTOffset");
  std::cout << "CalibResponseTOffsets: ";
  for(auto& x : fCalibResponseTOffset) { std::cout << x << " "; }
  std::cout << std::endl;

  for(size_t ktype=0;ktype<2;++ktype) {
    if(fDefaultDriftVelocity.size() != geo->Nplanes() ||
       fFieldResponseTOffset[ktype].size() != geo->Nplanes() )
      throw cet::exception(__FUNCTION__)
      << "\033[93m"
      << "Drift velocity vector and Field response time offset fcl parameter must have length = Nplanes!"
      << "\033[00m" << std::endl;
  }
  fNoiseFactVec =  pset.get<DoubleVec2>("NoiseFactVec");

  f3DCorrectionVec = pset.get<DoubleVec>("Drift3DCorrVec");

  fFieldRespAmpVec = pset.get<DoubleVec>("FieldRespAmpVec");

  fShapeTimeConst = pset.get<DoubleVec2 >("ShapeTimeConst");
  fDeconvPol = pset.get<std::vector<int> >("DeconvPol");

  fGetFilterFromHisto= pset.get<bool>("GetFilterFromHisto");

  // Construct parameterized collection filter function.
  if(!fGetFilterFromHisto) {

    fFilterFuncVec.resize(fNViews);
    mf::LogInfo("SignalShapingServiceMicroBooNE") << "Getting Filters from .fcl file" ;

    DoubleVec2 params = pset.get< DoubleVec2 >("FilterParamsVec");
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

  // Load 2D filters for induced charge deconvolution (M. Mooney)

  fFilterFuncVecICTime.resize(fNViews);
  fFilterFuncVecICWire.resize(fNViews);
  mf::LogInfo("SignalShapingServiceMicroBooNE") << "Getting 2D Filters from .fcl file" ;

  DoubleVec2 paramsICTime = pset.get< DoubleVec2 >("FilterParamsVecICTime");
  fFilterFuncVecICTime = pset.get<std::vector<std::string> > ("FilterFuncVecICTime");

  fFilterTF1VecICTime.resize(fNViews);
  fFilterICTimeMaxFreq.resize(fNViews);
  fFilterICTimeMaxVal.resize(fNViews);
  for(_vw=0;_vw<fNViews; ++_vw) {
    std::string name = Form("FilterICTime_vw%02i_wr%02i", (int)_vw, (int)_wr);
    fFilterTF1VecICTime[_vw] = new TF1(name.c_str(), fFilterFuncVecICTime[_vw].c_str() );
    for(_ind=0; _ind<paramsICTime[_vw].size(); ++_ind) {
      fFilterTF1VecICTime[_vw]->SetParameter(_ind, paramsICTime[_vw][_ind]);
      fFilterICTimeMaxFreq[_vw] = fFilterTF1VecICTime[_vw]->GetMaximumX();
      fFilterICTimeMaxVal[_vw] = fFilterTF1VecICTime[_vw]->GetMaximum();
    }
  }

  DoubleVec2 paramsICWire = pset.get< DoubleVec2 >("FilterParamsVecICWire");
  fFilterFuncVecICWire = pset.get<std::vector<std::string> > ("FilterFuncVecICWire");

  fFilterTF1VecICWire.resize(fNViews);
  fFilterICWireMaxFreq.resize(fNViews);
  fFilterICWireMaxVal.resize(fNViews);
  for(_vw=0;_vw<fNViews; ++_vw) {
    std::string name = Form("FilterICWire_vw%02i_wr%02i", (int)_vw, (int)_wr);
    fFilterTF1VecICWire[_vw] = new TF1(name.c_str(), fFilterFuncVecICWire[_vw].c_str() );
    for(_ind=0; _ind<paramsICWire[_vw].size(); ++_ind) {
      fFilterTF1VecICWire[_vw]->SetParameter(_ind, paramsICWire[_vw][_ind]);
      fFilterICWireMaxFreq[_vw] = fFilterTF1VecICWire[_vw]->GetMaximumX();
      fFilterICWireMaxVal[_vw] = fFilterTF1VecICWire[_vw]->GetMaximum();
    }
  }

  fFilterScaleVecICTime = pset.get<DoubleVec>("FilterScaleVecICTime");
  fFilterScaleVecICWire = pset.get<DoubleVec>("FilterScaleVecICWire");
  fFilterNormVecIC = pset.get<DoubleVec>("FilterNormVecIC");

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

  //Adding calibrated field response at 70kV
  fUseCalibratedResponses = pset.get<bool>("UseCalibratedResponses");

  mf::LogInfo("SignalShapingServiceMicroBooNE") << " using the field response provided from a .root file " ;

  // constructor decides if initialized value is a path or an environment variable
  std::string fileNameBase = pset.get<std::string>("FieldResponseFNameBase");
  std::vector<std::string> version      = pset.get<std::vector<std::string> >("FieldResponseFVersion");
  fDefaultEField                 = pset.get<double>("DefaultEField");
  fDefaultTemperature            = pset.get<double>("DefaultTemperature");

  fTimeScaleParams               = pset.get<DoubleVec>("TimeScaleParams");

  std::string histNameBase = pset.get<std::string>("FieldResponseHNameBase");
  cet::search_path sp("FW_SEARCH_PATH");

  DoubleVec tOffset(fNViews, 0.0);
  //tOffset2;

  fFieldResponseHistVec.resize(fNConfigs);
  for(size_t config=0;config<fNConfigs; ++config) {
    fFieldResponseHistVec[config].resize(2);
    for(size_t ktype=0;ktype<2;++ktype) {
      fFieldResponseHistVec[config][ktype].resize(fNViews);
      
      // calculate the time scale factor for this event
      if(config==0 && ktype==0 && !fUseCalibratedResponses) SetTimeScaleFactor();

      _vw = 0;
      //std::cout << "config " << config << " in " << fNConfigs << std::endl;
      for(auto& plane : fFieldResponseHistVec[config][ktype]) {
        std::string fname0 = Form("%s_vw%02i_%s.root", fileNameBase.c_str(), (int)_vw, version[ktype].c_str());
        std::string fname;
        sp.find_file(fname0, fname);
        plane.resize(fNResponses[ktype][_vw]);
        std::unique_ptr<TFile> fin(new TFile(fname.c_str(), "READ"));
        _wr = 0;
        // load up the response functions
        for(auto& resp : plane) {
          TString histName = Form("%s_vw%02i_wr%02i", histNameBase.c_str(), (int)_vw, (int)_wr);
          resp = (TH1F*)fin->Get(histName);

          auto Xaxis = resp->GetXaxis();
          fNFieldBins[ktype] = Xaxis->GetNbins();

          // internal time is in nsec
          fFieldBinWidth[ktype] = resp->GetBinWidth(1)*1000.;
          // get the offsets for each plane... use wire 0 and either peak or zero-crossing

          SetFieldResponseTOffsets(resp, ktype);

          _wr++;
        }
        fin->Close();
        _vw++;
      }
    }
  }
}

void util::SignalShapingServiceMicroBooNE::SetTimeScaleFactor()
{
  // get the scale factor between the bulk drift velocity used to generate the field response
  //   and that used for this simulation.

  art::ServiceHandle<util::LArProperties> larp;

  double defaultVelocity = larp->DriftVelocity(fDefaultEField, fDefaultTemperature);
  double thisVelocity    = larp->DriftVelocity( larp->Efield(0), larp->Temperature() );
  double vRatio = defaultVelocity/thisVelocity;
  double vDiff = vRatio -1.0;

  fTimeScaleFactor = 0.0;
  double term = 1.0;

  // the time scale params are from a fit to Garfield simulations at different E Fields
  for(size_t i = 0;i<fTimeScaleParams.size(); ++i) {
    fTimeScaleFactor += fTimeScaleParams[i]*term;
    term *= vDiff;
  }

  std::cout << "Current E field = " << larp->Efield(0) << " KV/cm, Ratio of drift velocities = " << vRatio << ", timeScaleFactor = " << fTimeScaleFactor << std::endl;

}

//-----------------------------
// Extract the time offsets from the field-response histograms
void util::SignalShapingServiceMicroBooNE::SetFieldResponseTOffsets(const TH1F* resp, const size_t ktype)
{
  double tOffset = 0.0;

  // this is for the calibrated response, for now we disable this so that the timing
  //   for the standard and BNL responses is the same.

  /*
  if(fUseCalibratedResponses) {
    double calibResponseTOffset = fCalibResponseTOffset[_vw];
    fFieldResponseTOffset[ktype].at(_vw) = (-tOffset + calibResponseTOffset)*1000.0;
    return;
  }
  */ 

  // this is for the standard response

  if(_wr==0 && _vw==fViewForNormalization) {
    // for the collection plane, find the peak
    int binMax = resp->GetMaximumBin();
    tOffset = (resp->GetXaxis()->GetBinCenter(binMax) - resp->GetXaxis()->GetBinCenter(1));
    // for later, to be a bit cleverer, take weighted average of 3 bins
    //          for(int bin=binMax-1; bin<=binMax+1; ++ bin) {
    //            content = resp->GetBinContent(bin);
    //            binVal = resp->GetXaxis()->GetBinCenter(bin);
    //            numer += content*binVal;
    //            denom += content;
    //          }
    //          tOffset[_vw] = numer/denom*delta - resp->GetXaxis()->GetBinCenter(1);
  } else {
    // for the other planes, find the zero-crossing
    // lets find the minimum, and work backwards!

    int binMin = resp->GetMinimumBin();
    for(int bin=binMin;bin>0; --bin) {
      double content = resp->GetBinContent(bin);
      bool found = false;
      if(content>0) {
        double binVal = resp->GetXaxis()->GetBinCenter(bin);
        tOffset = binVal - resp->GetXaxis()->GetBinCenter(1);
        found = true;
        // for later
        //            } else if (content>0) {
        //              // If it's already gone through zero, split the difference
        //              std::cout << resp->GetBinContent(bin) << " " << resp->GetBinContent(bin+1) << " " << bin << std::endl;
        //              binVal = resp->GetXaxis()->GetBinCenter(bin);
        //              numer = resp->GetBinContent(bin)*binVal - resp->GetBinContent(bin+1)*(binVal+delta);
        //              denom = resp->GetBinContent(bin) - resp->GetBinContent(bin+1);
        //              found = true;
      }
      if(found) {
        //              std::cout << numer << " " << denom << " " << delta << std::endl;
        //              tOffset[_vw] = numer/denom*delta - resp->GetXaxis()->GetBinCenter(1);
        break;
      }
    }
  }

  //std::cout << "view " << _vw << ", wire " << _wr << ", toffset " << tOffset << std::endl;
  tOffset *= f3DCorrectionVec[_vw]*fTimeScaleFactor;
  fFieldResponseTOffset[ktype].at(_vw) = (-tOffset+ fCalibResponseTOffset[_vw])*1000.;

}

//----------------------------------------------------------------------
// Accessor for single-plane signal shaper.
const util::SignalShaping&
util::SignalShapingServiceMicroBooNE::SignalShaping(size_t channel, size_t wire, size_t ktype) const
{
  if(!fInit)
    init();

  // Figure out plane type.

  art::ServiceHandle<geo::Geometry> geom;
  //geo::SigType_t sigtype = geom->SignalType(channel);

  auto view = (size_t)geom->View(channel);

  // Return appropriate shaper.

  if(view<fViewIndex[0]||view>fViewIndex[fNViews-1]) {
    throw cet::exception("SignalShapingServiceMicroBooNE")<< "can't determine"
    << " View\n";
  }

  size_t config = GetConfig(channel);

  // something very wrong here... why doesn't map behave the way it's supposed to???
  //thisView = fViewMap[view];

  //if(ktype==0&&view==0) std::cout << "SS Params " << channel << " " << config << " " << std::endl;

  return fSignalShapingVec[config][ktype][view][wire];
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



    // re-initialize the FFT service for the request size
    art::ServiceHandle<util::LArFFT> fFFT;
    std::string options = fFFT->FFTOptions();
    int fitbins = fFFT->FFTFitBins();
    int fftsize = fFFT->FFTSize();


    // Calculate field and electronics response functions.


    std::string kset[2] = { "Convolution ", "Deconvolution "};


    for(size_t config=0;config<fNConfigs;++config) {
      for(size_t ktype=0;ktype<2;++ktype) {
        if (fNFieldBins[ktype]*4>fftsize)
          fFFT->ReinitializeFFT( fNFieldBins[ktype]*4, options, fitbins);
        //std::cout << std::endl << kset[ktype] << "functions:" << std::endl;

        // call this first, so that the binning will be known to SetElectResponse
        SetFieldResponse(ktype);


        //std::cout << "Input field responses" << std::endl;

        for(_vw=0;_vw<fNViews; ++_vw) {
          SetElectResponse(ktype,fShapeTimeConst[config].at(_vw),fASICGainInMVPerFC[config].at(_vw));

          //Electronic response
          //std::cout << "Electonic response " << fElectResponse[ktype].size() << " bins" << std::endl;

          if(fPrintResponses) {
            for(size_t i = 0; i<100; ++i) {
              std::cout << "Electronic Response " << fElectResponse[ktype][i] << " " ;
              if((i+1)%10==0) std::cout << std::endl;
            }
            std::cout << std::endl;
          }

          for(_wr=0; _wr<fNResponses[ktype][_vw]; ++_wr) {

            if(fPrintResponses) {          std::cout << "Input field response for view " << _vw << " wire " << _wr
              << ", " << (fFieldResponseVec[config][ktype][_vw][_wr]).size() << " bins" << std::endl;
              for(size_t i = 0; i<(fFieldResponseVec[config][ktype][_vw][_wr]).size(); ++i) {
                std::cout << fFieldResponseVec[config][ktype][_vw][_wr][i] << " " ;
                if((i+1)%10==0) std::cout << std::endl;
              }
              std::cout << std::endl;
            }

            (fSignalShapingVec[config][ktype][_vw][_wr]).AddResponseFunction(fFieldResponseVec[config][ktype][_vw][_wr]);
            (fSignalShapingVec[config][ktype][_vw][_wr]).AddResponseFunction(fElectResponse[ktype]);
            (fSignalShapingVec[config][ktype][_vw][_wr]).save_response();
            (fSignalShapingVec[config][ktype][_vw][_wr]).set_normflag(false);
          }
        }
        SetResponseSampling(ktype, config);

        // Currently we only have fine binning "fFieldBinWidth"
        // for the field and electronic responses.
        // Now we are sampling the convoluted field-electronic response
        // with the nominal sampling.
        // We may consider to do the same for the filters as well.
        if (fftsize!=fFFT->FFTSize()){
          std::string options = fFFT->FFTOptions();
          int fitbins = fFFT->FFTFitBins();
          fFFT->ReinitializeFFT( fftsize, options, fitbins);
        }


        // Calculate filter functions.
        if(config==0) SetFilters();

        // Configure deconvolution kernels.

        for(_vw=0;_vw<fNViews; ++_vw) {
          // std::cout << "filtervec size" << fFilterVec[_vw].size() << std::endl;
          for(_wr=0; _wr<fNResponses[ktype][_vw]; ++_wr) {
            (fSignalShapingVec[config][ktype][_vw][_wr]).AddFilterFunction(fFilterVec[_vw]);
            (fSignalShapingVec[config][ktype][_vw][_wr]).SetDeconvKernelPolarity( fDeconvPol.at(_vw));
            (fSignalShapingVec[config][ktype][_vw][_wr]).CalculateDeconvKernel();
          }
        }
      }
    }
  }
}

void util::SignalShapingServiceMicroBooNE::SetDecon(size_t fftsize, size_t channel)
{
  art::ServiceHandle<geo::Geometry> geo;

  //std::cout << "enter SetDecon, init flag "  << fInit <<  " fftsize " << fftsize << " channel " << channel << std::endl;

  init();

  
  art::ServiceHandle<util::LArFFT> fFFT;
  
  // streamline this method:
  // if the deconvolution kernel is already appropriate for the datasize (aka fftsize) do nothing
  // otherwise, set it to the appropriate size
  // do this test for *every* ss
  // But it will in general only happen once per run!
  
  bool setDecon = false;
  
  size_t FFTSize = fFFT->FFTSize();
  if (fftsize>FFTSize||fftsize<=FFTSize/2){
    std::string options = fFFT->FFTOptions();
    int fitbins = fFFT->FFTFitBins();
    fFFT->ReinitializeFFT( fftsize, options, fitbins);
    setDecon = true;
  }
  
  if(!setDecon) return;
  
  size_t ktype = 1;
  //std::cout << "nconfigs/nviews " << fNConfigs << " " << fNViews << std::endl;
  
  for(size_t config=0; config<fNConfigs; ++config) {

    for (size_t view=0;view<fNViews; ++view) {
      
      //size_t config = GetConfig(channel);
      //geo::View_t view = geo->View(channel);
      //size_t ktype = 1;
      
      //std::cout << "view/_vw " << view << " " << _vw << std::endl;
      
      
      for(_wr=0; _wr<fNResponses[ktype][view]; ++_wr) {
        (fSignalShapingVec[config][ktype][view][_wr]).Reset();
      }
    }
    //for(size_t view=0; view<fNViews; ++view) {
      
      //std::cout << "about to call SetResponseSampling" << std::endl;
      int mode = 0;
      SetResponseSampling(ktype, config, mode, channel);
    //}

      //std::cout << "Xin2 " << std::endl;
      // Calculate filter functions.
      if(config==0) {
        //std::cout << "set the filters" << std::endl;
        SetFilters();
      }
      // Configure deconvolution kernels.
      //std::cout << "Xin3 " << std::endl;
      //std::cout << "FInish the SS" << std::endl;
     
    for(size_t view=0; view < fNViews; ++view) {
      for(_wr=0; _wr<fNResponses[ktype][view]; ++_wr) {
        //std::cout << "this wire " << _wr << std::endl;
        (fSignalShapingVec[config][ktype][view][_wr]).AddFilterFunction(fFilterVec[view]);
        (fSignalShapingVec[config][ktype][view][_wr]).SetDeconvKernelPolarity( fDeconvPol.at(view));
        (fSignalShapingVec[config][ktype][view][_wr]).CalculateDeconvKernel();
        //std::cout << "Xin4 " << std::endl;
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

  //  }


  ////////////////////////////////////////////////////

  // Ticks in nanosecond
  // Calculate the normalization of the collection plane
  double integral;
  double weight;

  for(size_t config=0; config<fNConfigs; ++config) {

    integral = fFieldResponseHistVec[config][ktype][fViewForNormalization][0]->Integral();
    weight = 1./integral;

    for(_vw=0; _vw<fNViews; ++_vw) {
      for(_wr=0; _wr<fNResponses[ktype][_vw]; ++_wr) {
        size_t nBins = fFieldResponseHistVec[config][ktype][_vw][_wr]->GetNbinsX();
        (fFieldResponseVec[config][ktype][_vw][_wr]).resize(nBins);
        for(_bn=1; _bn<=nBins; ++_bn) {
          fFieldResponseVec[config][ktype][_vw][_wr][_bn-1] =
            fFieldResponseHistVec[config][ktype][_vw][_wr]->GetBinContent(_bn);
          fFieldResponseVec[config][ktype][_vw][_wr][_bn-1]
            *= fFieldRespAmpVec[_vw]*weight;
        }
      }
    }
  }


  return;
}


//----------------------------------------------------------------------
// Calculate microboone field response.
void util::SignalShapingServiceMicroBooNE::SetElectResponse(size_t ktype,double shapingtime, double gain)
{
  // Get services.

  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArFFT> fft;

  LOG_DEBUG("SignalShapingMicroBooNE") << "Setting MicroBooNE electronics response function...";

  size_t nticks = fft->FFTSize();
  DoubleVec time(nticks,0.);

  fElectResponse.resize(2);
  for(auto& resp : fElectResponse) {
    resp.resize(nticks, 0.);
  }

  //Gain and shaping time variables from fcl file:
  double Ao = 1.0;//Gain
  double To = shapingtime;  //peaking time

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


  // For the cold electronics,  the gain (i.e. 4.7 mV/fC) represents the peak 
  // height. The shaping time will not affect the peak height, but make the 
  // peak broader

  double max = 0;

  for(size_t i=0; i<nticks;++i) {
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

  //Normalization are the following
  // Peak is firstly normalized to 1
  // thus we expect peak to be 1 * 9390 (fADCPerPCtAtLowestAsicGain) * 1.602e-7 * (1 fC) = 9.39 ADC
  // At 4.7 mV/fC, the ADC value should be 4.7 (mV/fC) * 2 (ADC/mV) ~ 9.4 ADC/fC
  // so the normalization are consistent

  

  for(auto& element : fElectResponse[ktype]){
    element /= max;
    element *= fADCPerPCAtLowestASICGain * 1.60217657e-7;
    element *= gain / 4.7;

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
void util::SignalShapingServiceMicroBooNE::SetResponseSampling(size_t ktype, size_t config, int mode, size_t channel)
{
  // Get services
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::LArProperties> larp;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArFFT> fft;

  /* This could be a warning, but in principle, there's no reason to restrict the binning

   // Operation permitted only if output of rebinning has a larger bin size
   if( fFieldBinWidth > samplingRate )
   throw cet::exception(__FUNCTION__) << "\033[93m"
   << "Invalid operation: cannot rebin to a more finely binned vector!"
   << "\033[00m" << std::endl;

   */
  //std::cout << "entering SetResponseSampling, ktype/config/mode/channel " << ktype << " " << config << " " << mode << " " << channel << std::endl;

  size_t view0, view1;
  if(mode==0) {
    view0 = 0;
    view1 = fNViews;
  } else {
    geo::View_t view = geo->View(channel);
    view0 = view;
    view1 = std::min(fNViews,(size_t)view+1);
  }

  //std::cout << "view0/1 " << view0 << " " << view1 << std::endl;

  size_t nticks = fft->FFTSize();
  DoubleVec SamplingTime( nticks, 0. );
  double deltaInputTime = fFieldBinWidth[ktype];
  //std::cout << "nticks = " << nticks << std::endl;
  for ( size_t itime = 0; itime < nticks; itime++ ) {
    SamplingTime[itime] = (1.*itime) * detprop->SamplingRate();
  }

  // Sampling


    for(size_t view=view0; view<view1; ++view) {
      //std::cout << "sampling view " << view  << " ktype/config/channel " << ktype << " " << config << " " << channel << std::endl;

      double timeFactor = fTimeScaleFactor*f3DCorrectionVec[_wr];
      // time factor is already included in the calibrated response
      if(fUseCalibratedResponses) timeFactor = 1.0;
      double timeFactorInv = 1./timeFactor;
      for(_wr=0; _wr<fNResponses[ktype][view]; ++_wr) {
        const DoubleVec* pResp = &((fSignalShapingVec[config][ktype][view][_wr]).Response_save());

        size_t nticks_input = pResp->size();
        DoubleVec InputTime(nticks_input, 0. );
        for (size_t itime = 0; itime < nticks_input; itime++ ) {
          InputTime[itime] = (1.*itime) * deltaInputTime*timeFactor;
        }
        //std::cout << "Input time vector done" << std::endl;

        DoubleVec SamplingResp(nticks, 0. );

        size_t SamplingCount = 0;

        size_t startJ = 1;
        SamplingResp[0] = (*pResp)[0];
        for ( size_t itime = 1; itime < nticks; itime++ ) {
          size_t low, high;
          for ( size_t jtime = startJ; jtime < nticks_input; jtime++ ) {
            if ( InputTime[jtime] >= SamplingTime[itime] ) {
              low  = jtime - 1;
              high = jtime;
              //            if(jtime<2&&itime<2) std::cout << itime << " " << jtime << " " << low << " " << up << std::endl;
              double interpolationFactor = ((*pResp)[high]-(*pResp)[low])/deltaInputTime;
              SamplingResp[itime] = ((*pResp)[low] + ( SamplingTime[itime] - InputTime[low] ) * interpolationFactor);
              // note: timeFactor = timeFactorInv =  1.0 for calibrated responses
              SamplingResp[itime] *= timeFactorInv;
              SamplingCount++;
              startJ = jtime;
              break;
            }
          } // for (  jtime = 0; jtime < nticks; jtime++ )
        } // for (  itime = 0; itime < nticks; itime++ )
        //std::cout << "SamplingResponse done " << std::endl;

        if(fPrintResponses) {
          size_t printCount = 0;
          int inc = 1;
          //std::cout << "Sampled response (ticks) for view " << view << " wire " << _wr << " nticks " << nticks << std::endl;
          for(size_t i = 0; i<nticks; i+=inc) {
            std::cout << SamplingResp[i] << " " ;
            if((printCount+1)%10==0) std::cout << std::endl;
            printCount++;
            if (printCount>=100) {inc = 100;}
          }
        }

        (fSignalShapingVec[config][ktype][view][_wr]).AddResponseFunction( SamplingResp, true);
        //std::cout << "Finished with wire " << _wr << ", view " << _vw << std::endl;
        
      }  //  loop over wires
    } // loop over views

  //std::cout << "Done with field responses" << std::endl;
  return;
}


//-----Give Gain Settings to SimWire-----//jyoti
double util::SignalShapingServiceMicroBooNE::GetASICGain(unsigned int const channel) const
{
  art::ServiceHandle<geo::Geometry> geom;
  geo::View_t view = geom->View(channel);
  double gain = 0;
  size_t config = GetConfig(channel);
  switch(view){
    case geo::kU:
      gain = fASICGainInMVPerFC[config].at(0);
      break;
    case geo::kV:
      gain = fASICGainInMVPerFC[config].at(1);
      break;
    case geo::kZ:
      gain = fASICGainInMVPerFC[config].at(2);
      break;
    default:
      throw cet::exception(__FUNCTION__) << "Invalid geo::View_t ... " << view << std::endl;
  }
  return gain;
}


//-----Give Shaping time to SimWire-----//jyoti
double util::SignalShapingServiceMicroBooNE::GetShapingTime(unsigned int const channel) const
{
  art::ServiceHandle<geo::Geometry> geom;
  geo::View_t view = geom->View(channel);

  size_t config = GetConfig(channel);

  double shaping_time(0);
  switch(view){
    case geo::kU:
      shaping_time = fShapeTimeConst[config].at(0);
      break;
    case geo::kV:
      shaping_time = fShapeTimeConst[config].at(1);
      break;
    case geo::kZ:
      shaping_time = fShapeTimeConst[config].at(2);
      break;
    default:
      throw cet::exception(__FUNCTION__) << "Invalid geo::View_t ... " << view << std::endl;
  }
  return shaping_time;
}

double util::SignalShapingServiceMicroBooNE::GetRawNoise(unsigned int const channel) const
{
  unsigned int plane;
  art::ServiceHandle<geo::Geometry> geom;
  geo::View_t view = geom->View(channel);

  size_t config = GetConfig(channel);
  switch(view){
    case geo::kU:
      plane = 0;
      break;
    case geo::kV:
      plane = 1;
      break;
    case geo::kZ:
      plane = 2;
      break;
    default:
      throw cet::exception(__FUNCTION__) << "Invalid geo::View_t ... " << view << std::endl;
  }

  double shapingtime = fShapeTimeConst[config].at(plane);
  double gain = fASICGainInMVPerFC[config].at(plane);
  int temp;
  if (shapingtime == 0.5){
    temp = 0;
  }else if (shapingtime == 1.0){
    temp = 1;
  }else if (shapingtime == 2.0){
    temp = 2;
  }else{
    temp = 3;
  }
  double rawNoise;

  auto tempNoise = fNoiseFactVec.at(plane);
  rawNoise = tempNoise.at(temp);

  rawNoise *= gain/4.7;
  return rawNoise;
}

double util::SignalShapingServiceMicroBooNE::GetDeconNoise(unsigned int const channel) const
{
  unsigned int plane;
  art::ServiceHandle<geo::Geometry> geom;
  geo::View_t view = geom->View(channel);

  size_t config = GetConfig(channel);
  switch(view){
    case geo::kU:
      plane = 0;
      break;
    case geo::kV:
      plane = 1;
      break;
    case geo::kZ:
      plane = 2;
      break;
    default:
      throw cet::exception(__FUNCTION__) << "Invalid geo::View_t ... " << view << std::endl;
  }

  double shapingtime = fShapeTimeConst[config].at(plane);
  int temp;
  if (shapingtime == 0.5){
    temp = 0;
  }else if (shapingtime == 1.0){
    temp = 1;
  }else if (shapingtime == 2.0){
    temp = 2;
  }else{
    temp = 3;
  }
  auto tempNoise = fNoiseFactVec.at(plane);
  double deconNoise = tempNoise.at(temp);

  deconNoise = deconNoise /4096.*2000./4.7 *6.241*1000/fDeconNorm;
  return deconNoise;
}

int util::SignalShapingServiceMicroBooNE::FieldResponseTOffset(unsigned int const channel, size_t ktype=0) const
{
  art::ServiceHandle<geo::Geometry> geom;
  geo::View_t view = geom->View(channel);

  double time_offset = 0;
  switch(view){
    case geo::kU:
      time_offset = fFieldResponseTOffset[ktype].at(0);
      break;
    case geo::kV:
      time_offset = fFieldResponseTOffset[ktype].at(1);
      break;
    case geo::kZ: 
      time_offset = fFieldResponseTOffset[ktype].at(2); 
      break;
    default:
      throw cet::exception(__FUNCTION__) << "Invalid geo::View_t ... " << view << std::endl;
  }
  auto tpc_clock = art::ServiceHandle<util::TimeService>()->TPCClock();
  return tpc_clock.Ticks(time_offset/1.e3);
}

//----------------------------------------------------------------------
// Get convolution kernel from SignalShaping service for use in CalWire's
// DeconvoluteInducedCharge() - added by M. Mooney
const std::vector<TComplex>& util::SignalShapingServiceMicroBooNE::GetConvKernel(unsigned int channel, unsigned int wire) const
{
  if(!fInit)
    init();

  art::ServiceHandle<geo::Geometry> geom;
  //geo::SigType_t sigtype = geom->SignalType(channel);

  auto view = (size_t)geom->View(channel);

  size_t config = GetConfig(channel);
  // Return appropriate shaper.

  if(view<fViewIndex[0]||view>fViewIndex[fNViews-1]) {
    throw cet::exception("SignalShapingServiceMicroBooNE")<< "can't determine"
    << " View\n";
  }

  return fSignalShapingVec[config][0][view][wire].ConvKernel();
}

//----------------------------------------------------------------------
// Evaluate 2D filter used in induced charge deconvolution (M. Mooney)

double util::SignalShapingServiceMicroBooNE::Get2DFilterVal(size_t planeNum, size_t freqDimension, double binFrac) const
{
  art::ServiceHandle<util::DetectorProperties> detprop;
  double ts = detprop->SamplingRate();

  double freq;
  double filtVal;
  double val;
  if(freqDimension == 1) {
    freq = fFilterScaleVecICTime.at(planeNum)*(500.0/ts)*2.0*(0.5-fabs(binFrac-0.5));
    filtVal = fFilterTF1VecICTime.at(planeNum)->Eval(freq);

    if(freq < fFilterICTimeMaxFreq.at(planeNum))
      val = fFilterICTimeMaxVal.at(planeNum);
    else
      val = filtVal;

    return val;
  }
  else if(freqDimension == 2) {
    freq = fFilterScaleVecICWire.at(planeNum)*(0.5-fabs(binFrac-0.5));
    filtVal = fFilterTF1VecICWire.at(planeNum)->Eval(freq);

    if(freq < fFilterICWireMaxFreq.at(planeNum))
      val = fFilterICWireMaxVal.at(planeNum);
    else
      val = filtVal;

    return val;
  }
  else {
    return 0.0;
  }
}

//----------------------------------------------------------------------
// Return the correct configuration for this channel
size_t util::SignalShapingServiceMicroBooNE::GetConfig(size_t channel) const
{
  if(fNConfigs<=1 || fConfigMap.size()==0) return 0;
  if(channel>fConfigMapLastChannel || channel<fConfigMapFirstChannel) return 0;

  // for a test with the special sim event, set n to 1000

  if(fConfigMap.find(channel)==fConfigMap.end()) return 0;

  return fConfigMap[channel];
}

//----------------------------------------------------------------------
// Return 2D filter normalization, used in induced charge deconvolution (M. Mooney)
double util::SignalShapingServiceMicroBooNE::Get2DFilterNorm(size_t planeNum) const
{
  return fFilterNormVecIC.at(planeNum);
}

namespace util {
  
  DEFINE_ART_SERVICE(SignalShapingServiceMicroBooNE)
  
}

