////////////////////////////////////////////////////////////////////////
// $Id: MicrobooneResp.cxx,v 1.0 2010/09/15  bpage Exp $
//
// MicrobooneResp class
//
// pagebri3@msu.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef MICROBOONERESP_H
#define MICROBOONERESP_H

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "Utilities/LArFFT.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"

// ROOT includes
#include <TFile.h>
#include <TH2D.h>
#include <TF1.h>
#include <TComplex.h>

namespace uboone{

  class MicrobooneResp : public art::EDProducer { 
  
  public:
 
    explicit MicrobooneResp(fhicl::ParameterSet const& pset);
  
    virtual ~MicrobooneResp();
  
    void reconfigure(const fhicl::ParameterSet &p);

    void produce(art::Event & evt);
 
    void beginJob();
    
  private:
    
    double fInd3DCorrection;
    double fCol3DCorrection;
    double fIndSimScale;
    double fColSimScale;
    double fAsymmetry;
    std::vector<double> fEParams;
    std::vector<double> fNParams;
    std::vector<double> fRIFParams;
    std::vector<double> fRCFParams;
    std::vector<double> fSIFParams;
    std::vector<double> fSCFParams;
    std::string fElectFunc;
    std::string fNoiseFunc;
    std::string fCDriftShape;
    std::string fIDriftShape;
    std::string fRIFFunc;
    std::string fRCFFunc;
    std::string fSIFFunc;
    std::string fSCFFunc;

  protected:
  }; //class MicrobooneResp
}

namespace uboone{

  MicrobooneResp::MicrobooneResp(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }
  
  MicrobooneResp::~MicrobooneResp()
  {
  }
 
  void MicrobooneResp::reconfigure(const fhicl::ParameterSet &p)
  {
    fInd3DCorrection = p.get< double >               ("Ind3DCorrection");
    fCol3DCorrection = p.get< double >               ("Col3DCorrection");
    fAsymmetry       = p.get< double >               ("Asymmetry");
    fIndSimScale     = p.get< double >               ("IndSimScale");
    fColSimScale     = p.get< double >               ("ColSimScale");
    fEParams         = p.get< std::vector< double > >("ElectParams");
    fNParams         = p.get< std::vector< double > >("NoiseParams");
    fRIFParams       = p.get< std::vector< double > >("RawIndFilterParams");
    fRCFParams       = p.get< std::vector< double > >("RawColFilterParams");
    fSIFParams       = p.get< std::vector< double > >("SimIndFilterParams");
    fSCFParams       = p.get< std::vector< double > >("SimColFilterParams");
    fElectFunc       = p.get< std::string >          ("ElectFunction");
    fNoiseFunc       = p.get< std::string >          ("NoiseFunction");
    fCDriftShape     = p.get< std::string >          ("ColDriftShape");
    fIDriftShape     = p.get< std::string >          ("IndDriftShape");
    fRIFFunc         = p.get< std::string >          ("RawIndFilter");
    fRCFFunc         = p.get< std::string >          ("RawColFilter");
    fSIFFunc         = p.get< std::string >          ("SimIndFilter");
    fSCFFunc         = p.get< std::string >          ("SimColFilter");
  }

  void MicrobooneResp::produce(art::Event & evt)
  {
  }

  void MicrobooneResp::beginJob()
  {
    art::ServiceHandle<util::LArFFT> FFT;
    art::ServiceHandle<geo::Geometry> geo;
    int signalSize = FFT->FFTSize();
    int nPlanes = geo->Nplanes();
    int nChannels = geo->Nchannels();
    int wires =0;
    int iRange =0;
    int nRanges = 0;
    std::vector<double> ranges;
    ranges.push_back(0);
    for(int p = 0; p < nPlanes; p++) {
      nRanges++;
      wires+=geo->Plane(p).Nwires();
      std::cout << wires << std::endl;
      if(p == nPlanes-2) iRange=wires;
      ranges.push_back(wires);
    }    
    std::vector<double> diffsig(signalSize);
    std::vector<double> ramp(signalSize);
    TComplex kernBin;
    int size = signalSize/2;
    int bin=0;
    std::vector<TComplex> freqSig(size+1);
    std::vector<double> bipolar(signalSize);
    TF1 rampf("rampf",fCDriftShape.c_str(),0,signalSize);
    rampf.SetParameters(8.0,16.0,.75);
    TF1 bipolf("bipolf",fIDriftShape.c_str(),0,signalSize);
    bipolf.SetParameters(.00843,.1534,bipolar.size()/2.0,1.77,fAsymmetry);
    TH1D * bipol = new TH1D("bipol","induction input shape",signalSize,0,signalSize);
    TH1D * rampc = new TH1D("ramp","collection input shape",signalSize,0,signalSize);
    for(int i = 0; i < signalSize; i++) {
      ramp[i]=rampf.Eval(i);
      rampc->Fill(i,ramp[i]);
      bipolar[i]=bipolf.Eval(i);
      bipol->Fill(i,bipolar[i]);
    }
    TF1 * noise = new TF1("noise",fNoiseFunc.c_str(),0,size+1);
    noise->SetParameters(&fNParams[0]);
    TH1D * noiseHist = new TH1D("noise","noise",size+1,0,size+1);
    for(int i = 0; i < size+1; i++) noiseHist->Fill(i,noise->Eval(i));
    TF1 * indFilter = new TF1("indFilter",fRIFFunc.c_str(),0,size+1);
    indFilter->SetParameters(&fRIFParams[0]); 
    TF1 * colFilter = new TF1("colFilter",fRCFFunc.c_str(),0,size+1);
    colFilter->SetParameters(&fRCFParams[0]);
    FFT->ShiftData(bipolar,(double)bipolar.size()/2.0);
    FFT->ShiftData(ramp,3.57);
    std::vector<double> electronics(signalSize);
    //integrated electronics shaping function with no RC tail decay
    TF1 * expt = new TF1("expt",fElectFunc.c_str(),0,signalSize+1); 
    TH1D * elect = new TH1D("electronics","electronics shaping",signalSize,0,signalSize);
    expt->SetParameters(&fEParams[0]); 
    for(int i = 0; i < signalSize; i++) {
      electronics[i] = expt->Eval(i);
      elect->Fill(i,electronics[i]);
    }
    TFile * out = new TFile("shape.root","NEW");
    bipol->Write();
    rampc->Write();
    elect->Write();
    noiseHist->Write();
    TDirectory *real = out->mkdir("real");
    real->AddDirectory();
    real->cd();
    
    TH1D* decayHistNew = new TH1D("decayHist","RC decay Constants",nChannels,0,nChannels);
    std::vector<double> induction(signalSize);
    std::vector<double> collection(signalSize);
    std::vector<double> delta(signalSize);
    delta[signalSize-1]=1.0;
    delta[0]=1.0;
    induction[0]=0;
    for(double i = 0; i < signalSize; i++) {
      induction[i]=electronics[i];
    }    
    FFT->Convolute(induction,bipolar);
    TH1D * indShape = new TH1D("indShape","Induction Shape",signalSize,0,signalSize);
  
    TH1D * colShape = new TH1D("colShape","Collection Shape",signalSize,0,signalSize);
    for(int i = 0; i < signalSize; i++) { 
      indShape->Fill(i,induction[i]);
    }
    TH1D * indPow = new TH1D("indPow","Induction power spectrum",size+1,0,size+1);
    TH1D * colPow = new TH1D("colPow","Collection power spectrum",size+1,0,size+1);
    TH1D * indFil = new TH1D("indFil","Induction Filter",size+1,0,size+1);
    TH1D * colFil = new TH1D("colFil","Collection Filter",size+1,0,size+1);
    for(bin = 0; bin < size+1; bin++) {
      indFil->Fill(bin,indFilter->Eval(bin));
      colFil->Fill(bin,colFilter->Eval(bin));
    }

    TH2D * RespRe = new TH2D("RespRe",
			     " Response Functions - Real Part",
			     nRanges,&ranges[0],size+1, 0, size+1);
    TH2D * RespIm = new TH2D("RespIm",
			     "ResponseFunctions - Imaginary Part",   
			     nRanges,&ranges[0],size+1, 0, size+1);

    std::vector<double> dataTemp;
    for(int i = 0; i < nRanges-1; i++) {
      if(i < iRange) {
	dataTemp = induction;
      }
      else {
	dataTemp[0]=0;
	for(double j = 1.0; j < signalSize; j++){ 
	  dataTemp[j]=electronics[j];
	}
      }
      FFT->AlignedSum(dataTemp,delta,false);
      FFT->DoFFT(dataTemp,freqSig);
      for(bin = 0; bin < 0.5*signalSize+1; bin++) {
	if(i<iRange) {
	  kernBin=indFilter->Eval(bin)/freqSig[bin];
	  indPow->Fill(bin,freqSig[bin].Rho());
	}
	else {
	  kernBin=colFilter->Eval(bin)/freqSig[bin];
	  colPow->Fill(bin,freqSig[bin].Rho());
	}
	RespRe->Fill(ranges[i],bin,kernBin.Re());     
	RespIm->Fill(ranges[i],bin,kernBin.Im());
      }
    }
    indShape->Write();
    colShape->Write();
    indPow->Write();
    colPow->Write();
    indFil->Write();
    colFil->Write();
    RespRe->Write();
    RespIm->Write();
    decayHistNew->Write();
    real->cd("../");
    //sim section
    TDirectory * sim = out->mkdir("sim");
    sim->cd();
    TF1* indFilterS = new TF1("indFilterS",fSIFFunc.c_str(),0,size+1);
    indFilterS->SetParameters(&fSIFParams[0]); 
    TF1 * colFilterS = new TF1("colFilterS",fSCFFunc.c_str(),0,size+1);
    colFilterS->SetParameters(&fSCFParams[0]); 
    TH2D * RespRe1 = new TH2D("RespRe","Response Functions - Real Part",nRanges,&ranges[0],signalSize/2+1,0,signalSize/2+1);
    TH2D* RespIm1 = new TH2D("RespIm","Response Functions - Imaginary Part",nRanges,&ranges[0],signalSize/2+1,0,signalSize/2+1);
    TH1D * indFil1 = new TH1D("indFil","Induction Filter",size+1,0,size+1);
    TH1D * colFil1 = new TH1D("colFil","Collection Filter",size+1,0,size+1);
    for(int i = 0; i < signalSize; i++) {
      induction[i]=electronics[i];
      collection[i]=electronics[i];
    } 
    FFT->Convolute(induction,bipolar);
    FFT->Convolute(collection,ramp);
    FFT->DoFFT(induction,freqSig);
    std::vector<double> indWFilter (signalSize);
    std::vector<double> colWFilter (signalSize);
    for(int i = 0 ; i < size+1; i++) indWFilter[i]=indFilterS->Eval(i);//*freqSig[i].Rho2()/(freqSig[i].Rho2()+noise->Eval(i)*noise->Eval(i));
    TH1D * indPow1 = new TH1D("indPow","Induction power spectrum",size+1,0,size+1);
    TH1D * colPow1 = new TH1D("colPow","Collection power spectrum",size+1,0,size+1);
    for(int i = 0; i <signalSize/2+1; i++) {
      indPow1->Fill(i,freqSig[i].Rho());
      kernBin=indWFilter[i]/freqSig[i];
      if(i==0 || i>1100) kernBin=0;
      indFil1->Fill(i,indWFilter[i]);
      RespRe1->Fill(0.0,i,kernBin.Re()); 
      RespRe1->Fill(iRange-1,i,kernBin.Re());
      RespIm1->Fill(0.0,i,kernBin.Im());
      RespRe1->Fill(iRange-1,i,kernBin.Im());
    }
    FFT->DoFFT(collection,freqSig);
    for(int i = 0 ; i < signalSize; i++) colWFilter[i]=colFilterS->Eval(i);//*exp(-0.5*pow((double)i/135.0,4.0));//*freqSig[i].Rho2()/(freqSig[i].Rho2()+0.5*noise->Eval(i)*noise->Eval(i));
    for(int i = 0; i <signalSize/2+1; i++) {
      colPow1->Fill(i,freqSig[i].Rho());
      kernBin=colWFilter[i]/freqSig[i];
      colFil1->Fill(i,colWFilter[i]);
      RespRe1->Fill(iRange+1,i,kernBin.Re()); 
      RespIm1->Fill(iRange+1,i,kernBin.Im());
    }
    indPow1->Write();
    colPow1->Write();
    indFil1->Write();
    colFil1->Write();
    RespRe1->Write();
    RespIm1->Write();
    decayHistNew->Write();
    out->cd("../");
    TH2D * shapes = new TH2D("shapes","Convolved drift and electronics shapes",nPlanes,0,nPlanes,signalSize,0,signalSize);
    for(int i = 0; i < signalSize; i++) {
      shapes->Fill(0.0,i,fIndSimScale*induction[i]);
      shapes->Fill(1.0,i,fIndSimScale*induction[i]);
      shapes->Fill(2.0,i,fColSimScale*collection[i]);
    }
    shapes->Write();
    out->Write();
     
  }
}


namespace uboone{

  DEFINE_ART_MODULE(MicrobooneResp)
  
} // end namespace uboone

#endif // MICROBOONERESP_H
