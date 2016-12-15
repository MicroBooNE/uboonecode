#include "uboone/CRT/CRTData.hh"
#include "uboone/CRT/CRTDetSim.hh"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

#include <vector>
#include <memory>
#include <sstream>
#include <string>
#include <math.h>

namespace crt{

  CRTDetSim::CRTDetSim(const fhicl::ParameterSet& pSet)
  {
    this->reconfigure(pSet);
    produces< std::vector<CRTData> >();
  }

  CRTDetSim::~CRTDetSim()
  {

  }

  void CRTDetSim::reconfigure(fhicl::ParameterSet const & pSet) {
    fTDelayNorm = pSet.get<double>("TDelayNorm");
    fTDelayShift = pSet.get<double>("TDelayShift");
    fTDelaySigma = pSet.get<double>("TDelaySigma");
    fTDelayOffset = pSet.get<double>("TDelayOffset");
    fTDelayRMSGausNorm = pSet.get<double>("TDelayRMSGausNorm");
    fTDelayRMSGausShift = pSet.get<double>("TDelayRMSGausShift");
    fTDelayRMSGausSigma = pSet.get<double>("TDelayRMSGausSigma");
    fTDelayRMSExpNorm = pSet.get<double>("TDelayRMSExpNorm");
    fTDelayRMSExpShift = pSet.get<double>("TDelayRMSExpShift");
    fTDelayRMSExpScale = pSet.get<double>("TDelayRMSExpScale");
    fPropDelay = pSet.get<double>("PropDelay");
    fPropDelayError = pSet.get<double>("fPropDelayError");
    fTResInterpolator = pSet.get<double>("TResInterpolator");
    fNpeScaleNorm = pSet.get<double>("NpeScaleNorm");
    fNpeScaleShift = pSet.get<double>("NpeScaleShift");
    fQ0 = pSet.get<double>("Q0");
    fQPed = pSet.get<double>("QPed");
    fQSlope = pSet.get<double>("QSlope");
    fQRMS = pSet.get<double>("QRMS");
    fAbsLenEff = pSet.get<double>("AbsLenEff");
    fProducerName = pSet.get<std::string>("ProducerName", "largeant");
  }

  double CRTDetSim::getChannelTriggerTicks(CLHEP::HepRandomEngine* engine,
                                           detinfo::ElecClock& clock,
                                           float t0, float npeMean, float r) {
    // Hit timing, with smearing and NPE dependence
    double tDelayMean = \
      fTDelayNorm *
        exp(-0.5 * pow((npeMean - fTDelayShift) / fTDelaySigma, 2)) +
      fTDelayOffset;

    double tDelayRMS = \
      fTDelayRMSGausNorm *
        exp(-pow(npeMean - fTDelayRMSGausShift, 2) / fTDelayRMSGausSigma) +
      fTDelayRMSExpNorm *
        exp(-(npeMean - fTDelayRMSExpShift) / fTDelayRMSExpScale);

    double tDelay = CLHEP::RandGauss::shoot(engine, tDelayMean, tDelayRMS);

    // Time resolution of the interpolator
    tDelay += CLHEP::RandGauss::shoot(engine, 0, fTResInterpolator);

    // Propagation time
    double tProp = CLHEP::RandGauss::shoot(fPropDelay, fPropDelayError) * r;

    double t = t0 + tProp + tDelay;

    // Get clock ticks
    clock.SetTime(t);
    return clock.Ticks();
  }

  void CRTDetSim::produce(art::Event& evt)
  {
    std::unique_ptr<std::vector<crt::CRTData> > crtHits(
        new std::vector<crt::CRTData>);

    art::ServiceHandle<geo::AuxDetGeometry> geoService;

    art::ServiceHandle<detinfo::DetectorClocksService> detClocks;

    detinfo::ElecClock trigClock = detClocks->provider()->TriggerClock();

    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine* engine = &rng->getEngine("crt");

    // Handle for (truth) AuxDetSimChannels
    art::Handle<std::vector<sim::AuxDetSimChannel> > channels;
    evt.getByLabel(fProducerName, channels);

    // I hate that you people made me do this.
    const geo::AuxDetGeometry* geometry = &*geoService;
    const geo::AuxDetGeometryCore* geoServiceProvider = geometry->GetProviderPtr();

    // Loop through truth AD channels
    for (auto& adsc : *channels) {

      const geo::AuxDetGeo& adGeo = geoServiceProvider->AuxDet(adsc.AuxDetID());

      const geo::AuxDetSensitiveGeo& adsGeo = adGeo.SensitiveVolume(adsc.AuxDetSensitiveID());

      // Simulate the CRT response for each hit
      for (auto ide : adsc.AuxDetIDEs()) {
        // Get the hit position in strip's local coordinates
        double x = (ide.entryX + ide.exitX) / 2;
        double y = (ide.entryY + ide.exitY) / 2;
        double z = (ide.entryZ + ide.exitZ) / 2;
        double world[3] = {x, y, z};
        double svHitPosLocal[3];
        adsGeo.WorldToLocal(world, svHitPosLocal);

        // Distance to the readout end ("outside") depends on module position
        // FIXME: FOR NOW ASSUME ALL THE SAME DIRECTION
        double distToReadout = abs(-adsGeo.HalfHeight() - svHitPosLocal[1]);

        // The expected number of PE
        double qr = ide.energyDeposited / fQ0;  // Scale linearly with charge
        double npeExpected = (fNpeScaleNorm / pow(distToReadout - fNpeScaleShift, 2) * qr);

        // Put PE on channels weighted by distance
        double d0 = abs(-adsGeo.HalfWidth1() - svHitPosLocal[0]);  // L
        double d1 = abs( adsGeo.HalfWidth1() - svHitPosLocal[0]);  // R
        double abs0 = exp(-d0 / fAbsLenEff);
        double abs1 = exp(-d1 / fAbsLenEff);
        double npeExp0 = npeExpected * abs0 / (abs0 + abs1);
        double npeExp1 = npeExpected * abs1 / (abs0 + abs1);

        // Observed PE
        long npe0 = CLHEP::RandPoisson::shoot(engine, npeExp0);
        long npe1 = CLHEP::RandPoisson::shoot(engine, npeExp1);

        // Time relative to trigger
        double tTrue = (ide.entryT + ide.exitT) / 2;
        uint32_t t0 = getChannelTriggerTicks(engine, trigClock, tTrue, npe0, distToReadout);
        uint32_t t1 = getChannelTriggerTicks(engine, trigClock, tTrue, npe1, distToReadout);

        // Time relative to PPS: Random for now
        uint32_t ppsTicks = CLHEP::RandFlat::shootInt(engine, trigClock.Frequency() * 1e6);

        // SiPM and ADC response: Npe to ADC counts
        short q0 = CLHEP::RandGauss::shoot(engine, fQPed + fQSlope * npe0, fQRMS * npe0);
        short q1 = CLHEP::RandGauss::shoot(engine, fQPed + fQSlope * npe1, fQRMS * npe1);

        // Adjacent channels on a strip are numbered sequentially
        uint32_t moduleID = adsc.AuxDetID();
        uint32_t stripID = adsc.AuxDetSensitiveID();
        uint32_t channel0ID = 32 * moduleID + 2 * stripID + 0;
        uint32_t channel1ID = 32 * moduleID + 2 * stripID + 1;

        // Write AuxDetDigit for each channel
        crtHits->push_back(CRTData(channel0ID, t0, ppsTicks, q0));
        crtHits->push_back(CRTData(channel1ID, t1, ppsTicks, q1));
      }
    }

    evt.put(std::move(crtHits));
  }

  DEFINE_ART_MODULE(CRTDetSim)
}
