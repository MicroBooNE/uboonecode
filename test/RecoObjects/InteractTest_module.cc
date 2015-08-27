//
// Name:  InteractTest_module.cc
//
// Purpose: InteractTest module.  Test interactors.
//
// Created:  12-Apr-2012  H. Greenlee

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <vector>

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"

#include "RecoObjects/KTrack.h"
#include "RecoObjects/SurfXYZPlane.h"
#include "RecoObjects/InteractPlane.h"
#include "RecoObjects/InteractGeneral.h"

#include "TMath.h"

namespace trkf
{
  class InteractTest : public art::EDAnalyzer
  {
  public:

    // Constructor, destructor.

    explicit InteractTest(fhicl::ParameterSet const& pset);
    ~InteractTest();

    // Overrides.

    void beginJob();
    void analyze(const art::Event& evt);

  };

  DEFINE_ART_MODULE(InteractTest)

  InteractTest::InteractTest(const fhicl::ParameterSet& pset)
  : EDAnalyzer(pset)
  {}

  void InteractTest::beginJob()
  {
    // Make sure assert is enabled.

    bool assert_flag = false;
    assert((assert_flag = true, assert_flag));
    if ( ! assert_flag ) {
      std::cerr << "Assert is disabled" << std::endl;
      abort();
    }

    // Make InteractPlane and InteractGenreal interactors to compare.

    const trkf::Interactor* int1 = new trkf::InteractPlane(10.);
    const trkf::Interactor* int2 = new trkf::InteractGeneral(10.);

    // Make some random surfaces and tracks.

    std::vector<trkf::KTrack> tracks;

    int ntrk = 10;
    for(int itrk = 0; itrk < ntrk; ++itrk) {

      // Make random surface.

      std::shared_ptr<const trkf::Surface> psurf;
      double x0 = 10.*double(rand()) / double(RAND_MAX) - 5.;  // (-5,5)
      double y0 = 10.*double(rand()) / double(RAND_MAX) - 5.;  // (-5,5)
      double z0 = 10.*double(rand()) / double(RAND_MAX) - 5.;  // (-5,5)
      double theta = std::acos(2. * double(rand()) / double(RAND_MAX) - 1.);  // (0, pi)
      double phi = TMath::TwoPi() * double(rand()) / double(RAND_MAX) - TMath::Pi();  // (-pi,pi)
      psurf = std::shared_ptr<const trkf::Surface>(new trkf::SurfXYZPlane(x0, y0, z0, theta, phi));

      // Make random track vector.

      double u = 100.*double(rand()) / double(RAND_MAX);  // (0,100)
      double v = 100.*double(rand()) / double(RAND_MAX) - 50.;  // (-50, 50)
      double dudw = 2.*double(rand()) / double(RAND_MAX) - 1.;  // (-1, 1)
      double dvdw = 2.*double(rand()) / double(RAND_MAX) - 1.;  // (-1, 1)
      double pinv = 0.9*double(rand()) / double(RAND_MAX) + 0.1;  // (0.1, 1.0)
      trkf::TrackVector vec(5);
      vec(0) = u;
      vec(1) = v;
      vec(2) = dudw;
      vec(3) = dvdw;
      vec(4) = pinv;

      // Make random track direction.

      trkf::Surface::TrackDirection dir = trkf::Surface::FORWARD;
      if(rand() % 2 == 0)
	dir = trkf::Surface::BACKWARD;

      // Make KTrack.

      tracks.emplace_back(psurf, vec, dir, 13);
    }

    // Loop over surface/track.

    for(size_t itrk = 0; itrk < tracks.size(); ++itrk) {

      std::cout << "\nInitial track " << itrk << std::endl;

      // Get track object.

      const trkf::KTrack& trk = tracks[itrk];

      // Choose a random distance, including negative distance.

      double s = 200.*double(rand()) / double(RAND_MAX) - 100.;

      // Calculate noise matrix using both interactors.

      TrackError noise1(5);
      TrackError noise2(5);
      int1->noise(trk, s, noise1); 
      int2->noise(trk, s, noise2); 

      // Check that the noise matrices are the same, within tolerance.

      for(size_t i=0; i<noise1.size1(); ++i) {
	for(size_t j=0; j<=i; ++j) {
	  std::cout << i << ' ' << j << ' ' << noise1(i,j) << ' ' << noise2(i,j) << std::endl;
	  assert(std::abs(noise1(i,j) - noise2(i,j)) <= 1.e-6*std::max(std::abs(noise1(i,j)), 1.));
	}
      }
    }

    // Done (success).

    std::cout << "InteractTest: All tests passed." << std::endl;
  }

  InteractTest::~InteractTest()
  {}

  void InteractTest::analyze(const art::Event& /* evt */)
  {}
}
