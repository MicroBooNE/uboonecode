//
// Name:  PropTest_module.cc
//
// Purpose: PropTest module.  Test propagators.
//
// Created:  12-Apr-2012  H. Greenlee

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <vector>

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"

#include "RecoObjects/KETrack.h"
#include "RecoObjects/SurfYZLine.h"
#include "RecoObjects/SurfYZPlane.h"
#include "RecoObjects/SurfXYZPlane.h"
#include "RecoObjects/PropAny.h"

#include "TMath.h"

namespace trkf
{
  class PropTest : public art::EDAnalyzer
  {
  public:

    // Constructor, destructor.

    explicit PropTest(fhicl::ParameterSet const& pset);
    ~PropTest();

    // Overrides.

    void beginJob();
    void analyze(const art::Event& evt);

  };

  DEFINE_ART_MODULE(PropTest)

  PropTest::PropTest(const fhicl::ParameterSet& pset)
  : EDAnalyzer(pset)
  {}

  void PropTest::beginJob()
  {
    // Make sure assert is enabled.

    bool assert_flag = false;
    assert((assert_flag = true, assert_flag));
    if ( ! assert_flag ) {
      std::cerr << "Assert is disabled" << std::endl;
      abort();
    }

    // Make a PropAny propagate to test.

    const trkf::Propagator* prop = new trkf::PropAny(10., true);

    // Make some random surfaces.
    // Also make initial tracks.

    std::vector<std::shared_ptr<const trkf::Surface> > surfaces;
    std::vector<trkf::KETrack> tracks;

    int nsurf = 30;
    for(int isurf = 0; isurf < nsurf; ++isurf) {

      // Make random surface.

      std::shared_ptr<const trkf::Surface> psurf;
      if(isurf < 10) {
	double x0 = 10.*double(rand()) / double(RAND_MAX) - 5.;  // (-5,5)
	double y0 = 10.*double(rand()) / double(RAND_MAX) - 5.;  // (-5,5)
	double z0 = 10.*double(rand()) / double(RAND_MAX) - 5.;  // (-5,5)
	double phi = TMath::TwoPi() * double(rand()) / double(RAND_MAX) - TMath::Pi();  // (-pi,pi)
	psurf = std::shared_ptr<const trkf::Surface>(new trkf::SurfYZLine(x0, y0, z0, phi));
	surfaces.push_back(psurf);
      }
      else if(isurf < 20) {
	double x0 = 10.*double(rand()) / double(RAND_MAX) - 5.;  // (-5,5)
	double y0 = 10.*double(rand()) / double(RAND_MAX) - 5.;  // (-5,5)
	double z0 = 10.*double(rand()) / double(RAND_MAX) - 5.;  // (-5,5)
	double phi = TMath::TwoPi() * double(rand()) / double(RAND_MAX) - TMath::Pi();  // (-pi,pi)
	psurf = std::shared_ptr<const trkf::Surface>(new trkf::SurfYZPlane(x0, y0, z0, phi));
	surfaces.push_back(psurf);
      }
      else {
	double x0 = 10.*double(rand()) / double(RAND_MAX) - 5.;  // (-5,5)
	double y0 = 10.*double(rand()) / double(RAND_MAX) - 5.;  // (-5,5)
	double z0 = 10.*double(rand()) / double(RAND_MAX) - 5.;  // (-5,5)
	double theta = std::acos(2. * double(rand()) / double(RAND_MAX) - 1.);  // (0, pi)
	double phi = TMath::TwoPi() * double(rand()) / double(RAND_MAX) - TMath::Pi();  // (-pi,pi)
	psurf = std::shared_ptr<const trkf::Surface>(new trkf::SurfXYZPlane(x0, y0, z0, theta, phi));
	surfaces.push_back(psurf);
      }

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

      // Make initial error matrix (diagonal, not random).

      trkf::TrackError err(5);
      err.clear();
      err(0,0) = 0.1;
      err(1,1) = 0.1;
      err(2,2) = 0.01;
      err(3,3) = 0.01;
      err(4,4) = 0.01;

      // Make random track direction.

      trkf::Surface::TrackDirection dir = trkf::Surface::FORWARD;
      if(rand() % 2 == 0)
	dir = trkf::Surface::BACKWARD;

      // Make KETrack.

      tracks.push_back(trkf::KETrack(psurf, vec, err, dir));    
    }

    // Loop over initial surface/track.

    for(int isurf = 0; isurf < nsurf; ++isurf) {

      std::cout << "\nInitial track " << isurf << std::endl;

      // Get initial surface and track object.

      const std::shared_ptr<const trkf::Surface>& psurf1 = surfaces[isurf];
      const trkf::KETrack& trk1 = tracks[isurf];
      const trkf::TrackVector& vec1 = trk1.getVector();
      const trkf::TrackError& err1 = trk1.getError();

      // Get initial position and momentum vector.

      double xyz1[3];
      double mom1[3];
      trk1.getPosition(xyz1);
      trk1.getMomentum(mom1);
      std::cout << "Initial position = "
		<< xyz1[0] << ", " << xyz1[1] << ", " << xyz1[2] << std::endl;
      std::cout << "Initial momentum = "
		<< mom1[0] << ", " << mom1[1] << ", " << mom1[2] << std::endl;

      // Loop over destination surface.

      for(int jsurf = 0; jsurf < nsurf; ++jsurf) {

	std::cout << "Destination " << jsurf << std::endl;

	// Get destination surface.

	const std::shared_ptr<const trkf::Surface>& psurf2 = surfaces[jsurf];

	// Make a copy of the initial track and propagate it to the 
	// destination surface.

	trkf::KETrack trk2 = trk1;
	boost::optional<double> dist12 = prop->err_prop(trk2, psurf2,
							trkf::Propagator::UNKNOWN, false);
	if(!!dist12)
	  std::cout << "Propagation distance = " << *dist12 << std::endl;
	else
	  assert(false);

	// Get final position and momentum vector.

	double xyz2[3];
	double mom2[3];
	trk2.getPosition(xyz2);
	trk2.getMomentum(mom2);
	std::cout << "Final position = "
		  << xyz2[0] << ", " << xyz2[1] << ", " << xyz2[2] << std::endl;
	std::cout << "Final momentum = "
		  << mom2[0] << ", " << mom2[1] << ", " << mom2[2] << std::endl;

	// Calculate displacement vector.

	double dx[3];
	dx[0] = xyz2[0] - xyz1[0];
	dx[1] = xyz2[1] - xyz1[1];
	dx[2] = xyz2[2] - xyz1[2];

	// Check the sign of the propagation distance.
	// The sign of the propatation distance should be the same as 
	// the sign of the dot product of the initial momentum and the
	// displacement vector.

	double pdotdx = mom1[0]*dx[0] + mom1[1]*dx[1] + mom1[2]*dx[2];
	assert(pdotdx * (*dist12) >= -1.e-10);

	// Check propagation distance against displacement vector.

	double dist = std::sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
	std::cout << "Displacement = " << dist << std::endl;
	assert(std::abs(dist - std::abs(*dist12)) <= 1.e-8);

	// Make sure final momentum is parallel to initial momentum.

	double m0 = mom1[1]*mom2[2] - mom1[2]*mom2[1];
	double m1 = mom1[2]*mom2[0] - mom1[0]*mom2[2];
	double m2 = mom1[0]*mom2[1] - mom1[1]*mom2[0];
	std::cout << "mom1 x mom2 = " << m0 << ", " << m1 << ", " << m2 << std::endl;
	assert(std::abs(m0) <= 1.e-8);
	assert(std::abs(m1) <= 1.e-8);
	assert(std::abs(m2) <= 1.e-8);

	// Make sure displacement vector is parallel to initial momentum.

	double d0 = mom1[1]*dx[2] - mom1[2]*dx[1];
	double d1 = mom1[2]*dx[0] - mom1[0]*dx[2];
	double d2 = mom1[0]*dx[1] - mom1[1]*dx[0];
	std::cout << "mom1 x dx = " << d0 << ", " << d1 << ", " << d2 << std::endl;
	assert(std::abs(d0) <= 1.e-8);
	assert(std::abs(d1) <= 1.e-8);
	assert(std::abs(d2) <= 1.e-8);

	// Test propagation matrix by numerical partial derivatives.

	trkf::TrackMatrix pm(vec1.size(), vec1.size());
	trkf::KTrack trk10 = trk1;
	boost::optional<double> stat =
	  prop->vec_prop(trk10, psurf2, trkf::Propagator::UNKNOWN, false, &pm, 0);
	assert(!!stat);

	double small = 1.e-5;

	for(unsigned int i = 0; i < vec1.size(); ++i) {
	  for(unsigned int j = 0; j < vec1.size(); ++j) {

	    // Calculate d(vec[i])/d(vec[j])

	    trkf::KTrack trk1a = trk1;
	    trkf::TrackVector vec1a = vec1;
	    vec1a(j) = vec1(j) - small;
	    trk1a.setVector(vec1a);
	    boost::optional<double> stata = prop->vec_prop(trk1a, psurf2,
							   trkf::Propagator::UNKNOWN,
							   false, 0, 0);
	    assert(!!stata);

	    trkf::KTrack trk1b = trk1;
	    trkf::TrackVector vec1b = vec1;
	    vec1b(j) = vec1(j) + small;
	    trk1b.setVector(vec1b);
	    boost::optional<double> statb = prop->vec_prop(trk1b, psurf2,
							   trkf::Propagator::UNKNOWN,
							   false, 0, 0);
	    assert(!!statb);

	    // Compare numerical and analytic partial derivative.

	    double dij = (trk1b.getVector()(i) - trk1a.getVector()(i)) / (2.*small);
	    std::cout << "(" << i << "," << j << "): " << dij << ", " << pm(i,j) << std::endl;
	    assert(std::abs(dij - pm(i,j)) <= 1.e-4*std::max(std::abs(dij), 1.));
	  }
	}

	// Now propagate back to the original surface.

	boost::optional<double> dist21 = prop->err_prop(trk2, psurf1, trkf::Propagator::UNKNOWN,
							false);
	assert(!!dist21);
	assert(std::abs(dist - std::abs(*dist21)) <= 1.e-10 * std::max(1., dist));

	// Check that state vector and error matrix returned to the original.
	// This will test that the forward and backward propagation matrices
	// are inverses.

	const trkf::TrackVector& vec2 = trk2.getVector();
	const trkf::TrackError& err2 = trk2.getError();
	int n = vec1.size();
	for(int i=0; i<n; ++i) {
	  assert(std::abs(vec1(i) - vec2(i)) <= 1.e-8);
	  for(int j=0; j<n; ++j) {
	    assert(std::abs(err1(i,j) - err2(i,j)) <= 1.e-4*std::max(std::abs(err1(i,j)), 1.));
	  }
	}
      }
    }

    // Test noise matrix calculation.  Test that error matrix
    // resulting from propagation of KETrack with zero initial error
    // is the same for propagation in a single step as compared to
    // propagation in many steps.

    // Make initial and final surfaces on planes of constant z.

    double z0 = 0.;
    double z1 = 100.;
    std::shared_ptr<const trkf::SurfYZPlane> psurf0(new trkf::SurfYZPlane(0., 0., z0, 0.));
    std::shared_ptr<const trkf::SurfYZPlane> psurf1(new trkf::SurfYZPlane(0., 0., z1, 0.));

    // Make initial KETrack (muon) with zero error.

    trkf::TrackVector vec(5);
    vec(0) = 0.;
    vec(1) = 0.;
    vec(2) = 1.;
    vec(3) = 2.;
    vec(4) = 1.;
    trkf::TrackError err(5);
    err.clear();
    trkf::KETrack tre0(psurf0, vec, err, trkf::Surface::FORWARD, 13);

    std::cout << "\nNoise matrix calculation test." << std::endl;
    std::cout << "\nInitial track:" << std::endl;
    std::cout << tre0 << std::endl;

    // Propagate track to destination surface in one step.

    trkf::KETrack tre1(tre0);
    prop->noise_prop(tre1, psurf1, trkf::Propagator::UNKNOWN, false);
    std::cout << "\nSingle step final track:" << std::endl;
    std::cout << tre1 << std::endl;

    // Propagate track to destination using multiple steps.

    trkf::KETrack tre2(tre0);
    int nprop = 100;
    for(int iprop=1; iprop <= nprop; ++iprop) {
      double z = z0 + (z1-z0) * double(iprop) / double(nprop);
      std::shared_ptr<const trkf::SurfYZPlane> psurf(new trkf::SurfYZPlane(0., 0., z, 0.));
      prop->noise_prop(tre2, psurf, trkf::Propagator::UNKNOWN, false);
      if(iprop == nprop) {
	std::cout << "\nStep " << iprop << " track: " << std::endl;
	std::cout << tre2 << std::endl;
      }
    }

    // Make sure two tracks are the same.

    trkf::TrackVector vec1 = tre1.getVector();
    trkf::TrackError err1 = tre1.getError();
    trkf::TrackVector vec2 = tre2.getVector();
    trkf::TrackError err2 = tre2.getError();
    for(unsigned int i=0; i<vec1.size(); ++i) {
      assert(std::abs(vec1(i) - vec2(i)) < 1.e-6);
      for(unsigned int j=0; j<=i; ++j)
	assert(std::abs(err1(i,j) - err2(i,j)) <= 1.e-6);
    }
    // Done (success).

    std::cout << "PropTest: All tests passed." << std::endl;
  }

  PropTest::~PropTest()
  {}

  void PropTest::analyze(const art::Event& /* evt */)
  {}
}
