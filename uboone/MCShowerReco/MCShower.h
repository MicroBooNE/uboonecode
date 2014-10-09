#ifndef MCSHOWER_H
#define MCSHOWER_H

#include <vector>
//#include <TMath.h>
#include "MCBase/MCLimits.h"

namespace sim {

  class MCShower {

  public:

    MCShower()
    { Clear(); }
    ~MCShower(){}

    //---- Mother info ----//
    int momPdgCode;                 ///< mother PDG code
    unsigned int momTrackId;        ///< mother G4 Track ID
    std::string momProcess;         ///< mother's creation process
    std::vector<double> vtxMother;  ///< mother position 4-vector @ generation
    std::vector<double> momMother;  ///< mother momentum 4-vector @ generation
    /// mother 3D angle phi (along shower angle definition, not ordinary coord. system)
    double phiMother;
    /// mother 3D angle theta (along shower angle definition, not ordinary coord. system)
    double thetaMother;
    /// mother 2D angle on U-plane
    double uAngleMother;
    /// mother 2D angle on V-plane
    double vAngleMother;
    /// mother 2D angle on W-plane
    double wAngleMother;

    //---- Daughter info ----//
    std::vector<unsigned int> daughterTrackId; ///< Daughters' track ID
    std::vector<double> vtxDaughter;           ///< Daughters' first energy deposition vtx
    std::vector<double> momDaughter;           ///< Daughters' deposit sum momentum 4-vector
    //std::vector<float> vtxDaughter;           ///< Daughters' first energy deposition vtx
    //std::vector<float> momDaughter;           ///< Daughters' deposit sum momentum 4-vector
    /// daughter 3D angle phi (along shower angle definition, not ordinary coord. system)
    float phiDaughter;
    /// daughter 3D angle theta (along shower angle definition, not ordinary coord. system)
    float thetaDaughter;
    /// daughter 2D angle on U-plane
    float uAngleDaughter;
    /// daughter 2D angle on V-plane
    float vAngleDaughter;
    /// daughter 2D angle on W-plane
    float wAngleDaughter;

    //---- Charge per plane ----//
    float qU; ///< Charge deposit on U plane
    float qV; ///< Charge deposit on V plane
    float qW; ///< Charge deposit on W plane

    /// Energy deposit point (x,y,z,dE) by daughters
    std::vector<std::vector<float> > vtxEdep;

    //#ifndef __GCCXML__
    
    void Clear() {

      momPdgCode = ::sim::kINVALID_INT;
      momTrackId = ::sim::kINVALID_UINT;      
      momProcess = "";
      vtxMother.clear();
      vtxMother.resize(4,0);
      momMother.clear();
      momMother.resize(4,0);
      //thetaMother = phiMother = TMath::Pi()*4;
      //uAngleMother = vAngleMother = wAngleMother = TMath::Pi()*4;
      thetaMother = phiMother = ::sim::kINVALID_FLOAT;
      uAngleMother = vAngleMother = wAngleMother = ::sim::kINVALID_FLOAT;

      daughterTrackId.clear();
      vtxDaughter.clear();
      vtxDaughter.resize(4,0);
      momDaughter.clear();
      momDaughter.resize(4,0);
      //phiDaughter = thetaDaughter = TMath::Pi()*4;
      //uAngleDaughter = vAngleDaughter = wAngleDaughter = TMath::Pi()*4;
      phiDaughter = thetaDaughter = ::sim::kINVALID_FLOAT;
      uAngleDaughter = vAngleDaughter = wAngleDaughter = ::sim::kINVALID_FLOAT;

      qU = qV = qW = 0;
      vtxEdep.clear();

    }
    //#endif
  };

}

#endif
