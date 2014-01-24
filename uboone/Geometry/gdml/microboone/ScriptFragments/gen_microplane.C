#include "Riostream.h"

void gen_microplane(){

  ofstream of("micro-plane.gdml");

  int nsw = 1;
  int ncw = 1;
  //define the solids
  of << "<solids>" << endl;

  //wires on either end of the tpc
  for(int i = 0; i < nsw; ++i){
    of << "<tube name=" << TString::Format("\"TPCWire%d\" ", i)             << endl;
    of << "  rmax=\"0.5*kTPCWireThickness\" "                               << endl;
    of << "  z=" << TString::Format("\"kTPCWireXPitch*(%d+1)/kCos30\" ", i) << endl;
    of << "  deltaphi=\"2*kPi\" "                                           << endl;
    of << "  aunit=\"rad\" "                                                << endl;
    of << "  lunit=\"cm\"/> "                                               << endl;
  }
  
  //the middle wires
  of << "<tube name=\"TPCWireCommon\""                                    << endl;
  of << "  rmax=\"0.5*kTPCWireThickness\" "                               << endl;
  of << "  z=\"kTPCWidth/kCos60\" "                                       << endl;
  of << "  deltaphi=\"2*kPi\" "                                           << endl;
  of << "  aunit=\"rad\" "                                                << endl;
  of << "  lunit=\"cm\"/> "                                               << endl;

  //the tpc wire plane
  of << "<box name=\"TPCPlane\" "                                         << endl;
  of << "  x=\"kTPCWirePlaneThickness\" "                                 << endl;
  of << "  y=\"kTPCWidth\" "                                              << endl;
  of << "  z=\"kTPCLength\" "                                             << endl;
  of << "  lunit=\"cm\"/> "                                               << endl;

  //end solids
  of << "</solids>" << endl;

  //now defined the structures
  of << "<structure>" << endl;
  
  //the wires at either end of the plane
  for(int i = 0; i < nsw; ++i){
    of << "  <volume name=" << TString::Format("\"volTPCWire%d\"> ", i)  << endl;
    of << "    <materialref ref=\"Titanium\"/>"                          << endl;
    of << "    <solidref ref=" << TString::Format("\"TPCWire%d\"/> ", i) << endl;
    of << "  </volume> "                                                 << endl;
  }

//   //the wires in the middle of the plane
  of << "  <volume name=\"volTPCWireCommon\"> "      << endl;
  of << "    <materialref ref=\"Titanium\"/>"    << endl;
  of << "    <solidref ref=\"TPCWireCommon\"/> " << endl;
  of << "  </volume> "                           << endl;

  //now the plane itself
  of << "  <volume name=\"volTPCPlane\"> "    << endl;
  of << "    <materialref ref=\"LAr\"/>"      << endl;
  of << "    <solidref ref=\"TPCPlane\"/> "   << endl;

  //the wires at the -z end
  for(int i = 0; i < nsw; ++i){
    of << "    <physvol>"                                                     << endl;
    of << "     <volumeref ref=" << TString::Format("\"volTPCWire%d\"/> ", i) << endl;
    of << "     <position name=" << TString::Format("\"posTPCWire%d\" unit=\"cm\" y=\"-0.5*kTPCWidth+0.5*(%d+1)*kTPCWireXPitch*kCos60/kCos30\" z=\"-0.5*kTPCLength+0.5*kTPCWireXPitch*(%d+1)\" x=\"0\"/> ", i, i, i)                                                                                          << endl;
    of << "     <rotationref ref=\"rPlus150AboutX\"/>"                        << endl;
    of << "    </physvol> "                                                   << endl;
  }

  for(int i = 0; i < ncw; ++i){
    of << "    <physvol>"                                                    << endl;
    of << "     <volumeref ref=\"volTPCWireCommon\"/>"                           << endl;
    of << "     <position name=" << TString::Format("\"posTPCWire%d\" unit=\"cm\" y=\"-0.5*kTPCWidth+0.5*%d*kTPCWireXPitch*kCos60/kCos30\" z=\"-0.5*kTPCWirePlaneLength+kTPCWireXPitch*(0.5*%d + %d+1)\" x=\"0\"/> ", nsw+i, nsw, nsw, i)                                                                    << endl;
    of << "     <rotationref ref=\"rPlus150AboutX\"/>"                       << endl;
    of << "    </physvol> "                                                  << endl;
  }

  //the wires at the +z end
  for(int i = 0; i < nsw; ++i){
    of << "    <physvol>"                                                         << endl;
    of << "     <volumeref ref=" << TString::Format("\"volTPCWire%d\"/> ", nsw-i-1) << endl;
    of << "     <position name=" << TString::Format("\"posTPCWire%d\" unit=\"cm\" y=\"0.5*kTPCWidth-0.5*(%d+1)*kTPCWireXPitch*kCos60/kCos30\" z=\"0.5*kTPCWirePlaneLength-0.5*kTPCWireXPitch*(%d+1)\" x=\"0\"/> ", ncw+nsw+i, nsw-i-1, nsw-i-1)                                                                       << endl;
    of << "     <rotationref ref=\"rPlus150AboutX\"/>"                            << endl;
    of << "    </physvol> "                                                       << endl;
  }

  of << "  </volume> "                            << endl;
  of << "</structure> "                           << endl;

  of.close();
}
