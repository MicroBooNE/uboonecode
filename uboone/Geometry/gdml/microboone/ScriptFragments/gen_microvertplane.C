#include "Riostream.h"

void gen_microvertplane(){

  ofstream of("micro-vertplane.gdml");

  int nwi = 1;

  //define the solids
  of << "<solids>" << endl;

  //wires for the plane
  of << "<tube name=\"TPCWireVert\" "        << endl;
  of << "  rmax=\"0.5*kTPCWireThickness\" "  << endl;
  of << "  z=\"kTPCWidth\" "                 << endl;
  of << "  deltaphi=\"2*kPi\" "              << endl;
  of << "  aunit=\"rad\" "                   << endl;
  of << "  lunit=\"cm\"/> "                  << endl;
  
  //the tpc wire plane
  of << "<box name=\"TPCPlaneVert\" "        << endl;
  of << "  x=\"kTPCWirePlaneThickness\" "    << endl;
  of << "  y=\"kTPCWidth\" "                 << endl;
  of << "  z=\"kTPCLength\" "                << endl;
  of << "  lunit=\"cm\"/> "                  << endl;

  //end solids
  of << "</solids>" << endl;

  //now defined the structures
  of << "<structure>" << endl;
  
  //the wires 
  of << "  <volume name=\"volTPCWireVert\"> "  << endl;
  of << "    <materialref ref=\"Titanium\"/>"  << endl;
  of << "    <solidref ref=\"TPCWireVert\"/> " << endl;
  of << "  </volume> "                         << endl;

  //now the plane itself
  of << "  <volume name=\"volTPCPlaneVert\"> "  << endl;
  of << "    <materialref ref=\"LAr\"/>"        << endl;
  of << "    <solidref ref=\"TPCPlaneVert\"/> " << endl;

  //the wires 
  for(int i = 0; i < nwi; ++i){
    of << "    <physvol>"                              << endl;
    of << "     <volumeref ref=\"volTPCWireVert\"/> "  << endl;
    of << "     <position name=" << TString::Format("\"posTPCWireVert%d\" unit=\"cm\" z=\"-0.5*kTPCLength+kTPCWirePitch*(%d+1)\" x=\"0\" y=\"0\"/> ", i, i)    << endl;
    of << "     <rotationref ref=\"rPlus90AboutX\"/> " << endl;
    of << "    </physvol> "                            << endl;
  }


  of << "  </volume> "                            << endl;
  of << "</structure> "                           << endl;

  of.close();
}
