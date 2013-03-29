// Script to place PMT's at specified coordinates 
// in the cryostat volume.  

#include "Riostream.h"
#include "vector"
#include "map"
#include "sstream"
#include "cmath"

void gen_cryostat(){

  float TankLength=1400;
  float PMTOffset=500;
  
  // Open input and output files

  ofstream of("micro-cryostat.gdml");
  ifstream pmts("pmt-coordinates.csv");


  // Read PMT positions from pmt-coordinates file

  double x,y,z,phi,r;
  std::vector<double> xs,ys,zs,phis;
  std::vector<std::string> RotationNames;
  std::map<double, int> RotationIndices;
  char d, d1[256];
  for(;;1)
    {
      pmts >> x >> d >> y >> d >> z >> d >> phi >> d >> r;
      xs.push_back(x);
      ys.push_back(y);
      zs.push_back(z);
      phis.push_back(phi);
      if(pmts.eof())
	break;
      pmts.getline(d1,256);
    }
  
  // Define rotations
  
  RotationNames.push_back("null");   // blank row, to distinguish index 0 from index null
  of << "<define>" << std::endl;
  for(int i=0; i!=phis.size(); i++)
    {
      std::stringstream TheName("");
      phi = phis[i];
      if(!RotationIndices[phi]) 
	{
	  TheName<<"rPMTRotation"<<RotationNames.size();
	  of<<"  <rotation name=\""<<TheName.str() << "\"  unit=\"deg\" x=\"0\"  y=\"90\"   z=\""<< phi << "\"/>"<<std::endl;
	  RotationIndices[phi]=RotationNames.size();
	  RotationNames.push_back(std::string(TheName.str()));
	}
    }
  of<<"</define>"<<std::endl;

 
 // Write top part of cryostat file (non PMT part)

  of<<"<solids>"                               <<std::endl;
  of<< "<tube name=\"Cryostat\""               <<std::endl;
  of<<"  rmax=\"(240)\"  "                     <<std::endl;
  of<<"  z=\""<<TankLength<<"+5\"  "           <<std::endl;
  of<<"  deltaphi=\"2*(3.1415926535897)\" "    <<std::endl;
  of<<"  aunit=\"rad\" "                       <<std::endl;
  of<<"  lunit=\"cm\"/> "                      <<std::endl;
  of<<"<tube name=\"SteelTube\" "              <<std::endl;
  of<<"  rmin=\"(200)\" "                      <<std::endl;
  of<<"  rmax=\"(240)\" "                      <<std::endl;
  of<<"  z=\"1400\" "                          <<std::endl;
  of<<"  deltaphi=\"2*(3.1415926535897)\" "    <<std::endl;
  of<<"  aunit=\"rad\" "                       <<std::endl;
  of<<"  lunit=\"cm\"/> "                      <<std::endl;
  of<<"</solids> "                             <<std::endl;

  of<<"<structure> "                           <<std::endl;
  of<<" <volume name=\"volSteelTube\"> "       <<std::endl;
  of<<"  <materialref ref=\"STEEL_STAINLESS_Fe7Cr2Ni\"/> " <<std::endl;
  of<<"  <solidref ref=\"SteelTube\"/> "       <<std::endl;
  of<<" </volume> "                            <<std::endl;
  of<<" <volume name=\"volCryostat\"> "        <<std::endl;
  of<<"  <materialref ref=\"LAr\"/> "          <<std::endl;
  of<<"  <solidref ref=\"Cryostat\"/> "        <<std::endl;

  // Write PMT part

  for(int i=0; i!=xs.size(); i++)
    {
      std::cout<<"Placing PMT at [x,y,z,phi]" << xs[i] << " " <<ys[i] << " " <<zs[i] <<" " <<phis[i]<<std::endl;
      // Find coordinate of centre of PMT vol, given lens position
      phi = phis[i];
      x = xs[i] + (2*2.54) * cos(phi); 
      y = ys[i] + (2*2.54) * sin(phi); 
      z = zs[i] - PMTOffset;

      of<<"  <physvol>"                                            <<std::endl;
      of<<"   <volumeref ref=\"volPMT\"/>"                         <<std::endl;
      of<<"   <position name=\"posPMT"<<i<<"\" unit=\"cm\" x=\""<<x<<"\" y=\""<<y<<"\" z=\""<<z<<"\"/>"        <<std::endl;
      of<<"   <rotationref ref=\""<<(RotationNames.at(RotationIndices[phi]))<<"\"/>"     <<std::endl;
      of<<"  </physvol>"                                           <<std::endl;

    }
  

  // Write footer part
  
  of<<"  <physvol> "                           <<std::endl;
  of<<"   <volumeref ref=\"volSteelTube\"/> "  <<std::endl;
  of<<"   <position name=\"posSteelTube\" unit=\"cm\" x=\"0\" y=\"0\" z=\"0\"/> " <<std::endl;
  of<<"  </physvol> "                          <<std::endl;
  of<<"  <physvol> "                           <<std::endl;
  of<<"   <volumeref ref=\"volTPC\"/> "        <<std::endl;
  of<<"   <position name=\"posTPC\" x=\"0.0\" y=\"0\" z=\"0\"/> " <<std::endl;
  of<<"  </physvol> "                          <<std::endl;
  of<<" </volume> "                            <<std::endl;
  of<<"</structure> "                          <<std::endl;



}
  



