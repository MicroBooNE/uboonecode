//---------------------------------------------------------------------------
// File:        gmtuple.cpp
// Description: Create a test root-tuple for testing GM browser - based
//              on hsimple.cxx from the Root distribution
// Created:     18-Oct-2002 Harrison B. Prosper
//---------------------------------------------------------------------------
#include <iostream>
#include <sstream>
#include <time.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TRandom.h>
#include <TSystem.h>

using namespace std;

int main(int argc, char **argv)
{
  cout << endl << "Welcome to gmtuple" << endl;

  int total = 3000;
      
  TFile   *file   = new TFile("tuple.root","RECREATE");
  TNtuple *ntuple = new TNtuple("ntuple","Test Tuple","x:y:z");

  cout << "Starting" << endl << endl;

  float x, y, z;
  for (int i = 0; i < total; i++) 
    {
      gRandom->Rannor(x, y);
      z = x*x + y*y;
      ntuple->Fill(x, y, z);
      if ( i % 10 == 0 )
	{
	  file->Write("ntuple", TObject::kOverwrite);
	  if (gSystem->ProcessEvents()) break;	  	    
	  gSystem->Sleep(1000);

	  time_t t = time(0);
	  cout << string(ctime(&t)).substr(0,24) 
	       << "\tgmtuple " << i << endl;
	}
     }
  file->Close();
}

