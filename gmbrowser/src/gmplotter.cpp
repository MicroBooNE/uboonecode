//---------------------------------------------------------------------------
// File:      gmrowser.cpp
// Created:   17-Oct-2002 Harrison B. Prosper
//---------------------------------------------------------------------------
#include <TApplication.h>
#include <TQObject.h>
#include <iostream>
#include "gmbrowser/GMPlot.hpp"

using namespace std;

int main(int argc, char **argv)
{
  if ( argc < 2 )
    {
      cout << "Usage:\n\tgmplotter <config-file>\n";
      return 1;
    }

  // Determine whether to make one pass or loop infinitely
  bool LOOP=kTRUE;
  if ( argc == 3) LOOP = kFALSE;
 
  const char *config = argv[1];

  // Instantiate GMBrowser.
  GMPlot *g = new GMPlot(kTRUE);
  cout << "\nWelcome to GM Plotter " << g->Version() << endl;
  cout << "Configuration file: " << config << endl;
  
  g->Configure(config);
  g->Run(LOOP);
  return 0;
}
