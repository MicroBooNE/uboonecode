//---------------------------------------------------------------------------
// File:      gmrowser.cpp
// Created:   17-Oct-2002 Harrison B. Prosper
//---------------------------------------------------------------------------
#include <TApplication.h>
#include <TQObject.h>
#include <iostream>
#include "gmbrowser/GMBrowser.hpp"

using namespace std;

int main(int argc, char **argv)
{
  if ( argc < 2 )
    {
      cout << "Usage:\n\tGMbrowser <config-file> [title] -width <width> "
           << " -height <height> \n";
      return 1;
    }

  const char *config = argv[1];
  const char *title;

  // Define title
  string str("Global Monitor Browser");
  if ( argc > 2 && argv[2][0] != '-')
    title  = argv[2];
  else
    title = str.c_str();

  // Get command line options
  int iarg;
  int width  = 0;
  int height = 0;

  for (iarg = 1; iarg < argc; ++iarg) {
    char chop0 = argv[iarg][0];
    char chopt = argv[iarg][1];

    if (chop0 == '-') {
      switch ( chopt ) {
        // Set width and height of screen from command line options.
        case 'w' : width = atoi(argv[++iarg]);  break;
        case 'h' : height = atoi(argv[++iarg]); break;
      }
    }
  }

  // Instantiate TApplication.
  int argtmp = 1;
  TApplication a("GM", &argtmp, argv);

  // ECC - get screen size and scale window to that size.
  Int_t x, y;  UInt_t w, h;
  gVirtualX->GetGeometry(-1, x, y, w, h);

//  if (width == 0)  width=(int)(w*0.95); 
  if (height == 0) {
    height=(int)(h*0.95);
    width =(int)(height*0.85 + 150);  // 150 is the width of the listbox.
  }
  

  // Default values if no values have been set.
  if (width == 0 || height == 0) {
    width=600;
    height=600;
  }

  cout << "Width, height: " << width << " " << height << endl;

  // Instantiate GMBrowser.
  cout << "Title: " << title << endl;
  GMBrowser g(title, (UInt_t)width, (UInt_t)height);
  cout << "\nWelcome to GM Browser " << g.Version() << endl;
  cout << "Configuration file: " << config << endl;
  
  // Configure gmbrowser
  g.Configure(config);

  // Start gmbrowser and TApplication.
  g.Run();
  a.Run();
  return 0;
}
