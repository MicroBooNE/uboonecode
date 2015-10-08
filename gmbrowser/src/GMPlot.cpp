/////////////////////////////////////////////////////////////////////////////
// File:         GMplotter.cpp
// Description:  Routines to generate .gif files based on .cfg file.
//               Based upon gmbrowser
// Created:      January 12, 2005 Elliott Cheu
// Modified:     February 3, 2006 E. Cheu
//                 Updated histogram options.
//               February 25, 2007 E. Cheu
//                 Improved file handling.
//                 Write out plots to .ps area also.
//               24 Aug 2010 GAF Separate run and subrun number.
//               07 Oct 2010 GAF Support for TH1D, TH2D, THProfile and THProfile2D.
//                               Do not show run and subrun number in histogram title.
//////////////////////////////////////////////////////////////////////////////
//$Revision: 1.8 $

#include "gmbrowser/GMPlot.hpp"

// Since there is no possibility of conflict, we expose the full std namespace.

using namespace std;


// FUNCTIONS
////////////

extern void Info(const char *location, const char *msgfmt, ...);
extern void Warning(const char *location, const char *msgfmt, ...);
extern void Erorr(const char *location, const char *msgfmt, ...);

string getTime()
{
  time_t t = time(0);
  return string(ctime(&t)).substr(0,24);
}


//---------------------------------------------------------------
// Remove leading and trailing space. (Couldn't get Root's 
// version to work)
//---------------------------------------------------------------
string strip(string line)
{
  int l = line.size();
  if ( l == 0 ) return string("");

  int n = 0;
  while (((line[n] == 0)    || 
          (line[n] == ' ' ) || 
          (line[n] == '\n') || 
          (line[n] == '\t')) && n < l) n++;

  int m = l-1;
  while (((line[m] == 0)    || 
          (line[m] == ' ')  || 
          (line[m] == '\n') || 
          (line[m] == '\t')) && m > 0) m--;
  return line.substr(n,m-n+1);
}
string strip(const char *line)
{
  return strip(string(line));
}

//---------------------------------------------------------------
// Split a string into two, at specified single character delimeter
// Default is to split at first tab or space.
//---------------------------------------------------------------
void split(string line, string &a, string &b, char symbol=0)
{
  int l = line.size();
  if ( l == 0 ) return;

  int n = 0;
  while ( n < l )
    {
      if (symbol == 0)
	{
	  if ( (line[n] == ' ') || (line[n] == '\t') )
	    break;
	  else
	    n++;
	}
      else 
	{
	  if (line[n] == symbol)
	    break;
	  else
	    n++;
	}
    }

  a = line.substr(0,n);
  if ( n < l )
    b = line.substr(n+1,l-n-1);
  else
    b = "";
}

//---------------------------------------------------------------
// send command to operating system and get reply
//---------------------------------------------------------------
string gmpopen(const char *cmd)
{  
  FILE *f = popen(cmd,"r");
  char s[4096];
  int n = fread(s,1,4096,f);
  pclose(f);
  return strip(string(s).substr(0,n));
}


//---------------------------------------------------------------
// Extract run number from histogram contents.
//---------------------------------------------------------------
int GMPlot::getRunNumber(const char *file, TFile *rootfile)
{
  string filename(file);

  // ECC - 4/1/03 - For PhysEx we can use the histogram, h_event_runno
  // ECC - 5/16/03 - now both TrigEx and PhysEx produce this histogram.
  if (fDebug) Info("getRunNumber","Read keys");
  rootfile->ReadKeys();
  if (fDebug) Info("getRunNumber","get histogram.");
  TH1D *physhist = (TH1D*)rootfile->Get(fRunNumberHist.c_str());
  if (fDebug) Info("getRunNumber","physhist = %d", physhist);
  if (physhist) {
    if (fDebug) Info("getRunNumber","Get bin contents");
    int irun = (int)physhist->GetBinContent(1);
    if (irun) return irun;
  }
  else return 0;
}

//---------------------------------------------------------------
// Get number of events from a specific histogram.
//---------------------------------------------------------------
int GMPlot::getTotalEvents(TFile *rootfile)
{

  // Read keys from root file.
  if (fDebug) Info("getTotalEvents","Read keys");
  rootfile->ReadKeys();

  // Determine the histogram.
  if (fDebug) Info("getTotalEvents","get histogram.");
  if (fTotalEventHist.c_str() != "") {
    TH1D *hist = (TH1D*)rootfile->Get(fTotalEventHist.c_str());

    // If histogram exists, then return integral of histogram..
    if (fDebug) Info("getRunNumber","hist = %d", hist);
    if (hist)  return (int)hist->GetEntries();    
  }

  // Return zero if histogram doesn't exist.
  return 0;
}


// CONSTANTS
////////////

string VERSION = "v1.00.00, 12 January 2005";

// E. Cheu - Used for determining histogram file types
enum HistFiles {
  K_DATA_FILE,
  K_REF1_FILE,
  K_REF2_FILE
};

int K_HISTS[K_MAX_FILE];              // histogram index number.

string configFile_curr;               // Current configuration file.

const float K_TITLE_SIZE     = 0.06;  // Size of title
const float K_TITLE_OFFSET   = 0.75;  // Offset of title from plot
const float K_PADYFRAC       = 0.95;  // Fractional height of pad
const float K_STAT_WIDTH     = 0.30;  // Width of stat button
const float K_STAT_HEIGHT    = 0.20;  // Height of stat button
const int K_LINE_WIDTH       = 2;     // in pixels
const int K_REF_COLOR        = 2;     // Color of reference plots
const int K_REF2_COLOR       = 4;     // Color of reference plots
const int K_REF2_STYLE       = 2;     // Dashed.
const int   K_SCREEN_WIDTH   = 1500;//900;   // width of canvas
const int   K_SCREEN_HEIGHT  = 1400;//925;   // height of canvas


//---------------------------------------------------------------
// Constructor
//---------------------------------------------------------------
GMPlot::GMPlot(Bool_t debug)
{
  // Initialize the error level.
  gErrorIgnoreLevel = 1001;

  // Make sure that no window is created.
  gROOT->SetBatch(kTRUE);


  // Initialize some internal variables

  fPageList.clear();
  fWWWDir = StrDup(".");
  fGMPlotPSDir = StrDup(".");
  fPad = 0;
  fCurrentPage = 0;
  fRunNumber   = 0;
  fRun         = 0;
  fSubRun      = 0;
  fTotalEventHist = "";
  fRunNumberHist  = "";

  // ECC - 25 Apr 2003 - need to loop over all available files.
  for (int i=0; i<K_MAX_FILE; i++) {
    fRootFile[i] = 0;
  }
  fRootFilename.clear();

  gStyle->SetPalette(1);
  // Create a Canvas
  fECanvas = new TCanvas("GM", "GM", K_SCREEN_WIDTH, K_SCREEN_HEIGHT);  

}



// METHODS
//////////

const char *GMPlot::Version(){return VERSION.c_str();}



//---------------------------------------------------------------
// Return a pointer to the TCanvas 
//---------------------------------------------------------------
TCanvas  *GMPlot::Canvas() 
{
  if ( fECanvas )
    return fECanvas->GetCanvas();
  else
    return (TCanvas *)0;
}

//---------------------------------------------------------------
// Return a pointer to specified page
//---------------------------------------------------------------
GMPage   *GMPlot::Page(int pageNumber) 
{
  if ( pageNumber < 0 ) pageNumber = fCurrentPage;

  GMPageList::iterator thePage = fPageList.begin();
  int count = 0;
  while ( thePage != fPageList.end() )
    {
      if ( count == pageNumber ) return &(*thePage);
      thePage++;
      count++;
    }
  return (GMPage *)0;
}

//---------------------------------------------------------------
// Return specified histogram
//---------------------------------------------------------------
TObject  *GMPlot::Histogram(const char *name, int which)
{
  TH1 *temp;
  // ECC - 13 May 2003 - Be sure to read keys before getting name.
  TFile *file = fRootFile[which];
  if ( file ) {
    int ntimes = 0;
    int nkeys  = 0;
    
    while (ntimes < 3 && nkeys == 0) {
      //cout << "ReadKeys" << endl;
      nkeys = file->ReadKeys();
      if (nkeys ==0) gSystem->Sleep(50);
    }

    // ECC - 1 Aug 2003 - Make further checks to try to protect against
    //                    bad reads.
   // cout << "Get(name)" << endl;

    ntimes = 0;
    while (ntimes < 3) {
    //  cout << "Get(name)" << endl;
      temp = (TH1*) file->Get(name);

      if (temp == 0 || (temp->TestBit(TObject::kNotDeleted) == 0) ) {
        ntimes++;
        gSystem->Sleep(50);
        nkeys = file->ReadKeys();
      }
      // Otherwise this should be a good pointer.
      else return temp;
    }

    // Failed to get good TObject.
    return (TObject *)0;

  }
  else return (TObject *)0;
} 

// ECC - 25 Feb 2007 - new version to clean up file handling.
Bool_t GMPlot::openFiles(int sumfile)
{
  // ECC - 25 Apr 2003 - Modify to allow multiple .root files.

  // Loop over all possible histogram files (only 3)
  for (int K=0; K<K_MAX_FILE; K++) {

    // Close files first. This is necessary because we might be reading files 
    // over the network and the file might be a soft link. In this case
    // if the softlink has been moved, we might actually be reading the
    // wrong file. So, we have to reopen the file each time we display
    // a new page.
    if ( fRootFile[K] ) {
      if (fDebug) {
        cout << "Closing File: " 
             << fRootFilename[K_HISTS[K]].c_str() << " " << fRootFile[K] << endl;
      }
      delete fRootFile[K];
      fRootFile[K] = 0;
    }

    // File to open is determined by K_HISTS array.
    int KREAD = K_HISTS[K];

    // Now try to open new file
    if (KREAD >= 0 && fRootFilename[KREAD].c_str() != "") {
      fRootFile[K] = TFile::Open(fRootFilename[KREAD].c_str(), "readonly"); 
      if (!fRootFile[K]) {
        Error("reopen","Unable to open file <%d>", 
              fRootFilename[KREAD].c_str());
        return kFALSE;
      }
      else {
	// ECC - 1 May 2003 -Check keys in file.
        if (!checkKeys(K)) return kFALSE;
        if (fDebug) 
          cout << "Opened: " << fRootFilename[KREAD].c_str() << endl;
      }
    }
  }


  // E. Cheu - close and open the files to be used for summing.
  int nhists = fsumFiles[sumfile-1].size();

  for (int ihist = 0; ihist < nhists; ihist++) {

    // Determine file to add.
    TFile *file = fsumFiles[sumfile-1][ihist];
    if (!file) continue;

    // Close the file. This allows gmbrowser to continue to run
    // even if the soft link has changed.
    if (file->IsOpen()) delete fsumFiles[sumfile-1][ihist];

    // Now open file again.
    fsumFiles[sumfile-1][ihist] = TFile::Open(fsumName[sumfile-1][ihist].c_str(),
      "readonly");
  }

  return kTRUE;
}


GMPageList *GMPlot::PageList(){return &fPageList;}


//---------------------------------------------------------------
// Read and parse configuration file and create pages specified
// therein.
//---------------------------------------------------------------
void GMPlot::Configure(const char *cFile)
{

  if ( fDebug )
    Info("Configure","START");

  // Set plot name prefix from name of configuation file.
  plotTitle(cFile);

  char *configFile = gSystem->ExpandPathName(cFile);

  // Read configuration file

  ifstream *input = new ifstream(configFile);
  if ( !input || (input && (input->fail() || input->bad()) ) )
    {
      Error("configure","Unable to read configuration file<%s>",configFile);
      return;
    }
  if ( fDebug )
    Info("configure", "Opened configuration file<%s>", configFile);

  write("Opened configuration file %s", configFile);

  configFile_curr = configFile;
  
  Bool_t FoundPageDivision = kFALSE;
  char line[fLineSize];
  int  page = 0;
  GMPage *thePage = 0;
  total_sumFiles = 0;


  // Loop over file and parse entries

  while ( input->getline(line, fLineSize) )
    {
      string str(strip(line));         // Strip away leading and trailing space

      string keyword, value;
      split(str, keyword, value, '#'); // Strip away comments
      str = keyword;

      if ( fDebug)
       	Info("configure","keyword<%s>", keyword.c_str());

      if ( str.size() == 0 ) continue; // Ignore blank lines

      split(str, keyword, value);
      value = strip(value);

      // ECC - 25 Apr 2003 - add option to turn on/off debugging.
      if ( keyword == "Debug.On:" )
	{
          int idebug = atoi(value.c_str());
          if (idebug) {
            fDebug = kTRUE;
            gErrorIgnoreLevel = 0;
	  }
          else
            fDebug = kFALSE;
          Info("configure","Setting Debug: %d ", fDebug);
	}
      else if      ( keyword == "Root.File:" )
	{
	  //	  if ( fDebug)
	  //	    Info("configure","Root.File<%s>", value.c_str());
	  //          cout << "Root File: " << value << endl;

	  load(value, K_DATA_FILE, &K_HISTS[K_DATA_FILE]);
	  // E. Cheu - need to initialize the addfile array.
          addfile(value, 0);
	  if ( fRootFile[K_DATA_FILE] != 0 ) {
  	    // ECC - new version.
	    //          fRunNumber = getRunNumber(value.c_str());
            fRunNumber = getRunNumber(value.c_str(), fRootFile[K_DATA_FILE]);
            // ECC - initialize reference files.
            K_HISTS[K_REF1_FILE] = -1;
            K_HISTS[K_REF2_FILE] = -1;
            if (fDebug) Info("configure","Run Number<%d>", fRunNumber);
	  }
	}

      // E. Cheu - 24 Mar 2006 - Adds the histograms from this line
      //                         to the ones from Root.File: <name>
      else if ( keyword == "Root.File.Add:" )
	{
	  addfile(value, 1);
	}

      else if ( keyword == "Root.File.Ref:" )
	{
	  //          cout << "Ref  File: " << value << endl;
	  //	  if ( fDebug)
	  //	    Info("configure","Root.File.Ref<%s>", value.c_str());
	  load(value, K_REF1_FILE, &K_HISTS[K_REF1_FILE]);
	}
      // E. Cheu - 25 Feb 2007 - Add code for second reference histogram.
      else if ( keyword == "Root.File.Ref2:" )
	{
	  load(value, K_REF2_FILE, &K_HISTS[K_REF2_FILE]);
	}


      else if ( keyword == "WWW.Dir:" )
	{
	  if ( fDebug)
	    Info("configure","WWW.Dir<%s>", value.c_str());

	  fWWWDir = StrDup(gSystem->ExpandPathName(value.c_str()));
	}

      //  Directory to store .ps files.
      else if ( keyword == "GMPlot.PS.Dir:" ) {
        if ( fDebug) Info("configure","GMPlot.PS.Dir<%s>", value.c_str());
        fGMPlotPSDir = StrDup(gSystem->ExpandPathName(value.c_str()));
      }

      else if ( keyword == "Load.Dir:" )
	{
	  if ( fDebug)
	    Info("configure","Load.Dir<%s>", value.c_str());

	  // Store load directory.
	  fFileLoadInfo.fIniDir = StrDup(gSystem->
                                         ExpandPathName(value.c_str()));
	}

      else if ( keyword == "Save.Dir:" )
	{
	  if ( fDebug)
	    Info("configure","Save.Dir<%s>", value.c_str());

	  fFileSaveInfo.fIniDir = StrDup(gSystem->
					 ExpandPathName(value.c_str()));
	}

      else if ( keyword == "AutoSave.On:" )
	{
	  //	  fAutoSave = atoi(value.c_str());
	  if ( fDebug)
	    Info("configure","AutoSave.On<%d>", fAutoSave);
	}
      else if ( keyword == "Update.On:" )
	{
	  fUpdate = atoi(value.c_str());
	  if ( fDebug)
	    Info("configure","Update.On<%d>", fUpdate);
	}
      else if ( keyword == "Cycle.On:" )
	{
	  fCycle = atoi(value.c_str());
	  if ( fDebug)
	    Info("configure","Cycle.On<%d>", fCycle);
	}

      else if ( keyword == "Cycle.Plots:" )
        {
          fCyclePlots = atoi(value.c_str());
          if ( fDebug)
            Info("configure","Cycle.Plots<%d>", fCycle);
	}

      else if ( keyword == "AutoSave.Period:" )
	{
	  //	  fAutoSavePeriod = atoi(value.c_str());
	  if ( fDebug)
	    Info("configure","AutoSave.Period<%d>", fAutoSavePeriod);

	  //          fAutoSavePeriod = 1000 * fAutoSavePeriod;
	}
      else if ( keyword == "Update.Period:" )
	{
	  fUpdatePeriod = atoi(value.c_str());
	  if ( fDebug)
	    Info("configure","Update.Period<%d>", fUpdatePeriod);

          fUpdatePeriod = 1000 * fUpdatePeriod;
	}
      else if ( keyword == "Cycle.Period:" )
	{
	  fCyclePeriod  = atoi(value.c_str());
	  if ( fDebug)
	    Info("configure","Cycle.Period<%d>", fCyclePeriod);

          fCyclePeriod = 1000 * fCyclePeriod;
	}
      else if ( keyword == "TotalEvent.Hist:" ) {
        fTotalEventHist = value;
      }
      else if ( keyword == "RunNumber.Hist:" ) {
        fRunNumberHist = value;
      }
      else if ( keyword == "Page.Title:" )
	{
	  FoundPageDivision = false;	  

	  page = AddPage(value.c_str());
	  //	  fListBox->AddEntry(value.c_str(), page);
          
	  if ( fDebug)
	    Info("configure","Page.Title<%s><%d>", value.c_str(), page);

	  thePage = Page(page); // Get pointer to current page
	  if ( ! thePage )
	    {
	      Warning("configure","Unable to get page<%d>", page);
	      input->close();
	      return;
	    }
          thePage->runMacro = 0;
          thePage->ShowRunNumber = 0;

          // ECC - 19 Sep 2006 - Initialize number of summed files.
          thePage->SumFileNum = 1;
	}
      // ECC - 6 Aug 2003 - Add code to allow people to run macros
      //                    from within a page.
      // This is the name of the file containing the macro.
      else if ( keyword == "Page.Macro.File:") {
        thePage->MacroFile = value;
        thePage->runMacro = 1;
        thePage->FileNumber=K_HISTS[K_DATA_FILE];
        thePage->RefFileNum=K_HISTS[K_REF1_FILE];
        thePage->RefFileNum2=K_HISTS[K_REF2_FILE];
	// Load macro at this time.
        string exp_value = StrDup(gSystem->ExpandPathName(value.c_str()));
        string tmpstr = ".L " + exp_value;
        gROOT->ProcessLine(tmpstr.c_str());
      }
      // ECC - 21 June 2004 - add in ability to associate help file with page.
      else if ( keyword == "Page.HelpFile:" ) {
        thePage->HelpFile = value;
      }
      // ECC - 6 Aug 2003 - This is the routine within MacroFile.
      //                    Just specify the name (don't include "()").
      else if ( keyword == "Page.Macro.Func:") {
        thePage->MacroFunc = value;
      }
      // Turn on option to display run number on page.
      else if ( keyword == "Page.ShowRunNumber:") {
        thePage->ShowRunNumber = atoi(value.c_str());
      }
      else if ( keyword == "Page.Division:" )
	{
	  FoundPageDivision = true;

	  if ( ! thePage )
	    {
	      Warning("configure","Unable to get page<%d>", page);
	      input->close();
	      return;
	    }
	  istringstream in(value.c_str()); int ix, iy; in >> ix >> iy;
	  thePage->DivX = ix;
	  thePage->DivY = iy;

	  if ( fDebug)
	    Info("configure", "Page.Division<%d, %d><%d>", ix, iy, page);
	}
      else if ( FoundPageDivision )
	{
	  if ( ! thePage )
	    {
	      Warning("configure","Unable to get page<%d>", page);
	      input->close();
	      return;
	    }

          // Get histogram name, title and Draw option 
	  // ECC - 6 Aug 2003 - Change delimiter to a "|"

	  string name, tmp, title, option;
	  split(str, name, tmp, '|');
	  name  = strip(name);
	  tmp   = strip(tmp);
          split(tmp, title, option, '|');
          title = strip(title);
          option= strip(option);

          // ECC - write out a message if we can't find the histogram.
          TH1* temp = (TH1*) fRootFile[K_DATA_FILE]->Get(name.c_str());
          if (temp == 0 || (temp->TestBit(TObject::kNotDeleted) == 0)) {
            Error("configure", 
                  "Histogram<%s> not found in RootFile<%s>.\n"
	          "Check your configuration file<%s>", 
	          name.c_str(), fRootFile[K_DATA_FILE]->GetName(), configFile);
	  }


	  thePage->Name.push_back(name);
	  thePage->Title.push_back(title);
          thePage->Option.push_back(option);
	  // ECC - 25 Apr 2003 - save number of histogram.
          thePage->FileNumber=K_HISTS[K_DATA_FILE];
          thePage->RefFileNum=K_HISTS[K_REF1_FILE];
          thePage->RefFileNum2=K_HISTS[K_REF2_FILE];
          thePage->SumFileNum=total_sumFiles;

	  if ( fDebug)
	    Info("configure", "\t<%s>Hist<%s>", thePage->PageName.c_str(),
		 name.c_str());		 
	}
    }
  input->close();

  if ( fDebug )
    Info("configure", "Closed configuration file<%s>", configFile);

  if (fDebug) {
      // First close all open files.
    for (int K=0; K<K_MAX_FILE; K++) {
      cout << "K: " << K << " fRootFile[K]: " << fRootFile[K] << endl;
    }
  }

  // ECC - The cycle mode has precedence over Update. So, turn
  //       off updating the current page if the cycle mode is on.
  if (fCycle && fUpdate) fUpdate = kFALSE;

  // Draw first page

  fCurrentPage = 0;
  //  fListBox->MapSubwindows();
  //  fListBox->Layout();
  //  fListBox->Select(fCurrentPage);
  //  gClient->NeedRedraw(fListBox); // Force immediate redraw


  if ( fDebug)
    Info("Configure","END");
}

//---------------------------------------------------------------
// Add new pages to any existing ones and return page number
//---------------------------------------------------------------
int GMPlot::AddPage(const char *name, int divx, int divy)
{
  if ( fDebug)
    Info("AddPage","START");

  if ( fDebug )
    Info("AddPage","\tPage.Title<%s>", name);

  int page = fPageList.size();

  if ( fDebug )
    Info("AddPage","\tPage.Number<%d>", page);

  fPageList.push_back(GMPage());
  GMPage *thePage  = &(fPageList.back());
  if ( ! thePage )
    {
      Warning("AddPage","Unable to get new page");
      return -1;
    }

  thePage->PageName= name; 
  thePage->DivX    = divx;
  thePage->DivY    = divy;
  thePage->Name.clear();
  thePage->Title.clear();
  thePage->Option.clear();

  if ( fDebug)
    Info("AddPage","END");

  return page;
}


//---------------------------------------------------------------
// Open root-file.
// Read histograms for this page one by one, from the root-file, 
// and add each to a different pad until all pads are used up or
// until we run out of histograms.
// Close root-file
//---------------------------------------------------------------
void GMPlot::DrawPage(int *TotalEvents)
{
  int page;
  TFile *fDataFile, *fRefFile;

  if ( fDebug)
    Info("DrawPage","START");

  if ( fPageList.size() == 0 ) return;

  // Default is to draw current page

  //  if ( page < 0 )
  page = fCurrentPage;

  // Get specified page

  GMPage *thePage = Page(page);
  if ( ! thePage )
    {
      Warning("DrawPage","Unable to get page number %d", page);
      return;
    }
  
  if ( fDebug) {
    Info("DrawPage","\tPage<%s>Page<%d>,Div<%d,%d>", 
	 thePage->PageName.c_str(), 
	 page, thePage->DivX, thePage->DivY);
    Info("DrawPage","\tPage<%s>Page<%d>,Files<%d,%d>", 
	 thePage->PageName.c_str(), 
	 page, thePage->FileNumber, thePage->RefFileNum);
  }

  // Get canvas
  TCanvas *canvas = Canvas();
  if ( ! canvas )
    {
      Warning("AddPage","Canvas pointer is zero!!!");
      return;
    }

  // ECC - Set up the correct files.
  K_HISTS[K_DATA_FILE] = thePage->FileNumber;
  K_HISTS[K_REF1_FILE] = thePage->RefFileNum;
  K_HISTS[K_REF2_FILE] = thePage->RefFileNum2;
  
  ///////////////////////
  // Open histogram files
  ///////////////////////

  if ( ! openFiles(thePage->SumFileNum) ) return;

  // ECC - 6 Aug 2003 - Handle macros here.
  if (thePage->runMacro) {

    // Check if file exists and read keys.
    if ( fRootFile[K_DATA_FILE] ) {
      fRootFile[K_DATA_FILE]->ReadKeys();
      if (fRootFile[K_REF1_FILE]) fRootFile[K_REF1_FILE]->ReadKeys();
      if (fRootFile[K_REF2_FILE]) fRootFile[K_REF2_FILE]->ReadKeys();

      // Clear canvas and set background to white.
      canvas->Clear();
      canvas->SetFillColor(0);

      // Run the macro.
      //   Need to pass the pointers to the canvas, data and ref files to
      //     the macro.
      //   Use: Form("function((TFile*)%p)", tfile)
      string mstr = thePage->MacroFunc + 
                    "((TCanvas*)%p, (TFile*)%p, (TFile*)%p)";
      gROOT->ProcessLine( Form(mstr.c_str(), canvas, 
                          fRootFile[K_DATA_FILE], fRootFile[K_REF1_FILE]) );
    }

    // Update canvas to show picture.
    canvas->Update();
    return;
  }

  // Make sure number of pads is at least 1
  int npads = thePage->DivX * thePage->DivY;
  if ( npads < 1 )
    {
      Error("DrawPage","Zero pads in page %d! Check config file", page);
      return;
    }


  // Update run number just in case file has changed
  fRunNumber = getRunNumber(fRootFilename[K_HISTS[K_DATA_FILE]].c_str(), 
                            fRootFile[K_DATA_FILE]);
  
  // GAF - 24 Aug 2010 - separate run and subrun number 
  //                     fRunNumber=1000*fRun+fSubRun
  fRun    = (int)floor(fRunNumber/10000.);
  fSubRun = (int)(fRunNumber-10000*fRun);

  // Get total number of events processed.
  *TotalEvents = getTotalEvents(fRootFile[K_DATA_FILE]);


  PageStyle(canvas, thePage);
 
  int pad   = 1; 
  string draw_opt;
  string config_file;

  static int run_exec = 0;
  
  if ( fDebug) Info("DrawPage","Loop over histograms");
  for (unsigned int hist = 0; hist < thePage->Name.size(); hist++)
    {
      canvas->cd(pad);
      // ECC - 20 Aug 2003 - change to default statistics display.
      //                     This has to be done after changing to the next
      //                     pad.
      if (option_stat_value > 0) {
        gStyle->SetOptStat(1111);
      }

      // Zero various option variables.
      config_file = "";

      // Get histogram
      string name = thePage->Name[hist];
      if ( fDebug) Info("DrawPage","Get histogram %s", name.c_str());
      TH1 *h0 = (TH1 *)Histogram(name.c_str(), K_DATA_FILE);
      if ( ! h0 ) {
        Warning("DrawPage", "Unable to read histogram<%s>", name.c_str());
	pad++;
	if ( pad > npads ) break;
	continue;
      }

      TH1 *hst = SumHists(thePage->SumFileNum, h0, name.c_str());

      if ( fDebug)
	Info("DrawPage","\tDraw: %s", name.c_str());

      // Set histogram title
      string atitle = thePage->Title[hist];
      if ( atitle != "" ) hst->SetTitle(atitle.c_str());

      char title[132];
      // ECC - 2 Oct 2003 - only add run number to hist title if
      //                    run number is valid.
      //if (fRunNumber > 0 && !thePage->ShowRunNumber)
        // GAF - 24 Aug 2010 - separate run and subrun number 
        //sprintf(title, "Run: %6.6d - %s", fRunNumber, hst->GetTitle());
        //sprintf(title, "Run: %4.4d SubRun: %4.4d - %s", fRun, fSubRun, hst->GetTitle());
      //else
        sprintf(title, "%s", hst->GetTitle());

      hst->SetTitle(&title[0]);    
      hst->SetTitleSize(K_TITLE_SIZE); 
      hst->SetTitleOffset(K_TITLE_OFFSET);

      // Set Draw option
      hst->SetLineWidth(K_LINE_WIDTH);

      // Parse the options. 
      string option = thePage->Option[hist];
      histOptions(option, config_file, hst);

      draw_opt = "E0";
      if (option_selected) draw_opt = "";

      // ECC - 11 Nov 2003 - change line width.
      if (option_lwid > 0 && option_lwid < 10)  hst->SetLineWidth(option_lwid);
      // ECC - 20 Aug 2003 - change statistics box.
      if (option_stat_value > 0) {
        // update canvas so that previous pads have the correct
        // style.
        if (pad > 1) {
          canvas->Modified();
          canvas->Update();
	}
        gStyle->SetOptStat(option_stat_value);
      }
     
      // ECC - 25 Apr 2003 - protect against no reference file.
      // If reference histogram exists then overlay it on histogram
      TH1 *ref1;
      if (K_HISTS[K_REF1_FILE] > 0 && option_PlotRef) 
        ref1 = (TH1 *)Histogram(name.c_str(), K_REF1_FILE);      
      else
        ref1 = 0;

      // Do the same for the second reference histogram.
      TH1 *ref2;
      if (K_HISTS[K_REF2_FILE] > 0 && option_PlotRef2) 
        ref2 = (TH1 *)Histogram(name.c_str(), K_REF2_FILE);      
      else
        ref2 = 0;

      // ECC - 03 Oct 2005 - Add option to set X axis range.
      if (option_range_max > option_range_min) {
        hst->SetAxisRange(option_range_min, option_range_max);
      }

      // ECC - 22 Dec 2005 - Add option to rebin histogram.
      if (option_rebin > 0) {
        hst->Rebin(option_rebin);
      }
      if ( ref1 ) {
	// ECC - 20 Aug 2003 - make line width smaller for reference.
        ref1->SetLineWidth(K_LINE_WIDTH/2);
        ref1->SetLineColor(K_REF_COLOR);

        // ECC - 03 Oct 2005 - Add option to set X axis range.
        if (option_range_max > option_range_min) {
          ref1->SetAxisRange(option_range_min, option_range_max);
        }
        // ECC - 22 Dec 2005 - Add option to rebin histogram.
        if (option_rebin > 0) {
          ref1->Rebin(option_rebin);
        }

  	// Do the same for the 2nd reference histogram
        if (ref2) {
          ref2->SetLineWidth(K_LINE_WIDTH/2);
          ref2->SetLineStyle(K_REF2_STYLE);
          ref2->SetLineColor(K_REF2_COLOR);
          if (option_range_max > option_range_min) 
            ref2->SetAxisRange(option_range_min, option_range_max);
          if (option_rebin > 0) ref2->Rebin(option_rebin);
	}

        // Compute goodness of fit for 1 and 2-D histograms only
	if ( option_goodness && ((string(hst->ClassName()) == "TH1D") ||
	                    (string(hst->ClassName()) == "TProfile") || 
	                    (string(hst->ClassName()) == "TProfile2D") || 
	                  (string(hst->ClassName()) == "TH2D")) )
	  {
	    double prob  = hst->KolmogorovTest(ref1);
	    char title[132];
	    sprintf(title, "%s (Prob: %1.3f) ", hst->GetTitle(), prob);
	    hst->SetTitle(&title[0]);
	  }

        // Superimpose reference histograms only for 1-D histograms
	if ( string(hst->ClassName()) == "TH1D" || string(hst->ClassName()) == "TProfile") {
          // Only rescale reference histogram if scaleRef is nonzero.
          if (option_scaleRef) {
	    double norm  = hst->Integral();
	    double rnorm = ref1->Integral();
	    double r2norm = 0.0; 
  	    double scale = 1.0;
            // ECC - 1 May 2003 - don't rescale if the display histogram
            // has no entries.
	    if ( rnorm > 0.0 && norm > 0.0) scale = norm / rnorm;
	    ref1->Scale(scale);

  	    // ECC - 12 Jan 2007 Rescale 2nd ref hist.
            if (ref2) {
              r2norm=ref2->Integral();
              double scale2;
  	      if ( r2norm > 0.0 && norm > 0.0) scale2 = norm/r2norm;
              ref2->Scale(scale2);
	    }
	  }

	  // Draw histograms.
          // ECC - 24 Apr 2003 - Allow one to set limits from .cfg file.
  	  // ECC - 31 Mar 2003 - determine maximum y-axis for plot.
          // ECC - 06 Jul 2005 - Allow y max to float if it is set to zero.
          float h_max = hst->GetMaximum();
          float r_max = ref1->GetMaximum();
          float h_min = 0.0;

          // Do automatic rescaling.
          if (option_y_max == 0.0) {
            if (h_max > r_max)
              option_y_max = 1.1 * h_max;
            else
              option_y_max = 1.1 * r_max;
          }
          hst->SetMaximum(option_y_max);

	  // Scaling overrides the option_y_min feature.
          if (option_scale_ymin != 0.0) {
            h_min = -option_y_max * option_scale_ymin;
	  }
          else if (option_y_min != 0.0) h_min = option_y_min;
          if (!option_logy) hst->SetMinimum(h_min);

	  // Draw histograms
          hst->Draw(draw_opt.c_str());
          ref1->Draw("SAME");
          if (ref2) ref2->Draw("SAME");

  	  // Draw second time to have data on top of reference.
          if (!option_selected) {
            string new_opt = draw_opt + " SAME";
            hst->Draw(new_opt.c_str());
	  }

        }
        // This is a 2D plot.
        else {
	  gStyle->SetOptStat(0);

  	  hst->Draw(draw_opt.c_str());
	}
      }
      // No reference histogram. Just plot it.
      else {
        // ECC - 24 Apr 2003 - Set limits by hand.
        // ECC - 06 Jul 2005 - Allow one to only set the minimum.
        float h_min = 0.0;
        float h_max = hst->GetMaximum()*1.1;
        if (option_y_max != 0.0) {
          h_max = option_y_max;
          hst->SetMaximum(option_y_max);
        }
        if (option_scale_ymin != 0.0) h_min = -option_scale_ymin*h_max;
        else if (option_y_min != 0.0) h_min = option_y_min;
        if (!option_logy) hst->SetMinimum(h_min);

	// GAF - 7 Oct 2010 - Allow modification of stats box when no ref
	gStyle->SetOptStat(option_stat_value);

        hst->Draw(draw_opt.c_str());
      }


      // ECC - 12 May 2003 - Write run number onto plot.
      if (thePage->ShowRunNumber && hist == 0) {
        TLatex l;
        char run_label[15];
        l.SetTextAlign(11);  // bottom left.
        l.SetTextSize(0.07); // Set size relative to pad size.
        l.SetNDC();          // Use universal coordinates
        // GAF - 24 Aug 2010 - separate run and subrun number 
        //sprintf(run_label,"Run: %6.6d", fRunNumber); // make label
        sprintf(run_label,"Run: %8.8d SubRun: %4.4d",fRun, fSubRun); // make label
        l.DrawLatex(0.1, 0.0, run_label); // Show run number.
      }

      // ECC - try to associate function with this plot.
      // ECC - associate new config file with histogram.
      if ( config_file.size() != 0) {
      // See whether we want to associate an exec with this plot.
        gPad->SetTitle(config_file.c_str());

      }

      // If the statistics box (or anything specific to the style has
      // been changed for this pad, then one has to update it. However,
      // this is not ideal since it slows the refresh.
      if (option_stat_value > 0) {
        gPad->Modified();
        gPad->Update();
      }

      // Don't forget to increment pad number
      pad++;
      if ( pad > npads ) break;

    }

  canvas->Modified();
  canvas->Update();

  if ( fDebug)
    Info("DrawPage","END");
}

//---------------------------------------------------------------
// Set default page style.
//---------------------------------------------------------------
void GMPlot::PageStyle(TCanvas *canvas, GMPage *thePage)
{
  // Set canvas attributes
  if ( fDebug) Info("DrawPage","Set canvas attributes");
  canvas->Clear();

  // Use the gmbrowser style.
  //  GMStyle->cd();
  gROOT->SetStyle("Plain");

  // Need to specify this since it seems that it gets written over...
  gStyle->SetTitleW(0);       // width  of title-box
  gStyle->SetTitleH(0);       // height of title-box

  // Make statistics box transparent.
  gStyle->SetOptStat(1);
  gStyle->SetStatW(K_STAT_WIDTH);
  gStyle->SetStatH(K_STAT_HEIGHT);
  gStyle->SetStatColor(fCanvasColor);
  gStyle->SetStatStyle(0);
  gStyle->SetStatBorderSize(1);


  gROOT->ForceStyle();

  canvas->SetFillColor(0);

  // Percent of canvas to use as whitespace.
  Float_t small = 1e-4;
  canvas->Divide(thePage->DivX, thePage->DivY, small, small);

  //  canvas->SetFillColor(fCanvasColor);
  canvas->SetTitle(thePage->PageName.c_str());

}

//---------------------------------------------------------------
// Cause windows to appear on screen
//---------------------------------------------------------------
void GMPlot::Run(bool LOOP)
{
  if ( fDebug)
    Info("Run","Create Windows");

  bool lcont = kTRUE;
  bool make_plots;
  int tot_evt = 0, prev_tot = 0, ifirst = 1;
  int prev_run = 0;

  // loop continuously.
  
  while(lcont) {
    // Loop over all pages and create a .gif file for each page.
    for (int ipage=0; ipage<fPageList.size(); ipage++) {
      fCurrentPage = ipage;
      DrawPage(&tot_evt);

      // Only generate new .gif files if there are new events.
      if (ipage == 0) {
        if (tot_evt != prev_tot || prev_run != fRunNumber) make_plots = kTRUE;
        else make_plots = kFALSE;

	// However, make plots at least one time.
        if (ifirst) {
          ifirst = 0;
          make_plots = kTRUE;
	}
        prev_tot = tot_evt;
        prev_run = fRunNumber;
      }

      if (make_plots) saveWWW();

      if (!LOOP) lcont = kFALSE;
      else gSystem->Sleep(500);
    }
  }

  // Loop over the pages. Sleep after all pages have been displayed.
}




// Checks whether keys can be accessed from file.
// Returns kFALSE if keys are bad.
Bool_t GMPlot::checkKeys(int K)
{
  // ECC - 1 May 2003 - Try to get keys to be sure that file is readable.
  Int_t nkeys  = 0;
  Int_t ntimes = 0;

  while (!nkeys) {
    nkeys = fRootFile[K]->ReadKeys();
    //    cout << "checkKeys: " << nkeys << endl;

    // Unable to find valid keys. Try again.
    if (!nkeys) {
      ntimes++;
      // cout << "DrawPage: Couldn't find keys. ntimes = " << ntimes << endl;

      // Stop after 10 tries.
      if (ntimes >= 10) {
        Error("DrawPage","Unable to read keys.");
        return kFALSE;
      }
      // Wait a little bit and try again.
      gSystem->Sleep(50);
    }
  }
  return kTRUE;
}

// Decode options for each histogram.
// Options are delimited by spaces.
// Valid options are:
//                    lego
//                    surf
//                    logx
//                    logy
//                    logz
//                    noscale              i.e. turn of automatic scaling
//                    noKolfit             i.e. turn off the Kol fitting.
//                    limits(<ymin>,<ymax>)
//                    stats(bitmask)       i.e. stats(1111) is the root default
void GMPlot::histOptions(string all_options, 
                         string &config_file, TH1 *hst)
{
  string remain, option;
  string value1, value2;

  // Set defaults.
  option_y_min    = 0.0;
  option_y_max    = 0.0;
  option_scale_ymin = 0.0;
  option_stat_value = 0;
  option_lwid     = 0;
  option_selected = 0;
  option_scaleRef = 1;
  option_goodness = 1;
  option_PlotRef  = 1;
  option_PlotRef2 = 1;
  option_logy     = 0;
  option_range_min = 0;
  option_range_max = 0;
  option_rebin     = 0;

  // Get the first option from the string.
  split(all_options, option, remain);
  option = strip(option);

  // Loop over the requested options.
  while (option.size() != 0) {

    // ECC - 24 Apr 2003 - allow log options on plots.
    if (option == "logx") {
      gPad->SetLogx();
    }
    else if (option == "logy") {
      gPad->SetLogy();
      option_logy = 1;
    }
    else if (option == "logz") {
      gPad->SetLogz();
    }
    else if (option == "noscale") {
      option_scaleRef = 0;  // Don't rescale reference histogram.
    }
    else if (option == "noKolfit") {
      option_goodness = 0;  // Don't do Kolmogorov test.
    }
    else if (option == "noRef") {
      option_PlotRef = 0; // do not plot reference histogram
    }
    else if (option == "noRef2") {
      option_PlotRef2 = 0; // do not plot reference histogram
    }
    // Set ymin and ymax
    else if (optionValue(option, "limits", 2, value1, value2)) {
      option_y_min = atoi(value1.c_str());
      option_y_max = atoi(value2.c_str());
      // Protect against weird values.
      if (option_y_min > option_y_max)  option_y_max=0.;
    }
    // Set ymin scaling fraction.
    // i.e. if option_scale_ymin = 0.1, then the minimum is
    // -1*0.1 of the maximum.
    else if (optionValue(option, "scale_ymin", 1, value1, value2)) {
      option_scale_ymin = atof(value1.c_str());
      if (option_scale_ymin < 0. || option_scale_ymin > 1.0)
        option_scale_ymin = 0.0;
    }
    // Set ymin and ymax
    else if (optionValue(option, "range", 2, value1, value2)) {
      option_range_min = atoi(value1.c_str());
      option_range_max = atoi(value2.c_str());
      // Protect against weird values.
      if (option_range_min > option_range_max)  option_range_max=0.;
    }
    // Set rebinning factor.
    else if (optionValue(option, "rebin", 1, value1, value2)) {
      option_rebin = atoi(value1.c_str());
      // Protect against weird values.
      if (option_rebin < 0)  option_rebin=0;
    }
    // Change histogram line width.
    else if (optionValue(option, "lwid", 1, value1, value2)) {
      option_lwid = atoi(value1.c_str());
    }
    // Change statistics box.
    else if (optionValue(option, "stats", 1, value1, value2)) {
      option_stat_value = atoi(value1.c_str());
    }
    // Associate a configuration file with this plot.
    else if (optionValue(option, "config", 1, value1, value2)) {
      config_file = value1;
    }
    else {
      // These are TH1D root options
      hst->SetOption(option.c_str());
      option_selected = 1;
    }

    // Get next option.
    string line = remain;
    if (line.size() == 0) option = line;
    else split(line, option, remain);
  }


}
// Routine to get values from a particular option.
// Format is either "option(value)" or "option(value1:value2)".
int GMPlot::optionValue(string option, string search, int nvalues,
                            string &value1, string &value2)
{
  int ic;
  string line_end;
  ic = option.find(search.c_str());
  // See if search string is contained in option.
  if (ic >= 0) {
    // strip off search string + "(".
    string temp1 = strip(option.substr(ic+search.size()+1));
    split(temp1, value1, line_end, ')');

    // perform another split if we need two values.
    if (nvalues > 1) {
      string low, high;
      low = value1;
      split(low, value1, line_end,':');
      split(line_end, value2, high, ')');
    }
    return 1;
  }
  else return 0;
}




//---------------------------------------------------------------
// Save plots in .gif file for viewing via www.
//---------------------------------------------------------------
void GMPlot::saveWWW(int page)
{
  static int file_open = 0;

  // Temporary .ps file name.
  string tmp_file(fGMPlotPSDir), tmp_name; 
  tmp_file  += "/"; tmp_file += fTitle.c_str(); tmp_file += "_tmp.ps";

  if ( page < 0 ) page = fCurrentPage;

  // Add stuff to store all pages in a single .ps file.
  // Open temporary file if it isn't already open.
  if (!file_open) {
    tmp_name = tmp_file + "(";
    //cout << "Opening file: " << tmp_name << endl;
    file_open = 1;
  }
  else tmp_name = tmp_file;

  // If this is the last page, then close temporary file
  // and move it to the correct name.
  if (page == fPageList.size()-1) {
    //cout << "Last file: page = " << page << endl;
    tmp_name = tmp_file + ")";
    fECanvas->Print(tmp_name.c_str());

    // Change the file name.
    char command_string[256];
    /*sprintf(command_string,"mv %s %s/%s_run%6.6d.ps",
      tmp_file.c_str(), fGMPlotPSDir, fTitle.c_str(), fRunNumber);*/
    sprintf(command_string,"mv %s %s/%s_run%8.8d_%4.4d.ps",
      tmp_file.c_str(),
      fGMPlotPSDir,
      fTitle.c_str(),
      fRun,
      fSubRun);

    // Call mv command.
    system(command_string);
    file_open = 0;
    
    // Compress the file
    sprintf(command_string,"gzip --force %s/%s_run%8.8d_%4.4d.ps",
      fGMPlotPSDir,
      fTitle.c_str(),
      fRun,
      fSubRun);
 
    // Call gzip command
    system(command_string);

  }
  // Just print the page.
  else {
    //    cout << "Printing page: " << page << " " << tmp_name << endl;
    fECanvas->Print(tmp_name.c_str());
  }

  // Generate name for this page.
  string name = makeName(fWWWDir, page);

  // Names for .ps and .gif files.
  string psfile(name.c_str());  psfile += ".eps";
  string giffile(name.c_str()); giffile += ".gif";

  //  sprintf(psfile, "%s.eps", name.c_str());
  //  sprintf(jpgfile,"%s.jpg",name.c_str());

  // Generate .ps file.
  fECanvas->SaveAs(psfile.c_str());

  // Now convert it to a .gif file.
  // have to use system calls to do this conversion...
  //string command("(pstopnm -ppm  -xsize 900 -stdout ");
  string command("(pstopnm -ppm  -xsize 1500 -ysize 1400 -stdout ");
  command += psfile.c_str(); command +=  " | ppmtogif  > "; 
  command += giffile.c_str(); command += ") > /dev/null";
  //  cout << "Command: " << command.c_str() << endl;
  write("Writing file %s", giffile.c_str());
  gSystem->Exec(command.c_str());

  // Delete .eps files.
  string command2("rm "); command2 += psfile.c_str();
  gSystem->Exec(command2.c_str());
}


//---------------------------------------------------------------
// Write to status bar
//---------------------------------------------------------------
void GMPlot::write(const char *fmt, const char *str)
{
  char MSG[512];
  sprintf(MSG, fmt, str);
  cout << MSG << endl;
  //  fStatusBar->SetText(MSG, 1);
}



// Generates plot title prefix from configuration file name.
void GMPlot::plotTitle(const char *cFile)
{

  // Set title from name of configuration file.
  string tmpstr(cFile), ext_str;

  // Strip off the path by searching for slashes.
  int ic1 = 1;
  while (ic1 > -1) {
    ic1 = tmpstr.find('/');
    if (ic1 > -1) tmpstr=tmpstr.substr(ic1+1, ic1+tmpstr.size());
  }

  // Store output in fTitle and get rid of the file extension.
  split(tmpstr, fTitle, ext_str, '.');
  
}


string GMPlot::makeName(const char *dir, int page)
{
  if ( page < 0 ) page = fCurrentPage;

  TString name(Page(page)->PageName.c_str());
  name = name.ReplaceAll(" ","_");
  name = name.ReplaceAll("/","_");  
  name = name.ReplaceAll(":","_");  
  name = name.ReplaceAll("(","_");  
  name = name.ReplaceAll(")","_");  

  char file[256];
  //sprintf(file,"%s/%s%2.2d_%s", dir, fTitle.c_str(), page, name.Data());
  sprintf(file,"%s/%s_run%8.8d_%4.4d_%2.2d", dir, fTitle.c_str(), fRun, fSubRun,page);
  return string(&file[0]);
}

//---------------------------------------------------------------
// Adds names of histogram files so that multiple files can be
// summed together.
//---------------------------------------------------------------
void GMPlot::addfile(string &name, int flag) {
  if (!flag) {
    // Increment the number of TFile vectors.
    int nFiles = fsumFiles.size() + 1;

    // Resize the appropriate vectors of TFiles and Names.
    fsumFiles.resize(nFiles);
    fsumName.resize(nFiles);

    // Store the number of TFiles
    total_sumFiles = fsumFiles.size();
  }
  else {
    // Get full name
    std::string filename = StrDup(gSystem->ExpandPathName(name.c_str()));

    // Open root file.
    TFile * rootfile = TFile::Open(filename.c_str(), "readonly");

    // Store in our vector.
    fsumFiles[total_sumFiles-1].push_back(rootfile);

    // Store the filename too.
    fsumName[total_sumFiles-1].push_back(name);
  }
}

// Sum up all histogram files and return pointer to histogram.
TH1* GMPlot::SumHists( int sumfile, TH1* hin, const char *name) {

  TH1* hist = hin;

  // Make sure errors are correct.
  hist->Sumw2();

  // loop over all available histograms.
  int nhists = fsumFiles[sumfile-1].size();

  for (int ihist = 0; ihist < nhists; ihist++) {

    // Determine file to add.
    TFile *file = fsumFiles[sumfile-1][ihist];
    if (!file) continue;

    if (file->IsOpen()) {

      // Read keys for this file.
      int nkeys = file->ReadKeys();

      // Check to see whether file is okay.
      if (nkeys) {

	// Get histogram from the file.
        TH1* addfile = (TH1*) file->Get(name);

	// Check that the histogram is valid
        if (addfile && (addfile->TestBit(TObject::kNotDeleted)) ) {
          hist->Add(addfile, 1.);
        }
      } // if (nkeys)
    }   // if (file)
  }     // for (ihist
  return hist;
}

//---------------------------------------------------------------
// Closes the current file (if open) and tries to open new file.
//---------------------------------------------------------------
void GMPlot::load(string &name, int file, int *index)
{
  // Close current version of file.
  if (fRootFile[file]) {
    delete fRootFile[file];
    fRootFile[file] = 0;
  }

  // Store the file name in the vector fRootFilename.
  fRootFilename.push_back(StrDup(gSystem->ExpandPathName(name.c_str())));

  // Return the index of this entry.
  *index = fRootFilename.size()-1;
  //  std::cout << "Expanded name: " << fRootFilename[*index].c_str() << std::endl;


  if ( fDebug )
    Info("load", "Trying to open root file<%s>", 
	 fRootFilename[*index].c_str());

  // ECC - 13 Mar 2003 - use generic open
  fRootFile[file] = TFile::Open(fRootFilename[*index].c_str(),"readonly"); 
  if ( ! fRootFile[file] ) {
      Error("load","Unable to open file <%s>", fRootFilename[*index].c_str()); 
      return;
  }

  if ( fDebug )
    Info("load", "Opened root file<%s>", fRootFilename[*index].c_str());
  write("Opened root file %s", fRootFilename[*index].c_str());

  // ECC - 1 May 2003 - Check keys to see if file is readable.
  // ECC - 29 Jan 2007 - Fix indexing bug.
  if (!checkKeys(file)) {
    cout << "Problem reading keys." << endl;
    return;
  }

  if ( fDebug )
    Info("load", "Done");
}

//---------------------------------------------------------------
// Clean up after us!
//---------------------------------------------------------------
GMPlot::~GMPlot()
{

}

