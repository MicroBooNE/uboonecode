//////////////////////////////////////////////////////////////////////////////
// File:         GM.cpp
// Description:  Relatively simple histogram browser for Global Monitor
// Created:      August 2002 Harrison B. Prosper, Pushpa Bhat
//               09-Oct-2002 Version 3.0
//               04-Dec-2002 HBP 
//                           1. Combine "Print All" output into a single
//                              postscript file
//                           2. Add Draw option field in config file
//               10-Dec-2002 HBP Add SaveAll
//               15-Dec-2002 HBP Add TPaveLabel to canvas
//               20-Dec-2002 HBP Update run number each time a page is
//                               drawn
//                               Add auto-save
//               26-Mar-2003 ECC Add in code to read .root files over Ethernet.
//                               Plot histo/ref depending upon which has the
//                                 most entries.
//                               Display message screen when autosaving.
//               31-Mar-2003 ECC Plot data before reference histogram, but
//                               set y maximum from greater of the two.
//                1-Apr-2003 ECC For PhysEx use histogram to get run number.
//                               Modify size and placement of title.
//               15-Apr-2003 ECC Force statistics box to be transparent.
//               24-Apr-2003 ECC Allow one to set logy, logx and logz.
//                               Allow one to force limits on plot.
//                               Remove the 80 character limit in .cfg file.
//               25-Apr-2003 ECC Allow gmbrowser to access more up to 5 .root
//                                 files.
//                               Add Debug.On: directive to .cfg file.
//                               Set gErrorIgnoreLevel to suppress warnings.
//               30-Apr-2003 ECC Include noscal option so that reference 
//                                 histograms are not rescaled.
//                1-May-2003 ECC Check ReadKeys before trying to get histogram.
//                               Check limits on display histogram before 
//                                 rescaling.
//               15-May-2003 ECC Add in code to reread keys before getting
//                                 histogram.
//                               Set file index to 0 when file is closed.
//               04-Jun-2003 ECC Draw with error bars even when no reference
//                                 histogram exists.
//               01-Aug-2003 ECC Add code to protect against histogram read
//                                 errors
//               06-Aug-2003 ECC Change histogram delimiter from "/" to "|".
//                               Plot reference histogram first so that it
//                                 doesn't obscure the data histogram.
//                               Add in capability to run a macro. See
//                                 Page.Macro
//               20-Aug-2003 ECC Fix plotting of histogram. Previous fix
//                                 screwed up the stats box.
//                               Add option stats to change statistics plotting
//                               Make reference histogram thinner.
//               21-Aug-2003 ECC Turn off updating (i.e. refreshing the current
//                                 page) if cycling is also on. And turn off
//                                 cycling if the updating button is pushed.
//                               Make pad transparent for stats option 
//                                 otherwise autosave plot gets screwed up.
//               26-Aug-2003 ECC Modify how saveAll works so that plots with
//                                 different statistics boxes look all right.
//               27-Aug-2003 ECC Add SavePS so that the .ps file is generated
//                                 whenever the "Cycle" button is highlighted.
//                               Remove the Autosave feature since this is
//                                 now incorporated in the Cycle feature.
//               02-Sep-2003 ECC Small bug fix to SaveAll function.
//                               Fix print() command.
//               03-Oct-2003 ECC Remove run number from hist title if it
//                               is <= 0.
//               08-Oct-2003 ECC Fix icon directory lookup.
//               24-Oct-2003 ECC Remove fHistMap since it really isn't
//                                 necessary.
//                               Remove some code for determining the run 
//                                 number.
//               04-Nov-2003 ECC Create PageStyle method for initializing
//                                 style for each page.
//                               Close all open files when a new configuration
//                                 file is read.
//               11-Nov-2003 ECC Add in code to allow one to resize page index.
//                               Add in lwid option for line width.
//               19-Nov-2003 ECC Change from .gif to .jpg files.
//                               Change name of files from GMPagesXXX to
//                               <title>XXX.
//               20-Nov-2003 ECC Go back to .gif files since .jpg files are
//                                 not supported by root. 
//                               Use .cfg file name to generate <title>.
//                               Only write out WWW and .ps plots if 
//                                 Cycle.Plots: 1 option is chosen.
//               22-Feb-2004 ECC Bug fix in lwid parsing.
//                               Allow histogram to be linked to a new
//                                 configuration file.
//                               Add button to return to previous config
//                                 file.
//                               Try to fix problem with running out of
//                                 file descriptors. Modify Close() to delete.
//               03-Mar-2004 ECC Reduce space around pads.
//                               Get rid of some extraneous messages.
//                               Change default style to plain.
//               09-Mar-2004 ECC Define default style to get rid of title 
//                                 color.
//                               Translate environment variable in load 
//                                 directory.
//               12-May-2004 ECC Clean up histOptions.
//                               Add option (Page.ShowRunNumber) to allow
//                                 one to display the run number once on a page
//                               Add ability (noRef) to turn off reference
//                                 histograms.
//               21-Jun-2004 ECC Add button for plot descriptions.
//               23-Jun-2004 ECC Add total number of events to status bar.
//               16-Dec-2004 ECC Change Save button to save as .gif file.
//               14-Jan-2005 ECC Add Cycle.PSPlots option to save .ps plots
//                               during cycling.
//               28-Jan-2005 ECC In Update mode, only create www plots if
//                               fCyclePlots is set.
//               08-Apr-2005 ECC Add SaveAll (.gif) so that all plots can
//                               be saved at once.
//               06-Jul-2005 ECC Allow y min to be fixed while allowing y max to float.
//               19-Aug-2005 ECC Fix bug with logy introduced last time.
//               22-Sep-2005 ECC Force saveGIF to print in portrait mode.
//               22-Dec-2005 ECC Add feature to allow rebinning.
//               03-Feb-2006 ECC Changed from E to E0 for histogram type.
//               24-Mar-2006 ECC Add ability to add together histograms from
//                               multiple files.
//               03-Nov-2006 ECC Initialize the number of summed files to 1.
//               11-Jan-2007 ECC Rewrite histogram file handling.
//                               Add code for 2nd reference histogram.
//               29-Jan-2007 ECC Bug fix in GMBrowser::load.
//               24 Aug 2010 GAF Separate run and subrun number.
//               07 Oct 2010 GAF Support for TH1D, TH2D, THProfile and THProfile2D.
//                               Show run and subrun number in status bar instead 
//                               of histogram title.
//////////////////////////////////////////////////////////////////////////////
//$Revision: 1.9 $

#include "gmbrowser/GMBrowser.hpp"

// Since there is no possibility of conflict, we expose the full std namespace.

using namespace std;

#include "TPad.h"
#include "TClass.h"

#include <cassert>

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
long int GMBrowser::getRunNumber(const char *file, TFile *rootfile)
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
    long int irun = (int)physhist->GetBinContent(1);
    if (irun) return irun;
    cout << " run number is " << irun << endl;
  }
  else return 0;
}

//---------------------------------------------------------------
// Get number of events from a specific histogram.
//---------------------------------------------------------------
int GMBrowser::getTotalEvents(TFile *rootfile)
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

string VERSION = "v3.21.01,  29 January 2007";

const char *K_ABOUT = 
"\n"
"v3.20.00 24 March 2006   Dzero Global Monitor\n"
"Harrison B. Prosper      harry@fnal.gov\n"
"Pushpalatha Bhat         bhat@fnal.gov\n"
"Elliott Cheu             echeu@fnal.gov\n";

const char *K_HELP = 
"\n"
"See $GMBROWSER_DIR/gmbrowser/doc/README.gmbrowser\n";

// E. Cheu - Used for determining histogram file types
enum HistFiles {
  K_DATA_FILE,
  K_REF1_FILE,
  K_REF2_FILE
};

// Each MenuItem corresponds to code that is called when
// that item is activated by the user.

enum MenuItem
{
  M_FILE_OPEN_FILE,

  M_FILE_SAVEAS,
  M_FILE_SAVEAS_PS,
  M_FILE_SAVEAS_EPS,
  M_FILE_SAVEAS_GIF,
  M_FILE_SAVE_ALL,
  M_FILE_SAVE_ALL_GIF,

  M_FILE_PRINT,
  M_FILE_PRINT_ALL,

  M_FILE_EXIT,

  M_HELP_CONTENTS,
  M_HELP_ABOUT
};

// Likewise for the buttons
enum Button
{
  B_OPEN,
  B_SAVE,
  B_PRINT,
  B_UPDATE,
  B_CYCLE,
  B_CONFIG,
  B_PLOTINFO
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

const int K_LISTBOX          = 1;
const int K_MENUBAR_HEIGHT   = 26;    // pixels
const int K_TOOLBAR_HEIGHT   = 32;    // pixels
const int K_STATUSBAR_HEIGHT = 26;    // pixels
const int K_LISTBOX_WIDTH    = 150;   // pixels

const int K_PORTRAIT         = 111;   // TPostScript option
const int K_LANDSCAPE        = 112;
const int K_EPS              = 113;


// Use in Open File menu item to specify which file types should be
// displayed.

const char *OpenFileTypes[] = { "Configuration files", "*.cfg",
				  "All files",           "*",
				  0,               0 };

// Displays available file formats

const char *SaveAsTypes[] = { "PostScript",   "*.ps",
			      "Encapsulated PostScript", "*.eps",
			      "Gif files",    "*.gif",
			      "All files",    "*",
			      0,              0 };

/*
  Structure of GUI 
  ----------------

  MainFrame = [Menubar,
	       VerticalFrame=
	          [ToolBar,
		   HorizontalFramce=[ListBox, Canvas]
		   StatusBar
		  ]
	      ]
*/


//---------------------------------------------------------------
// Class to read from .GMBrowser file in $HOME
//---------------------------------------------------------------
GMEnv::GMEnv()
{
  string gmfile( string(gSystem->ExpandPathName("$HOME")) + "/.GMBrowser");
  ifstream *input = new ifstream(gmfile.c_str());
  if ( !input || (input && (input->fail() || input->bad()) ) )
    {
      ofstream out(gmfile.c_str());
      out << "Created:\t" << getTime() << endl;
      out.close();
    }
  else
    input->close();
}

GMEnv::~GMEnv() {}

const char *GMEnv::Get(const char *var, const char *defValue)
{
  // Open hidden file $HOME/.GMBrowser

  string gmfile( string(gSystem->ExpandPathName("$HOME")) + "/.GMBrowser");
  ifstream *input = new ifstream(gmfile.c_str());
  if ( !input || (input && (input->fail() || input->bad()) ) )
    return defValue;

  // Return value of specified variable

  
  char line[fLineSize];
  while ( input->getline(line, fLineSize) )
    {
      string keyword, value;

      string str(strip(line));         // Strip away leading and trailing space
      split(str, keyword, value, '#'); // Strip away comments
      str = keyword;
      if ( str.size() == 0 ) continue; // Ignore blank lines

      split(str, keyword, value);
      keyword = strip(keyword);
      // Remove trailing ":"
      keyword = keyword.substr(0,keyword.size()-1);
      value   = strip(value);
      if ( string(var) == keyword )
	{
	  input->close();
	  return value.c_str();
	}
    }
  input->close();
  return defValue;
}

void GMEnv::Set(const char *var, const char *newValue)
{
  // Open hidden file $HOME/.GMBrowser

  string gmfile( string(gSystem->ExpandPathName("$HOME")) + "/.GMBrowser");
  ifstream *input = new ifstream(gmfile.c_str());
  if ( !input || (input && (input->fail() || input->bad()) ) )
    {
      ofstream out(gmfile.c_str());
      out << "Created:\t" << getTime() << endl;
      out.close();
      input = new ifstream(gmfile.c_str());
    }

  string gmtmpfile(gmfile + ".tmp");
  ofstream out(gmtmpfile.c_str());

  // Return value of specified variable

  char line[fLineSize];
  while ( input->getline(line, fLineSize) )
    {
      string keyword, value, comment;

      string str(strip(line));         // Strip away leading and trailing space
      split(str, keyword,comment, '#');// Strip away comments
      str = keyword;
      if ( str.size() == 0 ) continue; // Ignore blank lines

      split(str, keyword, value);
      keyword = strip(keyword);
      // Remove trailing ":"
      keyword = keyword.substr(0,keyword.size()-1);
      value   = strip(value);
      if ( string(var) == keyword )
	{
	  if ( comment != "" ) comment = string(" # ") + comment;
	  string newline(keyword+":\t"+string(newValue)+comment);
	  out << newline << endl;
	}
      else
	out << line << endl;
    }
  input->close();
  out.close();
  string cmd("mv "+gmtmpfile+" "+gmfile);
  system(cmd.c_str());
}


//---------------------------------------------------------------
// Constructor
//---------------------------------------------------------------
GMBrowser::GMBrowser(const char *name, UInt_t width, UInt_t height,
		     Bool_t debug)
  : TGMainFrame     (gClient->GetRoot(), width, height),
  fParent         (gClient->GetRoot()), // Pointer to Window Manager
  fDebug          (debug),
  fUpdate         (kFALSE),
  fCycle          (kFALSE),
  fCyclePlots     (kTRUE),
  fCyclePSPlots   (kTRUE),
  fCallPrintDialog(kTRUE),
  fLastPage       (-1),
  fUpdatePeriod   (5000),
  fCyclePeriod    (5000),
  fCurrentPage    (0),
  fHelpText       (K_HELP),
  fRunNumber      (0),
  fRun            (0),
  fSubRun         (0),
  fTotalEventHist (""),
  fRunNumberHist  (""),
  fAutoSavePeriod (300000)
{
  // Initialize the error level.
  gErrorIgnoreLevel = 1001;

  // Initialize some internal variables

  fPageList.clear();
  fWWWDir = StrDup(".");
  fPad = 0;
  // ECC - 25 Apr 2003 - need to loop over all available files.
  for (int i=0; i<K_MAX_FILE; i++) {
    fRootFile[i] = 0;
  }
  fRootFilename.clear();

  // ECC - 8 Oct 2003 - First look in $GMBROWSER_DIR/icons
  // Tool to manage icons. Assume icon files are in the directory
  // specified by $GM_ICON_DIR

  fPicturePool = 0;
  
  const char *scriptsDir = getenv("GMBROWSER_DIR");
  if ( scriptsDir )
    fPicturePool = new TGPicturePool(gClient, "$GMBROWSER_DIR/icons");
  else {
    scriptsDir = getenv("GM_ICON_DIR");
    if ( scriptsDir )
      fPicturePool = new TGPicturePool(gClient, "$GM_ICON_DIR");
    else
      //Warning("GMBrowser", "Can't find icon directory $GM_ICON_DIR");
      cout << "GMBrowser: GM_ICON_DIR not defined" << endl;
  }

  // Get colors

  fBackgroundColor = GetDefaultFrameBackground();
  fClient->GetColorByName("red",    fRed);
  fClient->GetColorByName("yellow", fYellow);
  fClient->GetColorByName("green",  fGreen);
  fClient->GetColorByName("white",  fWhite);

  // Get preferences from $HOME/.GMBrowser

  fGMEnv    = new GMEnv();
  fPrinterCommand = StrDup(fGMEnv->Get("Printer.Command","lpr -P"));
  fPrinter        = StrDup(fGMEnv->Get("Printer","ps5"));
  fCanvasColor    = atoi(fGMEnv->Get("Canvas.Color","42"));

  // Check whether to start update immediately

  fAutoSave = kFALSE;
  
  if ( string(fGMEnv->Get("Update","no")) == "yes" )
    fUpdate = kTRUE;
  else
    fUpdate = kFALSE;
  
  if ( string(fGMEnv->Get("Cycle","no")) == "yes" )
    fCycle = kTRUE;
  else
    fCycle = kFALSE;

  if ( fDebug )
    {
      Info("GMEnv","Printer.Command<%s>", fPrinterCommand);
      Info("GMEnv","Printer<%s>",  fPrinter);
    }
  
  // Build Interface
  //////////////////

  // E. Cheu - Add in TList.
  fWidgets = new TList;

  // Create Menu bar - Child of main frame

  if ( fDebug )
    Info("GMBrowser", "Create Menu Bar");
 
  fMenuBar = new TGMenuBar(this, 1, 1);
  // Args: hints, padleft, padright, padtop, padbottom
  fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft);
  fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);
  fMenuBarLayout = new TGLayoutHints(kLHintsTop |kLHintsLeft | kLHintsExpandX);
  fWidgets->Add(fMenuBarItemLayout);
  fWidgets->Add(fMenuBarHelpLayout);
  fWidgets->Add(fMenuBarLayout);
  AddFrame(fMenuBar, fMenuBarLayout);


  // File Menu - Child of window manager, even though it is inserted into
  // menu bar...a bit confusing.

  if ( fDebug )
    Info("GMBrowser", "Create File Menu");

  fMenuFile = new TGPopupMenu(fParent);
  fMenuFile->AddEntry("&Load Configuration",M_FILE_OPEN_FILE);
  fMenuFile->AddSeparator();

  fMenuFile->AddEntry("Save As...",      M_FILE_SAVEAS);
  fMenuFile->AddEntry("Save As ps",      M_FILE_SAVEAS_PS);
  fMenuFile->AddEntry("Save As eps",     M_FILE_SAVEAS_EPS);
  fMenuFile->AddEntry("Save As gif",     M_FILE_SAVEAS_GIF);
  fMenuFile->AddEntry("Save All (.ps)",  M_FILE_SAVE_ALL);
  fMenuFile->AddEntry("Save All (.gif)", M_FILE_SAVE_ALL_GIF);
  fMenuFile->AddSeparator();

  fMenuFile->AddEntry("&Print",          M_FILE_PRINT);
  fMenuFile->AddEntry("Print &All",      M_FILE_PRINT_ALL);
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry("E&xit",           M_FILE_EXIT);
  fMenuFile->Associate(this);
  fMenuBar->AddPopup("&File", fMenuFile,fMenuBarItemLayout);

  // Help menu

  if ( fDebug )
    Info("GMBrowser", "Create Help Menu");

  fMenuHelp = new TGPopupMenu(fParent);
  fMenuHelp->AddEntry("Contents",  M_HELP_CONTENTS);
  fMenuHelp->AddSeparator();
  fMenuHelp->AddEntry("About",     M_HELP_ABOUT);
  fMenuHelp->Associate(this);
  fMenuBar->AddPopup("&Help", fMenuHelp,fMenuBarHelpLayout);

  // Create a vertical frame to hold toolbar, horizontal frame and status bar

  fVFrameLayout = new TGLayoutHints(kLHintsTop | 
				    kLHintsExpandX | 
				    kLHintsExpandY);
  fWidgets->Add(fVFrameLayout);
  fVFrame = new TGVerticalFrame(this, 1, 1);
  AddFrame(fVFrame, fVFrameLayout);

 // Create Tool bar with buttons and a progress bar

  if ( fDebug )
    Info("GMBrowser", "Create Tool Bar");

  fToolBar = new TGHorizontalFrame(fVFrame, 1, 1);
  fToolBarLayout = new TGLayoutHints(kLHintsTop | kLHintsExpandX);
  fWidgets->Add(fToolBarLayout);
  fVFrame->AddFrame(fToolBar, fToolBarLayout);

  // Create picture buttons

  fButtonLayout = new TGLayoutHints(kLHintsLeft, 0, 2, 0, 2);
  fWidgets->Add(fButtonLayout);

  if ( fPicturePool )
    {

      fOpenIcon    = fPicturePool->GetPicture("file.xpm");
      fOpenButton  = new TGPictureButton(fToolBar, fOpenIcon, B_OPEN); 
      fOpenButton->SetToolTipText("Open configuration file");
      fOpenButton->Associate(this);
      fToolBar->AddFrame(fOpenButton, fButtonLayout);

      fSaveIcon    = fPicturePool->GetPicture("disk.xpm");
      fSaveButton  = new TGPictureButton(fToolBar, fSaveIcon, B_SAVE); 
      fSaveButton->SetToolTipText("Save current page");
      fSaveButton->Associate(this);
      fToolBar->AddFrame(fSaveButton, fButtonLayout);

      fPrintIcon   = fPicturePool->GetPicture("printer.xpm");
      fPrintButton = new TGPictureButton(fToolBar, fPrintIcon,B_PRINT); 
      fPrintButton->SetToolTipText("Print current page");
      fPrintButton->Associate(this);
      fToolBar->AddFrame(fPrintButton, fButtonLayout);
    }

  if ( fDebug )
    Info("GMBrowser", "Create Progress Bar");

  // Progress Bar

  fProgress = new TGHProgressBar(fToolBar, TGProgressBar::kFancy, 1);
  fProgress->SetBarColor("green");
  fProgressLayout = new TGLayoutHints(kLHintsTop   | 
				      kLHintsLeft | 
				      kLHintsExpandX);
  fWidgets->Add(fProgressLayout);
  fToolBar->AddFrame(fProgress, fProgressLayout);

  // Control buttons

  // ECC - June 21, 2004- new button for plot description.
  fPlotButton=new TGTextButton(fToolBar,"Plot Info", B_PLOTINFO); 
  fPlotButton->SetToolTipText("Provide information for current plot.");
  fPlotButton->Associate(this);
  fToolBar->AddFrame(fPlotButton, fButtonLayout);

  // ECC - new button to return to previous configuration.
  fConfigButton=new TGTextButton(fToolBar,"Prev Config", B_CONFIG); 
  fConfigButton->SetToolTipText("Return to previous configuration");
  fConfigButton->Associate(this);
  fToolBar->AddFrame(fConfigButton, fButtonLayout);

  fUpdateButton=new TGTextButton(fToolBar,"Update", B_UPDATE); 
  fUpdateButton->SetToolTipText("Start/Stop periodic update of canvas");
  fUpdateButton->Associate(this);
  fToolBar->AddFrame(fUpdateButton, fButtonLayout);

  fCycleButton = new TGTextButton(fToolBar, "Cycle", B_CYCLE);
  fCycleButton->SetToolTipText("Start/Stop cycling through canvases");
  fCycleButton->Associate(this);
  fToolBar->AddFrame(fCycleButton, fButtonLayout);

  // Create horizontal frame to hold list box and canvas

  fHFrameLayout = new TGLayoutHints(kLHintsTop | 
				    kLHintsExpandX | 
				    kLHintsExpandY);
  fWidgets->Add(fHFrameLayout);
  fHFrame = new TGHorizontalFrame(fVFrame, 1, 1);
  fVFrame->AddFrame(fHFrame, fHFrameLayout);

  // Create Frames.
  fV1 = new TGVerticalFrame(fHFrame, K_LISTBOX_WIDTH, 1, kFixedWidth);
  fV2 = new TGVerticalFrame(fHFrame, 1, 1);

  // Create List box.
  fListBox = new TGListBox(fV1, K_LISTBOX, kSunkenFrame);
  fListBox->Resize(K_LISTBOX_WIDTH, 1);
  fListBox->Associate(this);  // allows communication...

  // ListBox layout: have to allow it to expand in both X and Y.
  fListBoxLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | 
                                     kLHintsExpandX | kLHintsExpandY);
  fWidgets->Add(fListBoxLayout);
  fV1->AddFrame(fListBox, fListBoxLayout);

  // E. Cheu - add listbox frame to horizontal frame.
  fV1Layout = new TGLayoutHints(kLHintsLeft | kLHintsExpandY);
  fWidgets->Add(fV1Layout);
  fHFrame->AddFrame(fV1, fV1Layout);

  // E. Cheu - Add splitter which allows one change size of list box.
  TGVSplitter *splitter = new TGVSplitter(fHFrame);
  splitter->SetFrame(fV1, kTRUE);
  splitLayout = new TGLayoutHints(kLHintsLeft | kLHintsExpandY);
  fWidgets->Add(splitter);
  fWidgets->Add(splitLayout);
  fHFrame->AddFrame(splitter, splitLayout);

  // Define layout for canvas.  
  fV2Layout = new TGLayoutHints(kLHintsRight
                         | kLHintsExpandX | kLHintsExpandY);
  fWidgets->Add(fV2Layout);
  fHFrame->AddFrame(fV2, fV2Layout);

  // Create the style.
  //  GMStyle  = new TStyle("GMStyle","Style with no colors/fill areas");

  // ECC - 03 Mar 2004 - make pad frame smaller.
  /*
  GMStyle->SetFrameBorderMode(0);
  GMStyle->SetCanvasBorderMode(0);
  GMStyle->SetCanvasColor(0);
  GMStyle->SetPadBorderMode(0);
  GMStyle->SetPadColor(0);
  GMStyle->SetPalette(1);
  TGaxis::SetMaxDigits(3);

  GMStyle->SetFillColor(0);
  GMStyle->SetLabelSize(0.04,"X");
  GMStyle->SetLabelSize(0.04,"Y");
  GMStyle->SetLabelOffset(0.01,"X");
  GMStyle->SetLabelOffset(0.01,"Y");

  GMStyle->SetTitleColor(10);
  GMStyle->SetTitleFillColor(10);
  GMStyle->SetTitleBorderSize(1);
  */

  // ECC - 15 Apr 2003 - Change style to be transparent.
  /*
  GMStyle->SetStatColor(0);
  GMStyle->SetOptStat(1);
  GMStyle->SetStatW(K_STAT_WIDTH);
  GMStyle->SetStatH(K_STAT_HEIGHT);
  GMStyle->SetStatColor(fCanvasColor);
  GMStyle->SetStatStyle(0);
  GMStyle->SetStatBorderSize(1);
  */
  //set_color_palette();
  gStyle->SetPalette(1);
  // Create a Canvas

  int w = GetWidth() - K_LISTBOX_WIDTH;
  int h = GetHeight() 
    - K_MENUBAR_HEIGHT
    - K_TOOLBAR_HEIGHT
    - K_STATUSBAR_HEIGHT;
  fECanvas = new TRootEmbeddedCanvas("GM", fV2, w, h);

  fCanvasLayout = new TGLayoutHints(kLHintsRight | kLHintsExpandY 
                                                 | kLHintsExpandX);
  fWidgets->Add(fCanvasLayout);
  fV2->AddFrame(fECanvas, fCanvasLayout);

  // Create a label and a pad

  fPaveLabel = new TPaveLabel(0.0,K_PADYFRAC,1.0,1.0,"Global Monitor");
  fPaveLabel->SetBorderSize(0);
  //fPad       = new TPad("Pad","Pad", 0.0, 0.0, 1.0, K_PADYFRAC,0,0,0);

  // Create Status Bar - split 22% vs 78%
  // Split(0) for the clock; split(1) for status 

  if ( fDebug )
    cout << "GMBrowser - Create Status Bar\n";

  fStatusBar = new TGStatusBar(fVFrame, 1, 1);
  int parts[] = {15, 15, 25};
  fStatusBar->SetParts(parts, 3);
  fStatusBarLayout = new TGLayoutHints(kLHintsTop  | 
				       kLHintsLeft | 
				       kLHintsExpandX);
  fWidgets->Add(fStatusBarLayout);
  fVFrame->AddFrame(fStatusBar, fStatusBarLayout);

  
  // Create timers to handle synchronous operations:
  // 1. Clock
  // 2. Updating
  // 3. Cycling
  // 4. AutoSaving

  if ( fDebug )
    Info("GMBrowser", "Create Timers");

  fClock       = new GMTimer(this, &GMBrowser::setTime);
  fUpdateClock = new GMTimer(this, &GMBrowser::update);
  fCycleClock  = new GMTimer(this, &GMBrowser::cycle);
  //  fAutoSaveClock = new GMTimer(this, &GMBrowser::autoSave);

  if ( fDebug )
    Info("GMBrowser", "Set Window Name");

  SetWindowName(name);
  SetIconName(name);

 // connect canvas events to my process
  TCanvas *canvas = Canvas();
  canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 
          "GMBrowser",
          this, "catch_event(Int_t,Int_t,Int_t,TObject*)");
}



// METHODS
//////////

const char *GMBrowser::Version(){return VERSION.c_str();}


//---------------------------------------------------------------
// Help!
//---------------------------------------------------------------
void GMBrowser::AddHelp(const char *text) {fHelpText = string(text);}

//---------------------------------------------------------------
// Return a pointer to the TCanvas 
//---------------------------------------------------------------
TCanvas  *GMBrowser::Canvas() 
{
  if ( fECanvas )
    return fECanvas->GetCanvas();
  else
    return (TCanvas *)0;
}

//---------------------------------------------------------------
// Return a pointer to specified page
//---------------------------------------------------------------
GMPage   *GMBrowser::Page(int pageNumber) 
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
TObject  *GMBrowser::Histogram(const char *name, int which)
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

Bool_t GMBrowser::openFiles(int sumfile)
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


GMPageList *GMBrowser::PageList(){return &fPageList;}


//---------------------------------------------------------------
// Read and parse configuration file and create pages specified
// therein.
//---------------------------------------------------------------
void GMBrowser::Configure(const char *cFile)
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

  // Get rid of current pages
  DeleteAll();

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
	  if ( fRootFile[K_DATA_FILE] == 0 ) CloseWindow();
          else {
  	    // ECC - new version.
	    //          fRunNumber = getRunNumber(value.c_str());
            fRunNumber = getRunNumber(value.c_str(), fRootFile[K_DATA_FILE]);

            // ECC - 25 Apr 2003 - initialize K_REFFILE to -1.
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
	  //	  if ( fDebug)
	  //	    Info("configure","Root.File.Ref<%s>", value.c_str());
	  load(value, K_REF1_FILE, &K_HISTS[K_REF1_FILE]);
	}

      // E. Cheu - 12 Jan 2007 - Add code for second reference histogram.
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

      else if ( keyword == "Cycle.PSPlots:" )
        {
          fCyclePSPlots = atoi(value.c_str());
          if ( fDebug)
            Info("configure","Cycle.PSPlots<%d>", fCycle);
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
	  fListBox->AddEntry(value.c_str(), page);
          
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
    for (int K=0; K<K_MAX_FILE; K++) {
      cout << "K: " << K << " fRootFile[K]: " << fRootFile[K] << endl;
    }
  }

  // ECC - The cycle mode has precedence over Update. So, turn
  //       off updating the current page if the cycle mode is on.
  if (fCycle && fUpdate) fUpdate = kFALSE;

  // Draw first page

  fCurrentPage = 0;
  fListBox->MapSubwindows();
  fListBox->Layout();
  fListBox->Select(fCurrentPage);
  gClient->NeedRedraw(fListBox); // Force immediate redraw

  DrawPage();

  if ( fDebug)
    Info("Configure","END");
}

//---------------------------------------------------------------
// Add new pages to any existing ones and return page number
//---------------------------------------------------------------
int GMBrowser::AddPage(const char *name, int divx, int divy)
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
// Delete specified page
//---------------------------------------------------------------
void GMBrowser::DeletePage(int pageNumber)
{
  if ( fDebug)
    Info("DeletePage","START");

  if ( fDebug )
    Info("DeletePage", "\tPage<%d>", pageNumber);

  int page = 0;
  GMPageList::iterator p = fPageList.begin();
  while (p != fPageList.end())
    {
      if ( page == pageNumber )
	{
	  if ( fDebug )
	    Info("DeletePage", "\tFound Page<%d>", pageNumber);
	  fPageList.erase(p);
	  break;
	}
      page++;
      p++;
    }
  if ( fDebug)
    Info("DeletePage","END");
}

//---------------------------------------------------------------
// Delete all pages
//---------------------------------------------------------------
void GMBrowser::DeleteAll()
{
  if ( fDebug)
    Info("DeleteAll","START");
  
  for(int page = fPageList.size()-1; page >-1; page--)
    fListBox->RemoveEntry(page);
  fListBox->Layout();

  fPageList.clear();
  if ( fDebug )
    Info("DeleteAll","END");
}


//---------------------------------------------------------------
// Open root-file.
// Read histograms for this page one by one, from the root-file, 
// and add each to a different pad until all pads are used up or
// until we run out of histograms.
// Close root-file
//---------------------------------------------------------------
void GMBrowser::DrawPage(int page)
{
  TFile *fDataFile, *fRefFile;

  if ( fDebug)
    Info("DrawPage","START");

  if ( fPageList.size() == 0 ) return;

  // Default is to draw current page

  if ( page < 0 )
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

  // Update run number just in case file has changed
  fRunNumber = getRunNumber(fRootFilename[K_HISTS[K_DATA_FILE]].c_str(), 
                            fRootFile[K_DATA_FILE]);

  // GAF - 24 Aug 2010 - separate run and subrun number 
  //                     fRunNumber=1000*fRun+fSubRun
  fRun    = (int)floor(fRunNumber/10000.);
  fSubRun = (int)(fRunNumber-10000*fRun);

  // Get total number of events processed.
  int TotalEvents = getTotalEvents(fRootFile[K_DATA_FILE]);

  // Update status bar with total number of events.
  // GAF - 7 Oct 2010 - and Run and SubRun number
  if (TotalEvents > 0) {
    char line[80];
    sprintf(line, "  Run: %d    SubRun: %d    Total Events: %d", fRun, fSubRun, TotalEvents);
    fStatusBar->SetText(line, 2);
  }

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
      //hst->SetTitleSize(K_TITLE_SIZE); 
      //hst->SetTitleOffset(K_TITLE_OFFSET);

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

      // ECC - 12 Jan 2007 - do the same for the second histogram.
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
			    (string(hst->ClassName()) == "TH1F") || 
			    (string(hst->ClassName()) == "TH2F") ||
			    (string(hst->ClassName()) == "TProfile2D") || 
	                    (string(hst->ClassName()) == "TH2D")) ) {
	  double prob  = hst->KolmogorovTest(ref1);
	  char title[132];
	  sprintf(title, "%s (Prob: %1.3f) ", hst->GetTitle(), prob);
	  hst->SetTitle(&title[0]);
	}

        // Superimpose reference histograms only for 1-D histograms
        if ( string(hst->ClassName()) == "TH1D" || string(hst->ClassName()) == "TH1F" || string(hst->ClassName()) == "TProfile") {
	   gStyle->SetOptStat(option_stat_value);

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
        ///sprintf(run_label,"Run: %6.6d", fRunNumber); // make label
        sprintf(run_label,"Run: %4.4d SubRun: %4.4d",fRun, fSubRun); // make label
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
void GMBrowser::PageStyle(TCanvas *canvas, GMPage *thePage)
{
  // Set canvas attributes
  if ( fDebug) Info("DrawPage","Set canvas attributes");
  canvas->Clear();

  // Use the gmbrowser style.
  //  GMStyle->cd();
  gROOT->SetStyle("Plain");

  // Need to specify this since it seems that it gets written over...
  gStyle->SetTitleW(0);       // width  of title-box
  gStyle->SetTitleH(0.07);       // height of title-box

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
void GMBrowser::Run()
{
  if ( fDebug)
    Info("MapWindows","Create Windows");

  MapSubwindows();
  Layout();
  MapWindow();
  fClock->Start(1000);

  // Check whether to start auto save immediately
  /*
  if ( fAutoSave )
    {
      fAutoSaveButton->SetBackgroundColor(fGreen);
      fStatusBar->SetText("Auto save on",1);
      fAutoSaveClock->Start(fAutoSavePeriod);
    }
  */
  // Check whether to start updating immediately

  if ( fUpdate )
    {
      fUpdateButton->SetBackgroundColor(fGreen);
      fStatusBar->SetText("Updating...",1);
      fUpdateClock->Start(fUpdatePeriod);
    }

  // Check whether to start cycling immediately

  if ( fCycle )
    {
      fCycleButton->SetBackgroundColor(fGreen);
      fStatusBar->SetText("Cycling...",1);
      fCycleClock->Start(fCyclePeriod);
    }
}


//---------------------------------------------------------------
// Exit browser
//---------------------------------------------------------------
void GMBrowser::CloseWindow()
{
  fClock->Stop();
  gApplication->Terminate(0);
}


//---------------------------------------------------------------
// Decode user interactions
//---------------------------------------------------------------
Bool_t GMBrowser::ProcessMessage(Long_t msg, Long_t id, Long_t parm)
{
  switch (GET_MSG(msg))
    {
    case kC_COMMAND:

      if ( fDebug)
	Info("ProcessMessage","kC_COMMAND<%d>", id);

      switch (GET_SUBMSG(msg))
	{
	case kCM_BUTTON:
	  handleButton(id);
	  break;

	case kCM_MENU:
	  handleMenu(id);
	  break;

	case kCM_LISTBOX:
	  handleListBox(parm);
	  break;

	default:
	  break;
	}

    default:
      break;
    }
  return kTRUE;	  
}

// Checks whether keys can be accessed from file.
// Returns kFALSE if keys are bad.
Bool_t GMBrowser::checkKeys(int K)
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
//                    scale_ymin(fractin)
//                    stats(bitmask)       i.e. stats(1111) is the root default
void GMBrowser::histOptions(string all_options, 
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
int GMBrowser::optionValue(string option, string search, int nvalues,
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
// Handle all menu items of all menus
//---------------------------------------------------------------
void GMBrowser::handleMenu(Int_t id)
{
  if ( fDebug)
    Info("handleMenu","id<%d>", id);

  switch (id) 
    {
    case M_FILE_OPEN_FILE:
      configure();
      break;

    case M_FILE_SAVEAS:
      saveAs();
      break;

    case M_FILE_SAVEAS_PS:
      saveAs(".ps");
      break;
      
    case M_FILE_SAVEAS_EPS:
      saveAs(".eps");
      break;
      
    case M_FILE_SAVEAS_GIF:
      saveAs(".gif");
      break;

      //    case M_FILE_SAVEAS_JPG:
      //      saveAs(".jpg");
      //      break;

    case M_FILE_SAVE_ALL:
      saveAll();
      break;

    // Saves all pages as .gif files.
    case M_FILE_SAVE_ALL_GIF:
      saveGIF();
      break;

    case M_FILE_PRINT:
      print(kTRUE); // Call print dialog and print current page
      break;

    case M_FILE_PRINT_ALL:
      printAll(kTRUE);// Call print dialog and print all pages
      break;

    case M_FILE_EXIT:
      CloseWindow();  
      break;
 
    case M_HELP_CONTENTS:
      {
	TRootHelpDialog *hd = new TRootHelpDialog(fParent, 
						  "Help Contents", 
						  600, 400);
	hd->SetText(fHelpText.c_str());
	hd->Popup();
      }
      break;
      
    case M_HELP_ABOUT:
      {
	TRootHelpDialog *hd = new TRootHelpDialog(fParent, 
						  "About GMBrowser", 
						  375, 100);
	hd->SetText(K_ABOUT);
	hd->Popup();
      }
      break;

    default:
      break;
    }  
}

//---------------------------------------------------------------
// Handle Start/Stop operations
//---------------------------------------------------------------
void GMBrowser::handleButton(Int_t id)
{
  switch (id) 
    {
    // Make help screen for current plot.
    case B_PLOTINFO:
      {
	TRootHelpDialog *hd = new TRootHelpDialog(fParent, 
						  "Help Contents", 
						  600, 400);
        GMPage *thePage = Page(fCurrentPage);

	//        cout << "Help for page: " << fCurrentPage << endl;
	//        cout << "Title: " << thePage->PageName.c_str() << endl;
	//        cout << "File: " << thePage->HelpFile.c_str() << endl;

	// Open file with plot description. Expand any environment variables.
        ifstream *input = 
	  new ifstream(gSystem->ExpandPathName(thePage->HelpFile.c_str()));

        // Title line
        char line[fLineSize];
        sprintf(line, "Information for page: %s", thePage->PageName.c_str());
  	hd->SetText(line);

	// Make sure file is open
        if ( !input || (input && (input->fail() || input->bad()) ) ) {
          sprintf(line,"No plot information available.");
          hd->AddText(line);
	}
        else {
          // Read all lines in file.
          while ( input->getline(line, fLineSize) )  {

	    // Add text to help display.
            hd->AddText(line);
          }
          input->close();

	}

        // Pop up help window.
	hd->Popup();
	// Redraw button.
  	gClient->NeedRedraw(fPlotButton);
      }
      break;
    case B_CONFIG:
      {
	// Return to the previous configuration.
        if (ConfigList.size() > 0) {
          int last = ConfigList.size()-1;
	  //          cout << "ConfigList last: " << last << endl;
	  GMBrowser::configureFile(ConfigList[last].c_str());
          ConfigList.pop_back();

	  // Return to previous page.
          if (PrevPageList.size() > 0) {
            int isize = PrevPageList.size()-1;
            fCurrentPage = PrevPageList[isize];

	    // Clear page off of list.
            PrevPageList.pop_back();

            if (fCurrentPage > 0) {
              // Highlist page name in list box.
              fListBox->Select(fCurrentPage);
              gClient->NeedRedraw(fListBox); // Force immediate redraw

  	      // Draw the page.
              DrawPage();
	    }
	  }
          
	}
	gClient->NeedRedraw(fConfigButton);
      }
      break;
    case B_UPDATE:
      {
	fUpdate = fUpdate ? kFALSE : kTRUE;
	if ( fUpdate )
	  {
	    fUpdateButton->SetBackgroundColor(fGreen);
	    fStatusBar->SetText("Updating...",1);
	    fUpdateClock->Start(fUpdatePeriod);
            
            // ECC - 21 Aug 2003 - turn off cycling.
            if (fCycle) {
              fCycle = kFALSE;
  	      fCycleButton->SetBackgroundColor(fBackgroundColor);
	      fCycleClock->Stop();
    	      gClient->NeedRedraw(fCycleButton);
	    }
	  }
	else
	  {
	    fUpdateButton->SetBackgroundColor(fBackgroundColor);
	    fStatusBar->SetText("Updating OFF",1);
	    fUpdateClock->Stop();
	  }
	gClient->NeedRedraw(fUpdateButton); // Force immediate redraw
	fProgress->Reset();
	fProgressCount = 0;
	gClient->NeedRedraw(fProgress); // Force immediate redraw
      }
      break;

    case B_CYCLE:
      {
	fCycle = fCycle ? kFALSE : kTRUE;
	if ( fCycle )
	  {
	    fCycleButton->SetBackgroundColor(fGreen);
	    fStatusBar->SetText("Cycling...",1);
	    fCycleClock->Start(fCyclePeriod);

	    // ECC - 21 Aug 2003 - turn off updating button if it's on.
            if (fUpdate){
              fUpdate = kFALSE;
  	      fUpdateButton->SetBackgroundColor(fBackgroundColor);
	      fUpdateClock->Stop();
    	      gClient->NeedRedraw(fUpdateButton);
	    }
	  }
	else
	  {
	    fCycleButton->SetBackgroundColor(fBackgroundColor);
	    fStatusBar->SetText("Cycling OFF",1);
	    fCycleClock->Stop();
	  }
	gClient->NeedRedraw(fCycleButton);
      }
      break;

    case B_OPEN:
      configure();
      break;
      
    case B_SAVE:
      saveAs(".gif");
      break;
      
    case B_PRINT:
      {
	print(fCallPrintDialog);//Call print dialog once and print current page
	fCallPrintDialog = kFALSE;
      }
      break;
      
    default:
      break;
    }
}

//---------------------------------------------------------------
// Handle all menu items of all menus
//---------------------------------------------------------------
void GMBrowser::handleListBox(Int_t page)
{
  fCurrentPage = page;
  if ( fDebug)
    Info("handleListBox","page<%d>", page);

  DrawPage(page);
}

//---------------------------------------------------------------
// Save as specified file type
//---------------------------------------------------------------
void GMBrowser::saveAs(const char *type, const char *name)
{
  TString filename("");

  if ( type != 0 )
    {
      TString pagename;
      if ( name != 0 )
	pagename = TString(name);
      else	
	{
	  string name = makeName(fFileSaveInfo.fIniDir, fCurrentPage);
	  pagename  = TString(name.c_str());
	}
      filename = pagename + type;
    }
  else
    {
      fFileSaveInfo.fFileTypes = SaveAsTypes;
      new TGFileDialog(fParent, this, kFDSave, &fFileSaveInfo);
      if (!fFileSaveInfo.fFilename) return;      
      filename = TString(fFileSaveInfo.fFilename);
    }

  write("Saving to file %s", filename.Data());
  if (filename.Contains(".ps")  ||
      filename.Contains(".eps")) {
      //      filename.Contains(".jpg"))
    Canvas()->SaveAs(filename.Data());
  }
  else if (filename.Contains(".gif")){
      gVirtualX->SelectWindow(Canvas()->GetCanvasID());
      gVirtualX->WriteGIF((char *)filename.Data());
    }
  else
    Warning("handleMenu", 
	    "Format for file %s not available", filename.Data());
}

//---------------------------------------------------------------
// configure browser
//---------------------------------------------------------------
void GMBrowser::configure()
{
  fFileLoadInfo.fFileTypes = OpenFileTypes;
  new TGFileDialog(fParent, this, kFDOpen, &fFileLoadInfo);
  TString name(fFileLoadInfo.fFilename);
  if ( name.Contains(".cfg") )
    {
      // Store name of current config file.
      ConfigList.push_back(configFile_curr);
      PrevPageList.push_back(fCurrentPage);

      // ECC - close all open files and reset number of files to -1.
      for (int K=0; K<K_MAX_FILE; K++) {
        if ( fRootFile[K] ) {
          if (fDebug) {
            cout << "Closing File: " << fRootFilename[K_HISTS[K]].c_str() << 
                    " " << fRootFile[K] << endl;
          }
          delete fRootFile[K];
          fRootFile[K] = 0;
        }
      }

      // E. Cheu - 24 Mar 2006 - close all files used in summing.
      printf("Closing files used in summing.\n");
      int nfiles = fsumFiles.size();
      for (int ifile = 0; ifile < nfiles; ifile++) {

        int nhists = fsumFiles[ifile].size();
        for (int ihist = 0; ihist < nhists; ihist++) {

          // Determine file to add.
          TFile *file = fsumFiles[ifile][ihist];
          if (!file) continue;

          // Close the file. This allows gmbrowser to continue to run
          // even if the soft link has changed.
          if (file->IsOpen()) {
            printf("Deleting fsumFiles\n");
            delete file;
	  }
	} // Loop over histogram files.
      }   // Loop over file sets.
      fsumFiles.clear();
      fsumName.clear();


      write("Configuring browser with file %s", name.Data());
      Configure(name.Data());
    }
}

//---------------------------------------------------------------
// configure browser
//---------------------------------------------------------------
void GMBrowser::configureFile(const char *config_file)
{
  TString name(config_file);


  // Only do something if this is a configuration file.
  if ( name.Contains(".cfg") )
    {
      // ECC - close all open files and reset number of files to -1.
      for (int K=0; K<K_MAX_FILE; K++) {
        if ( fRootFile[K] ) {
          if (fDebug) {
            cout << "Closing File: " << fRootFilename[K_HISTS[K]].c_str() << 
                    " " << fRootFile[K] << endl;
          }
          delete fRootFile[K];
          fRootFile[K] = 0;
        }
      }

      // E. Cheu - 24 Mar 2006 - close all files used in summing.
      printf("Closing files used in summing.\n");
      int nfiles = fsumFiles.size();
      for (int ifile = 0; ifile < nfiles; ifile++) {

        int nhists = fsumFiles[ifile].size();
        for (int ihist = 0; ihist < nhists; ihist++) {

          // Determine file to add.
          TFile *file = fsumFiles[ifile][ihist];
          if (!file) continue;

          // Close the file. This allows gmbrowser to continue to run
          // even if the soft link has changed.
          if (file->IsOpen()) {
            printf("Deleting fsumFiles\n");
            delete file;
	  }
	} // Loop over histogram files.
      }   // Loop over file sets.
      fsumFiles.clear();
      fsumName.clear();

      write("Configuring browser with file %s", name.Data());
      Configure(name.Data());
    }
}

//---------------------------------------------------------------
// Write time to split(0) of status bar. Handled by timer fClock 
//---------------------------------------------------------------
void GMBrowser::setTime()
{
  fStatusBar->SetText(getTime().c_str(), 0);
}

//---------------------------------------------------------------
// Auto Save
//---------------------------------------------------------------
void GMBrowser::autoSave()
{
  /*
  fStatusBar->SetText("Auto saving...", 1);
  saveAll();
  fStatusBar->SetText("Saved all pages", 1);
  */
}

//---------------------------------------------------------------
// Update visible page only. No point updating those that are
// hidden! Also, update the progress bar to indicate that things
// are still happening. This is handled by timer fUpdateClock.
//---------------------------------------------------------------
void GMBrowser::update()
{
  fProgressCount++;
  if ( fProgressCount < 101 )
    fProgress->Increment(1);
  else
    {
      fProgressCount = 0;
      fProgress->Reset();
      gClient->NeedRedraw(fProgress); // Force immediate redraw
    }
  DrawPage();
  if (fCyclePlots)   saveWWW();
}

//---------------------------------------------------------------
// Count the number of objects in the given canvas and its subpads,
// not including pads. Used to determine whether a canvas is empty
//---------------------------------------------------------------
int countPadObjects(TPad* p)
{
  int ret=0;
  TIter iter(p->GetListOfPrimitives());
  TObject* obj=0;
  while ((obj = iter())){
    if(obj->IsA()->InheritsFrom(TPad::Class())){
      TPad* subpad=dynamic_cast<TPad*>(obj);
      assert(subpad);
      ret+=countPadObjects(subpad);
    }
    else{
      ++ret;
    }
  }
  return ret;
}

//---------------------------------------------------------------
// Cycle through to the next page. Handled by timer fCycleClock
//---------------------------------------------------------------
void GMBrowser::cycle()
{
  if ( fCurrentPage >= 0 )
    {
      bool foundNonEmptyPage=false;
      while(!foundNonEmptyPage){
        // Increment page to display.
        fCurrentPage++;
        if ( fCurrentPage == (int)fPageList.size() ) fCurrentPage = 0;

        // Highlist page name in list box.
        fListBox->Select(fCurrentPage);
        gClient->NeedRedraw(fListBox); // Force immediate redraw

        // Draw page and make .gif file.
        DrawPage();
 
        if(countPadObjects(Canvas())!=0){
          fListBox->GetEntry(fCurrentPage)->SetBackgroundColor(0xffffff);
          foundNonEmptyPage=true;
        }
        else{
          fListBox->GetEntry(fCurrentPage)->SetBackgroundColor(0xe0e0e0);
        }

        // Only make plots if Cycle.Plots: 1 option chosen. (Default=T)
        if (fCyclePlots)   saveWWW();
        if (fCyclePSPlots) savePS(fCurrentPage);
      }
    }
}

//---------------------------------------------------------------
// Cycle through to the next page. Handled by timer fCycleClock
//---------------------------------------------------------------
void GMBrowser::saveWWW(int page)
{
  if ( page < 0 ) page = fCurrentPage;


  string name = makeName(fWWWDir, page);
  saveAs(".gif", name.c_str());
  write("Saving file %s.gif", name.c_str());
}

//---------------------------------------------------------------
// Save pages to .ps file. This is handled by time fCycleClock.
// Writes to a temporary .ps file and then closes it and moves
// it to its final destination.
//---------------------------------------------------------------
void GMBrowser::savePS(int ipage)
{
  static int file_open = 0;
  //  string tmp_file(fFileSaveInfo.fIniDir); tmp_file  += "/GMPages_tmp.ps";
  string tmp_file(fFileSaveInfo.fIniDir); 
  tmp_file  += "/"; tmp_file += fTitle.c_str(); tmp_file += "_tmp.ps";
  string tmp_name;
  TCanvas *curr_canvas = Canvas();

  // Protect against a weird page number.
  if ( ipage < 0 ) ipage = fCurrentPage;

  // Reinitialize things.
  if (ipage == 0) {

    // Close temporary file and move it to the correct name.
    if (file_open) {
      file_open = 0;
      tmp_name = tmp_file + "]";
      //      cout << "Closing file: " << tmp_name << endl;
      curr_canvas->Print(tmp_name.c_str());

      // Change the file name.
      char command_string[256];
      // GAF - 24 Aug 2010 - separate run and subrun number 
      /*sprintf(command_string,"mv %s %s/%s_run%6.6d.ps",
        tmp_file.c_str(),
        fFileSaveInfo.fIniDir, 
        fTitle.c_str(),
        fRunNumber);*/
      sprintf(command_string,"mv %s %s/%s_run%4.4d_%4.4d.ps",
        tmp_file.c_str(),
        fFileSaveInfo.fIniDir, 
        fTitle.c_str(),
        fRun,
        fSubRun);
     //cout << "Executing command: " << command_string << endl;

      // Call mv command.
      system(command_string);

      // Compress the file
      sprintf(command_string,"gzip -v --force %s/%s_run%4.4d_%4.4d.ps",
        fFileSaveInfo.fIniDir, 
        fTitle.c_str(),
        fRun,
        fSubRun);
 
      // Call gzip command
      system(command_string);

    }
  }


  // Open temporary file if it isn't already open.
  if (!file_open) {
    tmp_name = tmp_file + "[";
    //    cout << "Opening file: " << tmp_name << endl;
    curr_canvas->Print(tmp_name.c_str());
    file_open = 1;
  }

  // Save current canvas to .ps file.
  curr_canvas->Print(tmp_file.c_str());
}

//---------------------------------------------------------------
// Write to status bar
//---------------------------------------------------------------
void GMBrowser::write(const char *fmt, const char *str)
{
  char MSG[512];
  sprintf(MSG, fmt, str);
  fStatusBar->SetText(MSG, 1);
}


//---------------------------------------------------------------
// Print current page
//---------------------------------------------------------------
void GMBrowser::print(Bool_t callDialog)
{
  if ( fDebug )
    Info("print","Entered");

  if ( fPageList.size() == 0 ) return;

  // Check whether to call print dialog

  if ( callDialog )
    {
      Int_t OK = kTRUE;
      new TGPrintDialog(fParent, this, 600, 150, 
			&fPrinter, &fPrinterCommand, &OK); 
      if ( ! OK ) return;
    }

  char filename[160];
  sprintf(filename,"%s/GMpage%2.2d.ps", fFileSaveInfo.fIniDir, 
          fCurrentPage);
            
  if ( fDebug )
    Info("print","file<%s>", filename);
		  
  // ECC - Use canvas->Print instead of TPostScript.
  TCanvas *curr_canvas = Canvas();
  curr_canvas->Print(filename);

  // Send print command to shell

  string command(string(fPrinterCommand)+" "+
                 string(fPrinter)+" "+
                 string(filename)+" &");
  write("%s",command.c_str());
  system(command.c_str());

  if ( fDebug )
    Info("print","Exit");
}


//---------------------------------------------------------------
// Print all pages to a single postscript file
//---------------------------------------------------------------
void GMBrowser::printAll(Bool_t callDialog)
{
  if ( fDebug )
    Info("printAll","Entered");

  if ( fPageList.size() == 0 ) return;

  if ( callDialog )
    {
      Int_t OK = kTRUE;
      new TGPrintDialog(fParent, this, 600, 150, 
			&fPrinter, &fPrinterCommand, &OK); 
      if ( ! OK ) return;
    }

  string filename = saveAll();
  if ( filename == "" ) return;

  // Send print command to shell

  string command(string(fPrinterCommand)+" "+
                 string(fPrinter)+" "+
                 filename+" &");
  write("%s",command.c_str());
  system(command.c_str());

  if ( fDebug )
    Info("printAll","Exit");
}
//---------------------------------------------------------------
// Save all pages to a single postscript file
//---------------------------------------------------------------
string GMBrowser::saveAll()
{
  if ( fDebug )
    Info("saveAll","Entered");

  if ( fPageList.size() == 0 ) return string("");

  char filename[160];
    
  // GAF - 24 Aug 2010 - separate run and subrun number 
  /*sprintf(filename,"%s/GMpages_run%6.6d.ps",fFileSaveInfo.fIniDir, 
          fRunNumber);*/
  sprintf(filename,"%s/GMpages_run%4.4d_%4.4d.ps",fFileSaveInfo.fIniDir, 
          fRun,fSubRun);

  if ( fDebug )
    Info("saveAll","file<%s>", filename);
		  
  // ECC - 26 Mar 2003 - create new window to write stuff in.
  //                   - window is 400x200 and located at (200,200).
  TCanvas *cmsg;
  cmsg = new TCanvas("cmsg", "GMbrowser Message", 200, 200, 400, 200);
  cmsg->cd();
  // Write message in canvas to indicate that we're updating.
  TLatex l;
  l.SetTextAlign(22);  // centered
  l.SetTextColor(9);   // blue
  l.SetTextSize(0.2);
  l.DrawLatex(0.5, 0.5, "Writing PostScript file...");// At center of canvas.
  cmsg->Update();


  // ECC - First open .ps file.
  TCanvas *curr_canvas = Canvas();
  string fname(filename); fname += "(";   // used to denote start of file.
  string fname2(filename); fname2 += ")"; // used to denote end of file.

  // ECC - 26 Mar 2003 - modify message on page when updating.
  int last = (int)fPageList.size();

  // Loop over all pages.
  for (int page = 0; page < last; page++) {
    DrawPage(page);

    // write to .ps file
    if (page == 0)           curr_canvas->Print(fname.c_str());
    else if (page == last-1) curr_canvas->Print(fname2.c_str());
    else                     curr_canvas->Print(filename);
  }

  // Draw the current page.
  DrawPage();

  if ( fDebug )
    Info("saveAll","Exit");

  delete cmsg;

  return string(filename);
}


//---------------------------------------------------------------
// Save all pages to a .gif files.
//---------------------------------------------------------------
void GMBrowser::saveGIF()
{
  if ( fDebug )
    Info("saveGIF","Entered");

  if ( fPageList.size() == 0 ) return;

  // ECC - 26 Mar 2003 - create new window to write stuff in.
  //                   - window is 400x200 and located at (200,200).
  TCanvas *cmsg;
  cmsg = new TCanvas("cmsg", "GMbrowser Message", 200, 200, 400, 200);
  cmsg->cd();
  // Write message in canvas to indicate that we're updating.
  TLatex l;
  l.SetTextAlign(22);  // centered
  l.SetTextColor(9);   // blue
  l.SetTextSize(0.2);
  l.DrawLatex(0.5, 0.5, "Writing GIF files...");// At center of canvas.
  cmsg->Update();


  TCanvas *curr_canvas = Canvas();

  // ECC - 26 Mar 2003 - modify message on page when updating.
  int last = (int)fPageList.size();

  // Loop over all pages.
  for (int page = 0; page < last; page++) {

    // Generate name for this plot.
    string name = makeName(fFileSaveInfo.fIniDir, page);

    // Names for .ps and .gif files.
    string psfile(name.c_str());  psfile += ".eps";
    string giffile(name.c_str()); giffile += ".gif";

    DrawPage(page);

    // Generate .ps file.
    curr_canvas->SaveAs(psfile.c_str());

    // Now convert it to a .gif file. Use system calls to do this conversion...
    string command("(pstopnm -ppm  -xsize 900 -ysize 925 -portrait -stdout ");
    command += psfile.c_str(); command +=  " | ppmtogif  > "; 
    command += giffile.c_str(); command += ") >& /dev/null ";
    gSystem->Exec(command.c_str());

    // Remove .eps file.
    string command2("rm "); command2 += psfile.c_str();
    gSystem->Exec(command2.c_str());
  }

  // Draw the current page.
  DrawPage();

  if ( fDebug )
    Info("saveGIF","Exit");

  delete cmsg;

  return;
}

// Generates plot title prefix from configuration file name.
void GMBrowser::plotTitle(const char *cFile)
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

string GMBrowser::makeName(const char *dir, int page)
{
  if ( page < 0 ) page = fCurrentPage;

  TString name(Page(page)->PageName.c_str());
  name = name.ReplaceAll(" ","_");
  name = name.ReplaceAll("/","_");  
  name = name.ReplaceAll(":","_");  
  name = name.ReplaceAll("(","_");  
  name = name.ReplaceAll(")","_");  

  char file[256];
  sprintf(file,"%s/%s%2.2d_%s", dir, fTitle.c_str(), page, name.Data());
  return string(&file[0]);
}

//---------------------------------------------------------------
// Adds names of histogram files so that multiple files can be
// summed together.
//---------------------------------------------------------------
void GMBrowser::addfile(string &name, int flag) {
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
TH1* GMBrowser::SumHists( int sumfile, TH1* hin, const char *name) {

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
void GMBrowser::load(string &name, int file, int *index)
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
GMBrowser::~GMBrowser()
{
  fPicturePool->FreePicture(fOpenIcon);
  fPicturePool->FreePicture(fSaveIcon);
  fPicturePool->FreePicture(fPrintIcon);

  delete fButtonLayout;
  delete fCanvasLayout;
  delete fClock;
  delete fCycleButton;
  delete fCycleClock;
  delete fGMEnv;
  delete fMenuBar;
  delete fMenuBarHelpLayout;
  delete fMenuBarItemLayout;
  delete fMenuBarLayout;
  delete fMenuFile;
  delete fMenuHelp;
  delete fOpenButton;
  delete fPicturePool;
  delete fPrintButton;
  delete fProgress;
  delete fProgressLayout;
  delete fListBox;
  delete fSaveButton;
  delete fStatusBar;
  delete fStatusBarLayout;
  delete fToolBar;
  delete fToolBarLayout;
  delete fUpdateButton;
  delete fUpdateClock;
  delete fVFrame;
  delete fVFrameLayout;
  delete fListBoxLayout;
  delete fPaveLabel;
  delete fPad;
  if (fWidgets) fWidgets->Delete();
  delete fWidgets;
}

// This method allows one to load a new configuration file when
// a histogram is clicked on.
void GMBrowser::catch_event(Int_t event, Int_t x, Int_t y, TObject *selected)
{
  TCanvas *c = (TCanvas *) gTQSender;

  // Select events with left mouse click.
  if (event != 61) return;

  // Determine which pad has been selected.
  TPad *pad = (TPad *)c->GetSelectedPad();

  // The pad name has been renamed to the name of the configuration file.
  TString config_file(pad->GetTitle());

  // Now load a new configuration file.
  if (config_file.Contains(".cfg")) {

    // Store name of current config file so that we can return to it.
    ConfigList.push_back(configFile_curr);
    PrevPageList.push_back(fCurrentPage);

    // Write out files in ConfigList
    //      for (int ifile=0; ifile<ConfigList.size(); ifile++) {
    //        cout << "ConfigList #" << ifile <<" :" << 
    //                ConfigList[ifile].c_str() << endl;
    //      }

    // Load new configuration file.
    GMBrowser::configureFile(config_file.Data());
  }
}
void GMBrowser::set_color_palette(){
  //const Int_t NRGBs = 5;
  const Int_t NRGBs = 6;
  const Int_t NCont = 150;

  //-- rainbow scale
  //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  //Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  //Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  //Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  
  //-- gray scale
  //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  //Double_t red[NRGBs]   = { 1.00, 0.84, 0.61, 0.34, 0.00 };
  //Double_t green[NRGBs] = { 1.00, 0.84, 0.61, 0.34, 0.00 };
  //Double_t blue[NRGBs]  = { 1.00, 0.84, 0.61, 0.34, 0.00 };
  
  //-- white to red scale
  //Double_t stops[NRGBs] = { 0.00, 0.12, 0.41, 0.68, 1.00 };
  //Double_t red[NRGBs]   = { 0.95, 0.70, 0.87, 1.00, 1.00 };
  //Double_t green[NRGBs] = { 0.95, 0.70, 1.00, 0.50, 0.00 };
  //Double_t blue[NRGBs]  = { 0.95, 0.70, 0.12, 0.00, 0.00 };
  
  //-- gray to red to purple scale
  Double_t stops[NRGBs] = { 0.00, 0.03, 0.12, 0.35, 0.70, 1.00 };
  Double_t red[NRGBs]   = { 0.95, 0.70, 0.87, 0.90, 1.00, 0.20 };
  Double_t green[NRGBs] = { 0.95, 0.70, 0.80, 0.50, 0.00, 0.10 };
  Double_t blue[NRGBs]  = { 0.95, 0.70, 0.12, 0.00, 0.00, 0.50 };
  
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

}

