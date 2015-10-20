#ifndef GMBROWSER_HPP
#define GMBROWSER_HPP
#define fLineSize 255
#define K_MAX_FILE 3
//////////////////////////////////////////////////////////////////////////////
// File:         GM.hpp
// Description:  Histogram browser for Global Monitor
// Created:      August 2002 Harrison B. Prosper
//               09-Oct-2002 V3.0
//               04-Dec-2002 HBP, Add TPostScript
//               15-Dec-2002 HBP, Add fPaveLabel and fPad
//               20-Dec-2002 HBP, Add AutoSave
//               24-Apr-2003 ECC, Define fLineSize
//////////////////////////////////////////////////////////////////////////////
//$Revision: 1.2 $

#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <list>
#include <vector>
#include <math.h>

#include <TObject.h>
#include <TEnv.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TError.h>
#include <TStyle.h>
#include <TPad.h>
#include <TFrame.h>
#include <TGFrame.h>
#include <TGWindow.h>
#include <TGPicture.h>
#include <TGButton.h>
#include <TGLayout.h>
#include <TGLabel.h>
#include <TGMenu.h>
#include <TGFileDialog.h>
#include <TGTextEditDialogs.h>
#include <TGTextEntry.h>
#include <TGTextBuffer.h>
#include <TGNumberEntry.h>
#include <TGStatusBar.h>
#include <TGProgressBar.h>
#include <TGTab.h>
#include <TGListBox.h>
#include <TG3DLine.h>
#include <TGToolBar.h>
#include <TGListView.h>
#include <TGCanvas.h>
#include <TGSplitter.h>
#include <TCanvas.h>
#include <TRootEmbeddedCanvas.h>
#include <TRootHelpDialog.h>
#include <TList.h>
#include <TString.h>
#include <TMap.h>
#include <TKey.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TTimer.h>
#include <TUnixSystem.h>
#include <TGFSContainer.h>
#include <TGMimeTypes.h>
#include <TPostScript.h>
#include <TPaveLabel.h>
#include <TPad.h>
#include <TImage.h>
#include <TExec.h>
#include <RQ_OBJECT.h>

class GMPage;
class GMTimer;
class GMEnv;

typedef std::vector<std::string>           vString;
typedef std::list<GMPage>                  GMPageList;
typedef std::map<std::string, std::string> HistMap;
typedef std::vector<int>                   vInt;


/** @name GMBrowser
    @version v1.0, August-2002
    @author Harrison B. Prosper
 */


//@{


/** Global Monitor Histogram Browser (GMBrowser). 
    This is a simple graphical user interface consisting of<br> 
    1. A menu bar<br> 
    2. A tool bar<br>
    3. A horizontal frame containing a listbox and a canvas<br>
    4. A status bar<br>
    <p>
    Each <b>GMBrowser</b> page pertains one or more (optionally) updating 
histograms. Only one 
page can be viewed at a time, however, one can cycle through all pages 
automatically.
A page is made visible by clicking on the name of the page in the
listbox.
<p>
The GMBrowser pages are configured using information provided in a text file
called a <b>configuration file</b> with file extension <i>.cfg</i>, while
general preferences are specified in the hidden file $HOME/.GMBrowser.
 */
class GMBrowser : public TGMainFrame
{

  RQ_OBJECT("GMBrowser")

public:

  /** Constructor.
      @param name   Name of application 
      @param width  Width of browser in pixels
      @param height Height of browser in pixels
  */

  GMBrowser(const char *name, UInt_t width=800, UInt_t height=600,
	    Bool_t debug=kFALSE);  
	    	  	
  ///
  ~GMBrowser();

  // Methods

  /**
   */
  const char       *Version();

  /** Configure browser.
      @param configFile Name of configuration file (see, for example,
      $GM_CONFIG_DIR/GM.cfg)
  */
  void              Configure(const char *configFile);

  /** Add a page to browser.
      @param name   Page name
      @return       Page number (counting from zero)
   */
  int               AddPage(const char *name, int divx=1, int divy=1);

  /** Draw the specified page.
      @param page   Page number
  */
  void              DrawPage(int page=-1);

  /** Delete a page from the browser.
      @param name   Page number
   */
  void              DeletePage(int page=-1);

  /** Delete all pagees from the browser.
   */
  void              DeleteAll();

  /** Set page style
   */
  void PageStyle(TCanvas *, GMPage *);

  /** Add help file as a single (large) string.
      @param text   Help text
   */
  void              AddHelp(const char *text);

  /// Let's rock!
  void              Run();

  ///
  GMPage           *Page(int page=-1);

  ///
  GMPageList       *PageList();

 
  /** Return pointer to the canvas.
   */
  TCanvas          *Canvas();

  /** Return pointer to the specified histogram. The name is the key by
      which the histogram is fetched from the root-file.
   */
  TObject          *Histogram(const char *name, int which=0);

  
  ///
  void              DebugOn(){fDebug = kTRUE;}
  ///
  void              DebugOff(){fDebug = kFALSE;}

  Bool_t            ProcessMessage(Long_t msg, Long_t id, Long_t parm);

  /// Close window and exit browser.
  void              CloseWindow();

  void catch_event(Int_t event, Int_t x, Int_t y, TObject *selected);

  void set_color_palette();

private:

  void   configure();
  void   configureFile(const char * config_file);
  void   print(Bool_t callPrintDialog);
  void   printAll(Bool_t callPrintDialog);
  void   write(const char *fmt, const char *str);
  void   saveAs(const char *type=0, const char *name=0);
  std::string saveAll();
  void   saveGIF();
  // Determines options for plotting histograms.
  void   histOptions(std::string all_options,  
                     std::string &config_file, TH1* h);
  int    optionValue(std::string option, std::string search, int nvalue,
                     std::string &value1, std::string &value2);
  long int    getRunNumber(const char *file, TFile *rootfile);
  int    getTotalEvents(TFile *rootfile);
  Bool_t checkKeys(int K);
  void   handleMenu(Int_t id);
  void   handleButton(Int_t id);
  void   handleListBox(Int_t page);
  void   handlePage(Int_t id);
  void   load(std::string &name, int file, int *index);
  void   addfile(std::string &name, int flag);
  TH1*   SumHists(int sumindex, TH1* h, const char *name);
  void   setTime();                         // Clock
  void   update();                          // Update current page
  void   cycle();                           // Cycle through pages
  void   autoSave();                        // Save all pages
  void   saveWWW(int page=-1);
  void   savePS(int page=-1);
  void   plotTitle(const char *cFile);
  std::string makeName(const char *dir, int page);
  Bool_t openFiles(int sumfile);
  void   closeFiles();

  GMTimer                *fClock;           // Handles setTime()
  GMTimer                *fUpdateClock;     // Handles update()
  GMTimer                *fCycleClock;      // Handles cycle()
  GMTimer                *fAutoSaveClock;   // Handles autosave()
  TStyle                 *GMStyle;          // Style for gmbrowser.


  const TGWindow         *fParent;          // Pointer to window manager
  Bool_t                  fDebug;
  Bool_t                  fUpdate;          // True if updating
  Bool_t                  fCycle;           // True if cycling
  Bool_t                  fCyclePlots;      // True if automatically saving
                                            //   plots when cycling.
  Bool_t                  fCyclePSPlots;    // True if automatically saving
                                            //   .ps plots when cycling.
  Bool_t                  fCallPrintDialog; 
  int                     fLastPage;        // Last page displayed
  int                     fUpdatePeriod;    // Measured in seconds
  int                     fCyclePeriod;     // Measured in seconds
  int                     fCurrentPage;
  std::string             fHelpText;
  long int                     fRunNumber;
  int                     fRun;
  int                     fSubRun;
  std::string             fTotalEventHist, fRunNumberHist;
  Bool_t                  fAutoSave;
  int                     fAutoSavePeriod;
  std::string             fTitle;

  // Histogram options:
  float option_y_min, option_y_max;  // min and max of histogram
  float option_range_min, option_range_max;  // X axis limits.
  int   option_rebin;                 // Set rebinning factor.
  float option_scale_ymin;           // ymin set to fraction of ymax.
  int   option_stat_value;           // statistics box bit map
  int   option_lwid;                 // line width of histogram
  int   option_scaleRef;             // scale reference histogram to data plot
  int   option_goodness;             // turn on/off Kol fit
  int   option_selected;             // root options selected
  int   option_PlotRef;              // turn on/off reference histogram
  int   option_PlotRef2;             // turn on/off reference 2 histogram
  int   option_logy;                 // Set when log of y axis is on.

  // Handle info to/from $HOME/.GMBrowser
  // Perculiar problem: The char * variables must be defined
  // here before other variables, otherwise one gets a seg-fault!
  // Why, I don't know
  GMEnv                  *fGMEnv;
  char                   *fPrinterCommand;
  char                   *fPrinter;
  char                   *fWWWDir;          // Web location for GIF files
  int                     fCanvasColor;

  TGFileInfo              fFileLoadInfo;    // Directory containing configs.
  TGFileInfo              fFileSaveInfo;    // Directory to save files
  GMPageList              fPageList;        

  vString                 ConfigList;
  vInt                    PrevPageList;

  TGVerticalFrame        *fVFrame;          // Frame for toolbar, fHFrame, etc.
  TGHorizontalFrame      *fHFrame;          // For listbox and canvas
  TGVerticalFrame        *fV1, *fV2;         // Individual frame for listbox
                                            // and canvas
  TGMenuBar              *fMenuBar;         // Menu objects
  TGPopupMenu            *fMenuFile;
  TGPopupMenu            *fMenuEdit;
  TGPopupMenu            *fMenuHelp;

  TGHorizontalFrame      *fToolBar;         // and toolbar

  TGPictureButton        *fOpenButton;
  TGPictureButton        *fSaveButton;
  TGPictureButton        *fPrintButton;
  TGTextButton           *fUpdateButton;
  TGTextButton           *fCycleButton;
  TGTextButton           *fConfigButton;
  TGTextButton           *fPlotButton;

  TGPicturePool          *fPicturePool;     // Manages icons
  const TGPicture        *fOpenIcon;
  const TGPicture        *fSaveIcon;
  const TGPicture        *fPrintIcon;

  TGHProgressBar         *fProgress;
  int                     fProgressCount;   // 0 <= fProcessCount < 100

  TGListBox              *fListBox;
  TRootEmbeddedCanvas    *fECanvas;
  TGStatusBar            *fStatusBar;
  TPaveLabel             *fPaveLabel;
  TPad                   *fPad;

  // E. Cheu - There are only three root files allowed.
  TFile                  *fRootFile[K_MAX_FILE];        // Root files

  // E. Cheu define a vector of TFiles for summing histograms.
  std::vector< std::vector<TFile *> >     fsumFiles;
  std::vector< std::vector<std::string> > fsumName;
  // Number of different pages that sum histograms.
  int total_sumFiles;

  // List of files used by gmbrowser.
  std::vector<std::string>    fRootFilename;

  // Layouts
 
  TGLayoutHints          *fToolBarLayout;
  TGLayoutHints          *fButtonLayout;
  TGLayoutHints          *fMenuBarLayout;
  TGLayoutHints          *fMenuBarItemLayout;
  TGLayoutHints          *fMenuBarHelpLayout;
  TGLayoutHints          *fListViewLayout;
  TGLayoutHints          *fCanvasLayout;
  TGLayoutHints          *fStatusBarLayout;
  TGLayoutHints          *fHFrameLayout;
  TGLayoutHints          *fVFrameLayout;
  TGLayoutHints          *fV1Layout, *fV2Layout, *splitLayout;
  TGLayoutHints          *fProgressLayout;
  TGLayoutHints          *fListBoxLayout;

  TList                  *fWidgets;
  // Well,...colors!

  ULong_t                 fBackgroundColor;
  ULong_t                 fRed;
  ULong_t                 fYellow;
  ULong_t                 fGreen;
  ULong_t                 fWhite;

};

//@}




//---------------------------------------------------------------
// This struct models a page
//---------------------------------------------------------------
struct GMPage
{
  std::string PageName; // Page Name
  int         DivX;     // Number of divisions in x
  int         DivY;     // Number of divisions in y
  vString     Name;     // Names of histograms
  vString     Title;    // Titles of histograms
  vString     Option;   // Draw options of histograms
  std::string MacroFile;    // Filename of Macro.
  std::string MacroFunc;  // Function in MacroFile.
  std::string HelpFile;  // Help file for this page.
  int         runMacro; // 1 if there is a macro.
  int         FileNumber; // File number associated with this page
  int         RefFileNum; // Reference File number
  int         RefFileNum2; // Reference 2 File number
  int         SumFileNum; // index number for summed histograms.
  int         Number;
  int         ShowRunNumber; // Displays run number for whole page.
};

//---------------------------------------------------------------
// Declare a pointer to a method of class GMBrowser.
// This complicated thing is used by the GMTimers. The latter manage
// synchronous operations such as updating.
//---------------------------------------------------------------
typedef void (GMBrowser::*GMmethod)();

//---------------------------------------------------------------
// We inherit from TTimer so that we can use the same semantics as
// it and avail ourselves of its HandleTimer method, which we overload. 
// GMTimer takes a GMBrowser object and a pointer to one of its private 
// methods, which is invoked each time GMTimer times out. 
// Note: Because GMTimer needs to access the private methods of
// GMBrowser we must regard the latter as a friend class of GMTimer.
//---------------------------------------------------------------
class GMTimer : public TTimer
{
public:
  GMTimer(GMBrowser *browser, GMmethod method) 
    : TTimer(1000),
      fBrowser(browser),
      fMethod(method)
  {SetObject(this); Stop();} 
  // Set object whose HandleTimer method is to be invoked
  // when GMTimer times out; here it is GTimer itself. That is, when GMTimer
  // times out it calls its own HandleTimer method, which in turn calls the
  // appropriate method of GMBrowser.

  ~GMTimer() {Stop();}


  // When GMTimer times out this method is called.
  // It calls the method of GMBrowser, specified in the constructor above.

  Bool_t HandleTimer(TTimer *timer) 
  { 
    (fBrowser->*fMethod)();
    return kTRUE;
  }

  friend class GMBrowser; // GMTimer needs access to guts of GMBrowser

private:
  GMBrowser *fBrowser;
  GMmethod   fMethod;
};

// Get/Set variables in $HOME/.GMBrowser

class GMEnv
{
public:
  GMEnv();
  ~GMEnv();

  const char *Get(const char *var, const char *defValue=0);
  void        Set(const char *var, const char *newValue);
};





#endif
