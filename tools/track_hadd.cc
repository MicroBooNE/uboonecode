//
// Name: track_hadd
//
// Purpose: Stand alone executable for merging histogram files produced by 
//          TrackAna and SeedAna analyzer modules.
//
// Usage: This executable has an invocation interface similar to root's 
//        built-in hadd utility:
//
//        track_hadd [options] targetfile source1 [source2 source3 ...]
//
//        where target is the name of a target root file and sources
//        can be a single root file or a list file (if preceded by "@").
//
// Options:
//
// -f - Overwrite existing target file.
// -k - Skip unreadable source files.
//
// Created: H. Greenlee 19-May-2014

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "TFile.h"
#include "TH1.h"

// Function to print help message.

void help()
{
  std::cout << "\nUsage: track_hadd [-f] [-k] targetfile [@]source1 [[@]source2 [@]source3 ...]\n"
	    << "\n"
	    << "Options:\n"
	    << "\n"
	    << "-f - Overwrite existing target file.\n"
	    << "-k - Skip unreadable source files.\n" << std::endl;
}

// Function to test whether file exists.

bool file_exists(const std::string& fn)
{
  bool result = false;
  std::ifstream file(fn);
  if(file) {
    file.close();
    result = true;
  }
  return result;
}

// Fill 1-dimensional multiplicity histogram.

void mulcalc(const TH1* hnum, const TH1* hden, TH1* hmul)
{
  int nbins = hnum->GetNbinsX();
  if (nbins != hden->GetNbinsX())
    throw std::runtime_error("mulcalc: incompatible histograms (I)");
  if (nbins != hmul->GetNbinsX())
    throw std::runtime_error("mulcalc: incompatible histograms (II)");

  // Loop over bins, including underflow and overflow.

  for(int ibin = 0; ibin <= nbins+1; ++ibin) {
    double num = hnum->GetBinContent(ibin);
    double den = hden->GetBinContent(ibin);
    if(den == 0.) {
      hmul->SetBinContent(ibin, 0.);
      hmul->SetBinError(ibin, 0.);
    }
    else {
      double mul = num / den;
      if(mul < 0.)
	mul = 0.;
      double err = std::sqrt((1. + mul) * mul / den);
      hmul->SetBinContent(ibin, mul);
      hmul->SetBinError(ibin, err);
    }
  }
  hmul->SetMinimum(0.);
}

// Fill 1-dimensional efficiency histogram assuming binomial errors.

void effcalc(const TH1* hnum, const TH1* hden, TH1* heff)
{
  int nbins = hnum->GetNbinsX();
  if (nbins != hden->GetNbinsX())
    throw std::runtime_error("effcalc: incompatible histograms (I)");
  if (nbins != heff->GetNbinsX())
    throw std::runtime_error("effcalc: incompatible histograms (II)");

  // Loop over bins, including underflow and overflow.

  for(int ibin = 0; ibin <= nbins+1; ++ibin) {
    double num = hnum->GetBinContent(ibin);
    double den = hden->GetBinContent(ibin);
    if(den == 0.) {
      heff->SetBinContent(ibin, 0.);
      heff->SetBinError(ibin, 0.);
    }
    else {
      double eff = num / den;
      if(eff < 0.)
	eff = 0.;
      if(eff > 1.)
	eff = 1.;
      double err = std::sqrt(eff * (1.-eff) / den);
      heff->SetBinContent(ibin, eff);
      heff->SetBinError(ibin, err);
    }
  }
  heff->SetMinimum(0.);
  heff->SetMaximum(1.);
}

// Function to recalculate efficiency and multiplicity histograms.

void recalculate(TDirectory* dir)
{
  // Loop over contents of directory.

  TList* objs = dir->GetList();
  TIter next(objs);
  while(TObject* obj = next()) {

    // Subdirectories handled here.
    // Recursively recalculate subdirectories.

    if(obj->InheritsFrom(TDirectory::Class()))
      recalculate((TDirectory*)obj);

    // Histograms handled here.

    if(obj->InheritsFrom(TH1::Class())) {
      TH1* h = (TH1*)obj;

      // Efficiency histograms.
      // eX = gX / mcX;

      const char* hname = h->GetName();
      if(hname[0] == 'e') {

	// Find numerator and denominator histograms.

	std::string hnumname = std::string("g") + std::string(&hname[1]);
	TH1* hnum = 0;
	dir->GetObject(hnumname.c_str(), hnum);
	std::string hdenname = std::string("mc") + std::string(&hname[1]);
	TH1* hden = 0;
	dir->GetObject(hdenname.c_str(), hden);
	if(hnum != 0 && hden != 0) {
	  //std::cout << "Recalculating " << h->GetName() << std::endl;
	  effcalc(hnum, hden, h);
	}
      }

      // Multiplicity histograms.
      // mulX = mX / mcX;

      hname = h->GetName();
      if(std::string(hname).substr(0,3) == std::string("mul")) {

	// Find numerator and denominator histograms.

	std::string hnumname = std::string("m") + std::string(&hname[3]);
	TH1* hnum = 0;
	dir->GetObject(hnumname.c_str(), hnum);
	std::string hdenname = std::string("mc") + std::string(&hname[3]);
	TH1* hden = 0;
	dir->GetObject(hdenname.c_str(), hden);
	if(hnum != 0 && hden != 0) {
	  //std::cout << "Recalculating " << h->GetName() << std::endl;
	  mulcalc(hnum, hden, h);
	  h->SetMinimum(0.);
	}
      }
    }
  }
}

// Function to merge contents of source directory into target directory.
// This function recursively checks subdirectories and merges histograms
// into target directory.  Trees are not merged.

void merge_dir(TDirectory* targetdir, TDirectory* sourcedir)
{
  // Loop over contents of source directory.

  TList* objs = sourcedir->GetList();
  if(objs->GetSize() == 0) {
    sourcedir->ReadAll();
    objs = sourcedir->GetList();
  }
  TIter next(objs);
  while(TObject* obj = next()) {

    // Subdirectories handled here.

    if(obj->InheritsFrom(TDirectory::Class())) {

      // Check whether similarly named subdirectory exists in target
      // directory.  If not, create it.

      TDirectory* subdir = 0;
      targetdir->GetObject(obj->GetName(), subdir);
      if(subdir == 0) {
	subdir = targetdir->mkdir(obj->GetName(), obj->GetTitle());
	//std::cout << "Add target path " << subdir->GetPath() << std::endl;
      }

      // Merge subdirectories recursively.

      merge_dir(subdir, (TDirectory*)obj);
    }

    // Histograms handled here.

    if(obj->InheritsFrom(TH1::Class())) {
      TH1* sourceh = (TH1*)obj;

      // Check whether a similarly named histogram exists in target
      // directory.

      TH1* targeth = 0;
      targetdir->GetObject(sourceh->GetName(), targeth);
      if(targeth == 0) {

	// Histogram does not exist in target directory.
	// Clone source histogram.

	targetdir->cd();
	targeth = (TH1*)sourceh->Clone(sourceh->GetName());
	//std::cout << "Cloning histogram " << sourceh->GetName() << std::endl;
      }
      else {

	// Histogram already exists.
	// Add contents of source histogram.

	targeth->Add(sourceh);
      }
    }
  }
}

int main(int argc, char** argv)
{
  // Parse arguments.

  bool overwrite = false;
  bool skip = false;
  std::string target;
  std::vector<std::string> sources;

  for(int i=1; i<argc; ++i) {

    std::string arg(argv[i]);

    // Parse options.

    if(arg == "-f")
      overwrite = true;
    else if(arg == "-k")
      skip = true;

    // Process -v option (for compatibility with hadd, not implemented here).

    else if(arg == "-v")
      ++i;

    // Ignore any other arguments that begin with "-".

    else if(arg[0] == '-');

    // Extract target file name.

    else if(target.empty()) {
      target = arg;

      // Target file must not already exist unless overwrite option is specified.

      if(!overwrite && file_exists(target)) {
	//std::cout << "Target file " << target << " already exists." << std::endl;
	return 1;
      }
    }

    // Everything else is treated as a source file.
    // Can be plain root file or file list.

    else {
      if(arg[0] != '@') {

	// Not a list file (assume plain root file).
	// Check that source file exists.

	if(!skip && !file_exists(arg)) {
	  //std::cout << "Source file " << arg << " does not exist." << std::endl;
	  return 1;
	}
	sources.push_back(arg);
      }
      else {

	// List files handled here.

	std::string listfile = arg.substr(1);

	std::ifstream listf(listfile);
	if(!listf) {
	  std::cout << "Failed to open list file " << listfile << std::endl;
	  return 1;
	}
	std::string source;
	while(listf >> source) {
	  if(!skip && !file_exists(source)) {
	    //std::cout << "Source file " << source << " does not exist." << std::endl;
	    return 1;
	  }
	  sources.push_back(source);
	}
	listf.close();
      }
    }
  }

  // Done parsing arguments.
  // First do some sanity checks.

  if(target.empty() || sources.size() == 0) {
    help();
    return 1;
  }

  // Open target file for writing.

  TFile* targetf = new TFile(target.c_str(), "RECREATE");

  // Loop over source files.
  // Only open one source file at a time.

  for(auto const& source : sources) {

    //std::cout << "Merging " << source << std::endl;

    // Open source file for reading.
    // Use static method TFile::Open instead of TFile constructor, which 
    // method can handle urls (can return derived classes such as TNetFile,
    // TWebFile, etc.).

    TFile* sourcef = TFile::Open(source.c_str(), "READ");
    if(sourcef == 0 || !sourcef->IsOpen() || sourcef->IsZombie()) {
      if(skip) {
	std::cout << "Skipping source file " << source << std::endl;
	continue;
      }
      else {
	std::cout << "Failed to open source file " << source << std::endl;
	return 1;
      }
    }
    else {

      // File successfully opened.
      // Merge contents of this file.

      merge_dir(targetf, sourcef);

      // Close this soure.

      sourcef->Close();
      delete sourcef;
    }
  }

  // Recalculate efficiency and multiplicity histograms.
  // This is the reason we need a special program.

  recalculate(targetf);

  // Write everything out and close target file.

  if(targetf) {
    targetf->Write();
    targetf->Close();
  }

  // Done.
    
  return 0;
}
