/**
 * @file   IFDHFileTransfer.h
 * @brief Utility to copy files with IFDH
 * Author: Gianluca Petrillo
 */

#ifndef UBOONE_EVENTWEIGHT_IFDHFILETRANSFER_H
#define UBOONE_EVENTWEIGHT_IFDHFILETRANSFER_H

// support libraries
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "ifdh.h" 

// C/C++ standard libraries
#include <memory> // std::make_unique()
#include <string>

namespace evwgh {
  
  /**
   * @brief Manager for copy of files by IFDH
   *
   * This class manages file transfers via IFDH.
   * It cleans up the temporary files on destruction.
   * Example of usage:
   *
   * 1. have an instance of IFDHFileTransfer in the class that needs
   *    the external files:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * class MyAlgo {
   *   IFDHFileTransfer IFDH;
   *
   *   // ...
   *
   * }; // class MyAlgo
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * 2. fetch the files
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * void MyAlgo::DoSomething() {
   *   
   *   // get the files we need!
   *   std::string file1 = IFDH.fetch("/pnfs/path/file1.root");
   *   std::string file2 = IFDH.fetch("/pnfs/path/file2.root");
   *   
   *   // do something!
   *   TFile F1(file1.c_str());
   *   if (!F1.IsOpen()) {
   *     // error!!
   *   }
   *   // do the rest
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * 3. clean up will automatically happen when `MyAlgo` is destroyed.
   *
   */
  class IFDHFileTransfer {
    
    /// local IFDH client
    std::unique_ptr<ifdh_ns::ifdh> fIFDHclient = nullptr;
    
    void createClient()
      {
        if (!fIFDHclient) {
          mf::LogDebug("IFDHFileTransfer")
            << "[" << ((void*) this) << "] Initializing a IFDH client";
          fIFDHclient = std::make_unique<ifdh_ns::ifdh>();
        }
      }
    
      public:
    /// Destructor: clean up
    ~IFDHFileTransfer()
      {
        if (fIFDHclient) {
          mf::LogDebug("IFDHFileTransfer")
            << "[" << ((void*) this) << "] Cleaning up";
          fIFDHclient->cleanup();
        }
      }

    /**
     * @brief Copies a file from IFDH
     * @param file full path of the file to be copied
     * @return path to the local copy of the file
     *
     * The local copy will be cleaned when this object is destructed. 
     */  
    std::string fetch(std::string file)
      {
        createClient();
        std::string localPath = fIFDHclient->fetchInput(file);//****fetchInput instead of fetch?
        if (localPath.empty()) {
          mf::LogError("IFDHFileTransfer")
            << "[" << ((void*) this) << "] failed to fetch '" << file << "'";
          throw art::Exception(art::errors::NotFound)
            << "IFDH failed to fetch '" << file << "'!\n";
        }
        mf::LogInfo("IFDHFileTransfer")
          << "[" << ((void*) this) << "] fetched '"
          << file << "' as '" << localPath << "'";
        return localPath;
      } // fetch()
    
    
  }; // class IFDHFileTransfer

} // namespace evwgh 


#endif //  UBOONE_EVENTWEIGHT_IFDHFILETRANSFER_H

