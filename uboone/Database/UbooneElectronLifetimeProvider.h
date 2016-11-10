/**
 * \file UbooneElectronLifetimeProvider.h
 *
 * \ingroup WebDBI
 * 
 * \brief Class def header for a class UbooneElectronLifetimeProvider
 *
 * @author eberly@slac.stanford.edu
 */

#ifndef UBOONEELECTRONLIFETIMEPROVIDER_H
#define UBOONEELECTRONLIFETIMEPROVIDER_H

// C/C++ standard libraries
#include <string>
#include <array>

// LArSoft libraries
#include "larevt/CalibrationDBI/IOVData/Snapshot.h"
#include "larevt/CalibrationDBI/IOVData/IOVDataConstants.h"
#include "larevt/CalibrationDBI/Providers/DatabaseRetrievalAlg.h"
#include "larevt/CalibrationDBI/Interface/ElectronLifetimeProvider.h"

#include "larevt/CalibrationDBI/IOVData/ElectronLifetimeContainer.h"

namespace lariov {

  /**
   * @brief Retrieves channel information: pedestal and RMS
   * 
   * Configuration parameters
   * =========================
   * 
   * - *DatabaseRetrievalAlg* (parameter set, mandatory): configuration for the
   *   database; see lariov::DatabaseRetrievalAlg
   * - *UseDB* (boolean, default: false): retrieve information from the database
   * - *UseFile* (boolean, default: false): retrieve information from a file;
   *   not implemented yet
   * - *DefaultExp* (real, default: 400.0): exponential parameter returned
   *   when /UseDB/ and /UseFile/ parameters are false
   * - *DefaultIndMean* (real, default: 2048.0): constant parameter returned
   *   when /UseDB/ and /UseFile/ parameters are false
   */
  class UbooneElectronLifetimeProvider : public DatabaseRetrievalAlg, public ElectronLifetimeProvider {
  
    public:
    
      /// Constructors
      UbooneElectronLifetimeProvider(const std::string& foldername, 
      			      const std::string& url, 
			      const std::string& tag="");
	
      UbooneElectronLifetimeProvider(fhicl::ParameterSet const& p);
      
      /// Reconfigure function called by fhicl constructor
      void Reconfigure(fhicl::ParameterSet const& p);
      
      /// Update Snapshot and inherited DBFolder if using database.  Return true if updated
      bool Update(DBTimeStamp_t ts);
      
      /// Retrieve lifetime information
      float Lifetime(float t) const override;
      float Purity() const override;
      float LifetimeErr(float t) const override;
      float PurityErr() const override;
      
      const ElectronLifetimeContainer& LifetimeContainer() const;
      float ExpOffset() const;
      float TimeConstant() const;
      float ExpOffsetErr() const;
      float TimeConstantErr() const;
           
      //hardcoded information about database folder - useful for debugging cross checks
      static constexpr unsigned int NCOLUMNS = 5;
      static constexpr const char* FIELD_NAMES[NCOLUMNS]
        = {"channel", "exponential_offset", "err_exponential_offset", "time_constant", "err_time_constant"};
      static constexpr const char* FIELD_TYPES[NCOLUMNS]
        = {"bigint", "real", "real", "real", "real"};
      
    private:
    
      DataSource::ds fDataSource;
          
      Snapshot<ElectronLifetimeContainer> fData;
      
      const unsigned int fLifetimeChannel = 0;
  };
}//end namespace lariov

#endif
