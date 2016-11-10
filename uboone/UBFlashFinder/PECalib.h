/**
 * \file PECalib.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class PECalib
 *
 * @author drinkingkazu
 */

/** \addtogroup UBFlashFinder 

    @{*/
#ifndef PECALIB_H
#define PECALIB_H

#include "FlashFinderTypes.h"
#include "FlashFinderFMWKInterface.h"
#include <iostream>
#include <numeric>
#include <functional>
#include <algorithm>

namespace pmtana{
/**
   \class PECalib
   User defined class PECalib ... these comments are used to generate
   doxygen documentation!
 */

  class PECalib {
    
  public:
    
    /// Default constructor
    PECalib();
    
    /// Default destructor
    ~PECalib(){}

    void Configure(const Config_t &pset);

    double CosmicPE(const size_t opdet, const double area, const double amp) const;
    double BeamPE(const size_t opdet, const double area, const double amp) const;

  protected:

    std::vector<double> _spe_area_gain_v;
    std::vector<double> _spe_amp_gain_v;

    std::vector<double> _cosmic_ophit_correction_v;
    std::vector<double> _cosmic_ophit_ratio_mean_v;
    std::vector<double> _cosmic_ophit_ratio_std_v;
    std::vector<double> _cosmic_ophit_amppe_cut_v;
    std::vector<double> _cosmic_ophit_areape_cut_v;

    std::vector<double> _relative_qe_v;

  };
  
} 

#endif
/** @} */ // end of doxygen group 

