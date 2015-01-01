/*!
 * Title:   ROIAlg base class
 * Author:  wketchum@lanl.gov
 * Inputs:  
 * Outputs: 
 *
 * Description:
 * This is a base class for ROIAlgs.
 */

#include "ROIAlg.h"
#include "ROIAlg_DigitAboveThreshold.h"

std::unique_ptr<ROIAlg> util::ROIAlg::MakeROIAlg(fhicl::Parameterset const& p){

  fAlgName = p.get<std::string>("AlgName");

  std::unique_ptr<ROIAlg> ptr;
  if(fAlgName.compare("DigitAboveThreshold"))
    ptr.swap(new ROIAlg_DigitAboveThreshold(p));
  else
    throw "ERROR in ROIAlg: No registered ROIAlg with that name.";

  return std::move(ptr);
}

void util::ROIAlg::InsertRegion(Region region, Region_t type, Tick hint){

  auto const& it = fRegions.emplace_hint( hint, std::make_pair(region,type));

  auto const& next_it = std::next(it);
  auto const& prev_it = std::prev(it);

  if(CheckRegionOverlap(*it,*next_it) || CheckRegionOverlap(*it,*prev_it))
    throw "ERROR in ROIAlg: Insterted an overlapping region into Region Set."

}
