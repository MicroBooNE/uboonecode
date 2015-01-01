/*!
 * Title:   Range class, for use with ROIAlg
 * Author:  kazuhiro@nevis.columbia.edu, wketchum@lanl.gov
 * Inputs:  
 * Outputs: 
 *
 * Description:
 * This is a simple class for using an "range". It's a std::pair between two objects,
 * like iterators or even "tick" values. Utilities here for merging as well.
 */

#include "Range.h"

template<class T>
size_t util::UniqueRangeSet::Insert(Range<T> const& a){

  if(a.RangeType()==RangeType_t::kUndefined)
    throw std::runtime_error("Cannot use undefined range types in UniqueRangeSet!");
  
  const size_t current_range_size = range_set.size();
  
  //simple case: the ranges don't overlap at all
  auto res = range_set.insert(a);
  if(res.second) return 1;
  
  //OK, now we need to handle an "equal" range case
  auto& iter = res.first;
  InsertOverlappingRange(a,*iter);

  return (range_set.size() - current_range_size);
}


template<class T>
void util::UniqueRangeSet::InsertOverlappingRanges(Range<T> const& a,Range<T> const& b){

  if(a.IsSubRangeOf(b) || b.IsSubRangeOf(a))
    InsertSubRange(a,b);
  else
    InsertPartiallyOverlappingRanges(a,b);   

}

template<class T>
void util::UniqueRangeSet::InsertSubRange(Range<T> const& a,Range<T> const& b){
  
  Range<T> const& sub_range = (a.IsSubRangeOf(b)) a ? b;
  Range<T> const& super_range = (a.IsSubRangeOf(b)) b ? a;

  //if the existing "super-range" is already signal type, nothing to do
  if(super_range.RangeType() == RangeType_t::kSignal)
    return;
  
  //if super-range is baseline, but sub-range is too, then nothing to do
  if(sub_range.RangeType() == RangeType_t::kBaseline)
    return;
  
  //so sub_range is Signal, and super_range is baseline!
  //We should erase super_range from range set, and
  //add in three new ranges: b1, sub_range, b2
  range_set.erase(super_range);
  range_set.emplace(super_range.Start(),sub_range.Start(),RangeType_t::kBaseline);
  range_set.insert(sub_range);
  range_set.emplace(sub_range.End(),super_range.End(),RangeType_t::kBaseline);
}

templave<class T>
void util::UniqueRangeSet::InsertPartiallyOverlappingRanges(Range<T> const& a, Range<T> const& b){

  if(a.RangeType()==b.RangeType())
    InsertOverlappingRangesSameTypes(a,b);
  else
    InsertOverlappingRangesDifferentTypes(a,b);

}

template<class T>
void util::UniqueRangeSet::InsertOverlappingRangesSameTypes(Range<T> const& a,Range<T> const& b){

  if(b.RangeType()!=a.RangeType())
    throw std::runtime_error("Ranges need to be same types in this insert call!");

  Range<T> tmp_a = a;
  tmp_a.Merge(b);
  range_set.insert(tmp_a);
}

  template<class T>
void util::UniqueRangeSet::InsertOverlappingRangesDifferentTypes(Range<T> const& a,Range<T> const& b){

  if(b.RangeType()==a.RangeType())
    throw std::runtime_error("Ranges need to be different types in this insert call!");
  
  Range<T> const& signal_range = (a.RangeType()==RangeType_t::kSignal) a ? b;
  Range<T> const& baseline_range = (a.RangeType()==RangeType_t::kBaseline) a ? b;

  if(signal_range.Start() < baseline_range.Start() && signal_range.End() < baseline_range.End() ){
    range_set.emplace(signal_range.Start(),signal_range.End(),RangeType_t::kSignal);
    range_set.emplace(signal_range.End(),baseline_range.End(),RangeType_t::kBaseline);
  }
  else if(signal_range.Start() > baseline_range.Start() && signal_range.End() > baseline_range.End() ){
    range_set.emplace(baseline_range.Start(),signal_range.Start(),RangeType_t::kBaseline);
    range_set.emplace(signal_range.Start(),signal_range.End(),RangeType_t::kSignal);
  }
  else
    throw std::runtime_error("Ranges should overlap but they don't. This is messed up.");
  
}
