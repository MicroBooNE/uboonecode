#ifndef RANGE_H
#define RANGE_H

/*!
 * Title:   Range class and UniqueRangeSet, for use with ROIAlg
 * Author:  kazuhiro@nevis.columbia.edu, wketchum@lanl.gov
 * Inputs:  
 * Outputs: 
 *
 * Description:
 * This is a simple class for using an "range". It's a std::pair between two objects,
 * like iterators or even "tick" values. Utilities here for merging as well.
 */

#include <vector>
#include <string>
#include <exception>

#include "fhiclcpp/ParameterSet.h"

namespace util{

  typedef enum{
    kUndefined,
    kBaseline,
    kSignal
  } RangeType_t;
  
  template <class T>
    class Range : public std::pair<T,T> {
    
  public:
    Range(const T& start,
	  const T& end,
	  const RangeType_t& rt==RangeType_t::kUndefined)
      { 
	if(start>=end)
	  throw std::runtime_error("Inserted invalid range: end before start.");
	(*this).first = start;
	(*this).second = end;
	rangetype = rt;
      }
    
    /// Intuitive accessor
    const T& Start()               const { return (*this).first;  }
    const T& End()                 const { return (*this).second; }
    const RangeType_t& RangeType() const { return rangetype; }
    
    /// Ordering
    inline bool operator< (const Range& rhs) const
    {return ( (*this).second < rhs.first ); }

    ///Test if this range is subset of a given range
    bool IsSubRangeOf(const Range& a){
      return ( (*this.first) >= a.first && (this.second) <= a.second );
    }

    /// Merging utility
    void Merge(const Range& a) {
      (*this).first  = std::min( (*this).first,  a.first  );
      (*this).second = std::max( (*this).second, a.second );
    }
    
  private:
    RangeType_t rangetype;
    
  };
}

// Implement pointer comparison in case it's useful
template <class T>
class std::less<util::Range<T>*>
{
 public:
  bool operator()( const util::Range<T>* lhs, const util::Range<T>* rhs )
  { return (*lhs) < (*rhs); }
};

namespace util{

  //this class should provide a set of time-ordered, non-overlapping ranges.
  //Note that the constructer that takes a "start" and "end" point actually
  //creates a range for that whole region, set to be "Baseline". Inserting
  //signal ranges will carve out signal ranges from that larger baseline range.
  template <class T>
    class UniqueRangeSet {
  public:
    UniqueRangeSet(){}
    UniqueRangeSet(const Range<T>& a)
      { range_set.insert(a); }
    UniqueRangeSet(const T& start,const T& end)
      { range_set.emplace(start,end,RangeType_t::kBaseline); }
    
    /// Modified insert that merges overlapping range.
    //  Return number of added ranges
    size_t Insert(const Range<T>& a);
    
  private:
    std::set< Range<T> >  range_set;

    void InsertOverlappingRanges(Range<T> const& a,Range<T> const& b);
    void InsertSubRange(Range<T> const& a,Range<T> const& b);
    void InsertPartiallyOverlappingRanges(Range<T> const& a, Range<T> const& b);
    void InsertOverlappingRangesSameTypes(Range<T> const& a,Range<T> const& b);
    void InsertOverlappingRangesDifferentTypes(Range<T> const& a,Range<T> const& b);
    
  };

}

#endif
