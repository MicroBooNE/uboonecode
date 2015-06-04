#ifndef __UBOpChannelTypes__
#define __UBOpChannelTypes__

#include <boost/preprocessor.hpp>

#define X_DEFINE_ENUM_WITH_STRING_CONVERSIONS_TOSTRING_CASE(r, data, elem)    \
  case elem : return BOOST_PP_STRINGIZE(elem);

#define DEFINE_ENUM_WITH_STRING_CONVERSIONS(name, enumerators)   \
  enum name {                                                    \
    BOOST_PP_SEQ_ENUM(enumerators)				 \
  };                                                             \
                                                                 \
  inline const char* UBOpChannelEnumName(name v)                 \
  {                                                              \
  switch (v)                                                     \
    {                                                            \
      BOOST_PP_SEQ_FOR_EACH(					 \
	X_DEFINE_ENUM_WITH_STRING_CONVERSIONS_TOSTRING_CASE,     \
	name,							 \
	enumerators )						 \
								 \
    default: return "[Unknown " BOOST_PP_STRINGIZE(name) "]";	 \
    }								 \
  }

namespace opdet {


  DEFINE_ENUM_WITH_STRING_CONVERSIONS( UBOpticalChannelCategory_t, \
				       (Uncategorized)\
				       (UnspecifiedLogic)\
				       (OpdetCosmicHighGain)\
				       (OpdetCosmicLowGain)\
				       (OpdetBeamHighGain)\
				       (OpdetBeamLowGain)\
				       (BNBLogicPulse) \
				       (NUMILogicPulse) \
				       (FlasherLogicPulse)\
				       (StrobeLogicPulse)\
				       (NumUBOpticalChannelCategories) )

  DEFINE_ENUM_WITH_STRING_CONVERSIONS( UBOpticalChannelGain_t, \
				       (Undefined)\
				       (HighGain)\
				       (LowGain)\
				       (LogicChannel)\
				       (NumUBOpticalChannelGains) )
  
    UBOpticalChannelGain_t GetUBTypeFromCategory( UBOpticalChannelCategory_t cat ) {

    UBOpticalChannelGain_t chtype = Undefined;
    switch ( cat ) {
    case UnspecifiedLogic:
    case BNBLogicPulse:
    case NUMILogicPulse:
    case FlasherLogicPulse:
    case StrobeLogicPulse:
      chtype = LogicChannel;
      break;
    case OpdetCosmicHighGain:
    case OpdetBeamHighGain:
      chtype = HighGain;
      break;
    case OpdetCosmicLowGain:
    case OpdetBeamLowGain:
      chtype = LowGain;
      break;
    default:
      chtype = Undefined;
      break;
    }
    return chtype;
  }
}

#endif
