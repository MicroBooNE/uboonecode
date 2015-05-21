#ifndef __UBOpChannelTypes__

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
				       (Undefined)\
				       (FEMCosmicHighGain)\
				       (FEMCosmicLowGain)\
				       (FEMBeamHighGain)\
				       (FEMBeamLowGain)\
				       (FEMCosmicLogicPulse)\
				       (FEMBeamLogicPulse)\
				       (FEMBeamTriggerBNB)\
				       (FEMBeamTriggerNUMI)\
				       (FEMFlasherLogicPulse) )

  DEFINE_ENUM_WITH_STRING_CONVERSIONS( UBOpticalChannelGain_t, \
				       (HighGain)\
				       (LowGain)\
				       (LogicChannel) )
  
}

#endif
