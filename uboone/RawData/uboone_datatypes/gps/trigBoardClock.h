#ifndef _UBOONE_TRIGBOARDCLOCK_UB_H
#define _UBOONE_TRIGBOARDCLOCK_UB_H  1


// This struct will be a key in a map, so I must define "<".
typedef struct tbclkub
{
  uint32_t frame;
  uint16_t sample;
  uint16_t div;

    bool operator<(const tbclkub& mk) const 
    {
        if (frame < mk.frame)       
        {            
	  return true; // meaning: keep sorting.
	}
	return false;
    }

  tbclkub (uint32_t f=0, uint16_t s=0, uint16_t d=0): frame(f),sample(s),div(d) {  }

} tbclkub_t;

#endif // #define _UBOONE_TRIGBOARDCLOCK_UB_H  
