#ifndef LOGICPULSEFINDER_CXX
#define LOGICPULSEFINDER_CXX

#include "LogicPulseFinder.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

////////////////////////////////////////////////////////////////////////////////////
///Checks for TTL pulses and returns a vector of times at which they start.      ///
///Operating on the assumption that all TTL pulses are similar (they should be)  ///
///it uses the max value found as a standard by which to measure all potential   ///
///TTL pulses. It then returns the position associated with the 50% rising edge  ///
///of each TTL and outputs this list of rising edge positions to a vector 	 ///
////////////////////////////////////////////////////////////////////////////////////

template <class T>
std::vector<unsigned short> LogicPulseFinder<T>::Get_TTL_Starts(std::vector<T> wfm, T threshold )
{

	std::vector<unsigned short> TTL_times;

	auto max_amp = max_element(std::begin(wfm) , std:: end(wfm));
	auto minimum = min_element(std::begin(wfm) , std:: end(wfm));

	if ( (max_amp-minimum)<threshold )
	  return TTL_times; // no logic pulse exists!

	int in_pulse = 1;

	for (unsigned int i=0; i<wfm.size(); i++) {
		
	        if (float(wfm[i]) < *minimum+0.5*( *max_amp - *minimum)) {
			in_pulse = 0;
		}

		else if (in_pulse == 0) {
			
		        TTL_times.push_back(i);
		        in_pulse = 1;				
		}
	
	}

	return TTL_times;

}

////////////////////////////////////////////////////////////////////////////////////
/// Checks the TTL waveform to verify that it is a good one to use. It checks for///
/// 1) Presence of peak at the beginning of the waveform                         ///
/// 2) Makes sure there is no great variation in TTL amplitude                   ///
/// 3) MAYBE IMPLEMENT MORE CHECKS??						 ///
////////////////////////////////////////////////////////////////////////////////////

template <class T>
std::map<std::string,float> LogicPulseFinder<T>::TTL_Health(std::vector<T> wfm)
{

        std::map<std::string,float> TTL_Health_Flags;

////////////////////////////////////////////////////////////////////////////////////
/// The following segment does a simple average calculation of the first 10 bins ///
/// of the waveform and compares the average to the maximum amplitude found in   ///
/// waveform. If the average is above some threshold (chosen somewhat arbitrarily///
/// chosen to be 25% of max) then the waveform is flagged as likely having a peak///
/// at the beginning of the waveform						 ///
////////////////////////////////////////////////////////////////////////////////////

	auto max_amp = max_element(std::begin(wfm) , std:: end(wfm));

	float avg = 0;

	for (unsigned int i=0; i<10 ; i++){
	  
		avg += wfm[i];
	}	

	avg /= 10;

	if (avg/ *max_amp > 0.25){

		TTL_Health_Flags["peak_at_start"] = 1;
	}
 	else{
		TTL_Health_Flags["peak_at_start"] = 0;
	}


////////////////////////////////////////////////////////////////////////////////////
/// The following segment breaks the waveform into four pieces and checks out    ///
/// the maximum amplitudes in each segment. It then does a basic chi squared     ///
/// analysis to see how well a simple constant fits the four peak points to help ///
/// identifysituations which have anomalously high, low/ or absent pulses. The	 ///
/// chi2/df is read out								 ///
////////////////////////////////////////////////////////////////////////////////////

	std::vector<T> wfm_tmp1;
	std::vector<T> wfm_tmp2;
	std::vector<T> wfm_tmp3;
	std::vector<T> wfm_tmp4;
	
	int start2 = floor(wfm.size()/4);
	int start3 = floor(wfm.size()/2);
	int start4 = floor(3*wfm.size()/4);

	for (unsigned int i = 0; i<=wfm.size()/4 ; i++) {
		wfm_tmp1.push_back(wfm[i]);
		wfm_tmp2.push_back(wfm[i+start2]);
		wfm_tmp3.push_back(wfm[i+start3]);
		wfm_tmp4.push_back(wfm[i+start4]);
	}

	auto sub_max_amp1 = max_element(std::begin(wfm_tmp1) , std:: end(wfm_tmp1));
	auto sub_max_amp2 = max_element(std::begin(wfm_tmp2) , std:: end(wfm_tmp2));
	auto sub_max_amp3 = max_element(std::begin(wfm_tmp3) , std:: end(wfm_tmp3));
	auto sub_max_amp4 = max_element(std::begin(wfm_tmp4) , std:: end(wfm_tmp4));
	
	avg = 1.0/4.0 *(*sub_max_amp1 + *sub_max_amp2 + *sub_max_amp3 + *sub_max_amp4);

	float chi =(std::pow((*sub_max_amp1 - avg),2) + std::pow((*sub_max_amp2 - avg),2)+std::pow((*sub_max_amp3 - avg),2) + std::pow((*sub_max_amp4 - avg),2))/(4.0*avg);

	TTL_Health_Flags["amp_chi"] = chi;
	
	return TTL_Health_Flags;

}

#endif
