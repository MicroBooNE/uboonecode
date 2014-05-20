/**
 * \file UBOpticalException.h
 *
 * \ingroup OpticalDetector
 * 
 * \brief Class def header for exception classes in OpticalDetector package
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

    @{*/
#ifndef UBOPTICALEXCEPTION_H
#define UBOPTICALEXCEPTION_H

#include <string>
#include <exception>

namespace opdet {
  /**
     \class UBOpticalException
     Simple exception class for OpticalDetector
  */
  class UBOpticalException : public std::exception{

  public:

    UBOpticalException(std::string msg="") : std::exception(), _msg(msg)
    {}

    virtual ~UBOpticalException() throw(){};
    virtual const char* what() const throw() 
    {return _msg.c_str(); }

  private:

    std::string _msg;
  };

}
#endif
/** @} */ // end of doxygen group 

