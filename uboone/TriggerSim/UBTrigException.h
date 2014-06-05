/**
 * \file UBTrigException.h
 *
 * \ingroup UBTriggerSim
 * 
 * \brief Class def header for exception classes in UBTriggerSim  package
 *
 * @author kazuhiro
 */

/** \addtogroup UBTriggerSim

    @{*/
#ifndef UBTRIGGEREXCEPTION_H
#define UBTRIGGEREXCEPTION_H

#include <string>
#include <exception>

namespace trigger {
  /**
     \class UBTrigException
     Simple exception class for UBTriggerSim
  */
  class UBTrigException : public std::exception{

  public:

    UBTrigException(std::string msg="") : std::exception(), _msg(msg)
    {}

    virtual ~UBTrigException() throw(){};
    virtual const char* what() const throw() 
    {return _msg.c_str(); }

  private:

    std::string _msg;
  };

}
#endif
/** @} */ // end of doxygen group 

