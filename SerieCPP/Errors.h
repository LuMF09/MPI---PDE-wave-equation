//Errors.h
/*!
* \file Errors.h
* \brief All potential errors which will be thrown in the program and catch un the main.
* \author Lucien MAMAN
* \date 05/12/2016
*/
#ifndef HEADER_ERRORS
#define HEADER_ERRORS
#include <exception>
#include <string>

/*! \class DivideZero
* \brief Class inherit from exception
*
*  Used to signal a division by 0
*/
class DivideZero : public std::exception
{
private:
	std::string s; /*!< string used to personalise each throw to make the debbugage easier */

public:

	/*!
	*  \brief Overloaded constructor of the class
	*
	*  \param a: the message to help debuggage
	*/
	DivideZero(std::string a);

	/*!
	*  \brief Virtual function used to display the error thrown
	*/
	virtual char const* what() const throw();

	virtual ~DivideZero() throw(){}
};
#endif
