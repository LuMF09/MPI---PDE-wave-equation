#include "Errors.h"
#include <string>
#include <iostream>

//DivideZero to signal a division by 0

//CONSTRUCTOR
DivideZero::DivideZero(std::string a) : s("ERROR! Division by 0" + a){}

//FUNCTION WHAT()
char const* DivideZero::what() const throw() // display a personalized message
{
	return s.c_str();
}
