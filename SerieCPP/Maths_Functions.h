//Maths_Functions.h
#ifndef HEADER_MATHS_FUNCTIONS
#define HEADER_MATHS_FUNCTIONS
/*!
* \file Maths_Functions.h
* \brief Set of maths functions used in the construction of schemes
* \author Lucien MAMAN
* \date 05/12/2016
*/

class Maths_Functions
{
public:
	/*!
	*  \brief Const method used during the calculation of each point of both numerical and analytical solutions with the first set of equation
	*
	*  \param x : a double
	*
	*  \return -1 if x<0
	*  else if 0 x=0
	*  else 1
	*/
	static double sign(double x);
};
#endif