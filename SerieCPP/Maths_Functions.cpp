//Maths_Functions.cpp
/*
Author: Lucien MAMAN
Date: 05/12/2016

Description: Implementation of the maths functions used in the construction of schemes
*/
#include "Maths_Functions.h"

//FUNCTION sign

double Maths_Functions::sign(double x)
{
	if (x<0.00){
		return -1.00;
	}
	if (x>0.00){
		return 1.00;
	}
	else return 0.00;
}