//Schemes.h
#ifndef HEADER_SCHEMES
#define HEADER_SCHEMES
/*!
* \file Schemes.h
* \brief Structure of the different schemes to study a linear advection equation
* \author Lucien MAMAN
* \date 05/12/2016
*/

#include<vector>

/*! \class Schemes
* \brief Abstract class Schemes
*
*  Create and define different functions and attribute which will be use by the other inherit classes
*/
class Schemes
{
	protected:

		std::vector<double>x; /*!< Vector used to implement the different space positions  */
		std::vector<double> vANAL; /*!< Vector used to implement the analytical solution */
		std::vector<std::vector <double> > vNUM; /*!< Vector of vector used to implement the numerical solution */

		int nbIt; /*!< Number of iterations to go to the correct time step */

		double step;/*!< Time of simulation */
		double deltaT;/*!< Delta T (time), equals to (time of simulation)  / (grid Size - 1) */
		double deltaX;/*!< Delta X (space), equals to abs(Left boudary - right boundary)/ (grid Size - 1) */
		double cfl; /*!< Courant Friedrich Number, usefull to study the stability, equals to U*delta T/delta X */

	public:

		/*!
		*  \brief Constructor of the abstract class Schemes
		*
		*  \param step: time to run
		*  \param size: number of discretisation points
		*  \param cfl: the choosen CFL number
		*/
		Schemes(double step, int size, double cfl);

		//ACCESSORS

		/*!
		*  \brief Const Accessor method wich provide the size of the grid 
		*
		*  \return The size of the grid
		*/
		int getSizeGrid() const;

		/*!
		*  \brief Const Accessor method wich provide the time step
		*
		*  \return The time step 
		*/
		double getStep() const;

		//FUNCTIONS

		/*!
		*  \brief Const method used to calculate the analytical solution
		*/
		void calculAnalytical();

		/*!
		*  \brief Method used to calculate norms 1 and 2
		*/
		void calculNorms();

		/*!
		*  \brief Method used to intialize a vector of vector with the initial conditions
		*
		*  \param &v : a copy of a vector of vector
		*/
		void initNumerical(std::vector<std::vector <double> >&v);

		/*!
		*  \brief Method called in the main function to show the analytical solution, the numerical solution, norms and to know is the scheme is stable or not
		*/
		void showScheme();

		//PURE VIRTUAL FUNCTIONS

		/*!
		*  \brief Pure virtual function used to construct the numerical solutions
		*/
		virtual void calculNumerical() = 0;	

		/*!
		*  \brief Pure virtual function Test if a scheme is stable or not
		*
		*/
		virtual void isStable() = 0;
};

/*! \class ExplicitUpWind
* \brief Class inherit from Schemes class
*
*  Used to make numerical simulation of the explicit upwind scheme: FTBS
*/
class ExplicitUpWind : public Schemes
{

public:

	/*!
	*  \brief Constructor of the ExplicitUpWind scheme
	*
	*  \param step: time to run
	*  \param size: number of discretisation points
	*  \param cfl: the choosen CFL number
	*/
	ExplicitUpWind(double step, int size, double cfl);

	/*!
	*  \brief Definition of the pure virtual function: construct the numerical solution of this scheme
	*/
	void calculNumerical();

	/*!
	*  \brief Definition of the pure virtual function: Test if a scheme is stable or not
	*
	*/
	void isStable();
};

/*! \class ImplicitUpWind
* \brief Class inherit from Schemes class
*
*  Used to make numerical simulation of the implicite upwind scheme: FTBS
*/
class ImplicitUpWind : public Schemes
{
	public:

		/*!
		*  \brief Constructor of the ImplicitUpWind scheme
		*
		*  \param step: time to run
		*  \param size: number of discretisation points
		*  \param cfl: the choosen CFL number
		*/
		ImplicitUpWind(double step, int size, double cfl);

		/*!
		*  \brief Definition of the pure virtual function: construct the numerical solution of this scheme
		*/
		void calculNumerical();

		/*!
		*  \brief Definition of the pure virtual function: Test if a scheme is stable or not
		*
		*/
		void isStable();
};

/*! \class Crank_Nicolson
* \brief Class inherit from Schemes class
*
*  Used to make numerical simulation of the Crank_Nicolson scheme : ????
*/
class Crank_Nicolson : public Schemes
{

	public:
		/*!
		*  \brief Constructor of the Crank-Nicolson scheme
		*
		*  \param step: time to run
		*  \param size: number of discretisation points
		*  \param cfl: the choosen CFL number
		*/
		Crank_Nicolson(double step, int size, double cfl);

		/*!
		*  \brief Definition of the pure virtual function: construct the numerical solution of this scheme
		*/
		void calculNumerical();

		/*!
		*  \brief Definition of the pure virtual function: Test if a scheme is stable or not
		*
		*/
		void isStable();
};

#endif
