//*******************************************************************************************//
//* TWINKLE																				 	*//
//* Version: 0.2																		 	*//
//* copyright (c) 2018 Instituto Tecnológico de Aragón (ITAINNOVA) (España)				 	*//
//* Date: 4 December 2018																 	*//
//* Author: Instituto Tecnológico de Aragón (ita@itainnova.es) - R.Rodríguez and V.Zambrano *//
//* All rights reserved																	 	*//
//*******************************************************************************************//

#include "Term.h"

// standard constructor for the Term class, it creates terms using a specific initialization mode
Term::Term(std::vector<std::vector<double>> discret, double avg, int initMode)
{
	m_alpha = 1;
	for (unsigned int i = 0; i < discret.size(); i++)
	{
		f1D newf1D(discret[i], avg, initMode);
		m_fnD.push_back(newf1D);
	}
}

// default constructor for the Term class
Term::Term() {}

// evaluate the ROM's term at the values of the independent variables in xs
double Term::evalTerm(std::vector<double> xs)
{
	double yterm = m_alpha;

	for (unsigned int i = 0; i < m_fnD.size(); i++)
	{
		yterm *= m_fnD[i].evalf1D(xs[i]);
	}

	return yterm;
}

// evaluate the ROM's term (without the weighting factor alpha) at the values of the independent variables in xs
double Term::evalUnweightTerm(std::vector<double> xs)
{
	double yterm = 1;

	for (unsigned int i = 0; i < m_fnD.size(); i++)
	{
		yterm *= m_fnD[i].evalf1D(xs[i]);
	}

	return yterm;
}
