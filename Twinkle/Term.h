//*******************************************************************************************//
//* TWINKLE																				 	*//
//* Version: 0.2																		 	*//
//* copyright (c) 2018 Instituto Tecnológico de Aragón (ITAINNOVA) (España)				 	*//
//* Date: 4 December 2018																 	*//
//* Author: Instituto Tecnológico de Aragón (ita@itainnova.es) - R.Rodríguez and V.Zambrano *//
//* All rights reserved																	 	*//
//*******************************************************************************************//

#pragma once
#ifndef H__TERM__H
#define H__TERM__H

#include "f1D.h"

// class Term description
class ROM::Term
{
	friend class Twinkle;

private:
	Term(std::vector<std::vector<double>>, double, int);
	Term();

	double evalTerm(std::vector<double>);
	double evalUnweightTerm(std::vector<double>);

	double m_alpha;
	std::vector<f1D> m_fnD;
};

#endif // !H__TERM__H
