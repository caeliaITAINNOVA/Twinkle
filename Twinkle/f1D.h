//*******************************************************************************************//
//* TWINKLE																				 	*//
//* Version: 0.2																		 	*//
//* copyright (c) 2018 Instituto Tecnológico de Aragón (ITAINNOVA) (España)				 	*//
//* Date: 4 December 2018																 	*//
//* Author: Instituto Tecnológico de Aragón (ita@itainnova.es) - R.Rodríguez and V.Zambrano *//
//* All rights reserved																	 	*//
//*******************************************************************************************//

#pragma once
#ifndef H__F1D__H
#define H__F1D__H

#define _USE_MATH_DEFINES

#include <vector>
#include <armadillo>
#include <string>
#include <time.h>
#include <chrono>
#include <cmath>

// common namespace definition
namespace ROM
{
	class Twinkle;
	class Term;
	class f1D;
}

using namespace ROM;
using namespace arma;

extern std::string tlog;

// class f1D description
class ROM::f1D
{
	friend class Term;
	friend class Twinkle;

private:
	f1D(std::vector<double>, double, int);
	f1D();

	double evalf1D(double);
	vec getf1DinterData(double);
	vec setRandInit(double);
	vec setLinInit(double);

	vec m_discretization;
	vec m_zs;
};

#endif // !H__F1D__H
