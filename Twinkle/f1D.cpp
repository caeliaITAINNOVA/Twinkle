//*******************************************************************************************//
//* TWINKLE																				 	*//
//* Version: 0.2																		 	*//
//* copyright (c) 2018 Instituto Tecnológico de Aragón (ITAINNOVA) (España)				 	*//
//* Date: 4 December 2018																 	*//
//* Author: Instituto Tecnológico de Aragón (ita@itainnova.es) - R.Rodríguez and V.Zambrano *//
//* All rights reserved																	 	*//
//*******************************************************************************************//

#include "f1D.h"

using namespace std::chrono;

std::string tlog = "";

// standard constructor for the f1D class, it creates one-dimensional functions using a specific initialization mode
f1D::f1D(std::vector<double> discret1D, double avg, int initMode)
{
	unsigned int Mk = discret1D.size();
	vec W1(Mk);
	W1.zeros();
	m_discretization = W1;
	for (unsigned int i = 0; i < Mk; i++)
	{
		m_discretization[i] = discret1D[i];
	}
	switch (initMode)
	{
		case 0:
			m_zs = setRandInit(avg); // random initialization
			break;
		case 1:
			m_zs = setLinInit(avg); // linear initialization with random slope
			break;
		default:
			m_zs = setRandInit(avg); // random initialization in case of external wrong usage of the discretization mode parameter
			tlog += "Wrong initialization mode selected, random initialization will be performed.\n";
			std::cout << "Wrong initialization mode selected, random initialization will be performed." << std::endl;
			break;
	}
}

// default constructor for the f1D class
f1D::f1D() {}

// set a random shape function initialization
vec f1D::setRandInit(double avg)
{
	unsigned int Mk = m_discretization.size();
	vec W2(Mk);
	W2.zeros();
	nanoseconds ms = duration_cast<nanoseconds>(system_clock::now().time_since_epoch());
	srand(ms.count());
	for (unsigned int i = 0; i < Mk; i++)
	{
		W2[i] = (double)rand() / RAND_MAX * avg;
	}
	return W2;
}

// set a random-slope linear shape function initialization
vec f1D::setLinInit(double avg)
{
	unsigned int Mk = m_discretization.size();
	vec W2(Mk);
	W2.zeros();
	nanoseconds ms = duration_cast<nanoseconds>(system_clock::now().time_since_epoch());
	srand(ms.count());
	double slope =  (double)rand() / RAND_MAX * M_PI;
	while (fabs(slope - M_PI / 2.0) < 0.1)
	{
		slope = (double)rand() / RAND_MAX * M_PI;
	}
	double b = avg - tan(slope) * Mk / 2.0;
	for (unsigned int i = 0; i < Mk; i++)
	{
		W2[i] = (double)i * tan(slope) + b;
	}
	return W2;
}

// evaluate the unidimensional function value at a given point
double f1D::evalf1D(double x)
{
	vec interpData = getf1DinterData(x);
	return interpData[1] * m_zs[interpData[0] - 1] + interpData[2] * m_zs[interpData[0]];
}

// get the interpolation data at a given value x of the independent variable
vec f1D::getf1DinterData(double x)
{
	// check that the value of x falls within the discretazation net's boundaries taking into account numerical approximation
	std::vector<double> num_precision {0.99999999999, 1.00000000001};

	double lconstraint = m_discretization[0] + 1; // inizialization of the precision interval's left boundary
	double rconstraint = m_discretization[m_discretization.size() - 1] + 1; // inizialization of the precision interval's right boundary

	if (x != 0) // redefine precision boundaries when x!=0
	{
		lconstraint = m_discretization[0] / x;
		rconstraint = m_discretization[m_discretization.size() - 1] / x;
	}
	else if (m_discretization[0] < 0 && m_discretization[m_discretization.size() - 1] <= 0)
	{
		lconstraint = m_discretization[m_discretization.size() - 1] + 1;
		rconstraint = m_discretization[0] + 1;
	}

	if ((m_discretization[0] >= 0 && m_discretization[m_discretization.size() - 1] > 0))
	{
		if (lconstraint > num_precision[1] || rconstraint < num_precision[0])
		{
			return 0;
		}
	}
	else if (m_discretization[0] < 0 && m_discretization[m_discretization.size() - 1] <= 0)
	{
		if (lconstraint < num_precision[0] || rconstraint > num_precision[1])
		{
			return 0;
		}
	}
	else
	{
		if (x >= 0)
		{
			if (lconstraint > num_precision[1] || rconstraint < num_precision[0])
			{
				return 0;
			}
		}
		else
		{
			if (lconstraint < num_precision[0] || rconstraint > num_precision[1])
			{
				return 0;
			}
		}
	}

	// search for the segment in the discretization net in which the value x falls
	int i = 0;
	for (unsigned int j = 1; j < m_discretization.size(); j++)
	{
		i = j;
		if (((x > m_discretization[j - 1]) && (x <= m_discretization[j])) || (x == m_discretization[0]))
		{
			break;
		}
	}

	// normalized value of x
	double xp =  i + (x - m_discretization[i - 1]) / (m_discretization[i] - m_discretization[i - 1]) - 1;
	// value of the shape function evaluated at the input point's independent variable (corresponding to the discretization net's point i-1)
	double yr = xp - i + 1;
	// value of the shape function evaluated at the input point's independent variable (corresponding to the discretization net's point i)
	double yl = -xp + i;

	vec interpData(3);
	interpData[0] = i;
	interpData[1] = yl;
	interpData[2] = yr;

	return interpData;
}
