//*******************************************************************************************//
//* TWINKLE																				 	*//
//* Version: 0.2																		 	*//
//* copyright (c) 2018 Instituto Tecnológico de Aragón (ITAINNOVA) (España)				 	*//
//* Date: 4 December 2018																 	*//
//* Author: Instituto Tecnológico de Aragón (ita@itainnova.es) - R.Rodríguez and V.Zambrano *//
//* All rights reserved																	 	*//
//*******************************************************************************************//

#include "Twinkle.h"

// standard constructor for the Twinkle class
Twinkle::Twinkle()
{
	cube W(2,2,2);
	W.zeros();
	m_Ntensor = W;
	m_Tol = 0;
	m_bConverg = true;
	m_globalTol = 1E-3;
	m_termTol = 1E-4;
	m_maxIter = 150;
	m_maxZIter = 5;
	m_errorMin = 1;
	m_precision = 6;
	m_rows = 1;
	m_cols = 1;
}

// overloaded constructor for the Twinkle class, it fills member variables from a previously computed ROM file
Twinkle::Twinkle(std::string filename)
{
	// initialize member variables
	cube W(2,2,2);
	W.zeros();
	m_Ntensor = W;
	m_Tol = 0;
	m_bConverg = true;
	m_globalTol = 1E-3;
	m_termTol = 1E-4;
	m_maxIter = 150;
	m_maxZIter = 5;
	m_errorMin = 1;
	m_precision = 6;
	m_rows = 1;
	m_cols = 1;
	//

	std::ifstream file;
	file.open(filename);

	try
	{
		std::string line;
		int T;
		int N;
		std::vector<std::vector<int>> D;

		line.clear();
		std::getline(file, line);
		std::getline(file, line);
		line.clear();
		std::getline(file, line);
		std::stringstream lineStream(line);
		std::string value;
		std::getline(lineStream, value, '\t');
		std::getline(lineStream, value, '\t');
		T = std::stoi(value);
		std::getline(file, line);
		std::stringstream lineStream2(line);
		std::getline(lineStream2, value, '\t');
		std::getline(lineStream2, value, '\t');
		N = std::stoi(value);
		for (int i = 0; i < N; i++)
		{
			std::getline(file, line);
			std::stringstream lineStream3(line);
			std::getline(lineStream3, value, '\t');
			std::vector<int> di;
			for (int j = 0; j < T; j++)
			{
				std::getline(lineStream3, value, '\t');
				di.push_back(std::stoi(value));
			}
			D.push_back(di);
		}
		std::getline(file, line);
		std::getline(file, line);
		m_tmatrix.clear();
		for (int i = 0; i < T; i++)
		{
			Term Ti;
			std::getline(file, line);
			std::getline(file, line);
			std::stringstream lineStream4(line);
			std::getline(lineStream4, value, '\t');
			std::getline(lineStream4, value, '\t');
			Ti.m_alpha = std::stod(value);
			for (int j = 0; j < N; j++)
			{
				std::getline(file, line);
				f1D fij;
				vec fij_disc(D[j][i]);
				fij_disc.zeros();
				vec fij_zs(D[j][i]);
				fij_zs.zeros();
				for (int k = 0; k < D[j][i]; k++)
				{
					std::getline(file, line);
					std::stringstream lineStream4(line);
					std::getline(lineStream4, value, '\t');
					fij_disc[k] = std::stod(value);
					std::getline(lineStream4, value, '\t');
					fij_zs[k] = std::stod(value);
				}
				fij.m_discretization = fij_disc;
				fij.m_zs = fij_zs;
				Ti.m_fnD.push_back(fij);
			}
			m_tmatrix.push_back(Ti);
		}
		file.close();

		tlog += "ROM file: reading completed!\n";
		std::cout << "ROM file: reading completed!" << std::endl;
	}
	catch (...)
	{
		std::cout << "Please check ROM file format, it must be UTF-8." << std::endl;
		throw;
	}
}


// read the data from filename and assign it to the member variable m_data
void Twinkle::setData(std::string filename)
{
	std::ifstream file;
	file.open(filename);

	try
	{
		std::vector<std::vector<double>> edata;
		std::string line;
		while(!file.eof())
		{
			line.clear();
			std::getline(file, line);
			if (line != "" && line != "\r" && line != "\n")
			{
				std::stringstream lineStream(line);
				std::string value;
				std::vector<double> row;
				while(std::getline(lineStream, value, ';'))
				{
					row.push_back(std::stod(value));
				}
				edata.push_back(row);
			}
		}
		file.close();


		tlog += "Data file : reading completed!\n";
		std::cout << "Data file: reading completed!" << std::endl;

		m_data = edata;

		// number of points of the input data
		m_rows = m_data.size();

		// number of independent variables
		m_cols = m_data[0].size() - 1;
	}
	catch (...)
	{
		std::cout << "Please check data file format, it must be UTF-8." << std::endl;
		throw;
	}
}

// set data for rom calculation from a previously obtained matrix
void Twinkle::setData(std::vector<std::vector<double>> data)
{
	m_data = data;

	// number of points of the input data
	m_rows = m_data.size();

	// number of independent variables
	m_cols = m_data[0].size() - 1;
}

// set user parameters
void Twinkle::setParams(double globalTol, double termTol, int maxIter, int maxZIter, int precision)
{
	m_globalTol = globalTol; // global normalized tolerance value of the overall ROM error. the calculation of additional terms stops when the overall error is smaller than the global tolerance.
	m_termTol = termTol; // tolerance of the alternate least squares (ALS) computed for each term separately in the calculation of the shape functions. if the change in succesive iterations of the ALS method is lower than this tolerance, the process stops.
	m_maxIter = maxIter; // maximum number of initializations attempts of the shape functions if the convergence of the ALS method is not achieved.
	m_maxZIter = maxZIter; // maximum number of iterations in the ALS method.
	m_precision = precision;
}

// generate a homogeneous discretization net
void Twinkle::genDiscret(int nPoints)
{
	std::vector<std::vector<double>> output;
	std::vector<std::vector<double>> bounds; // minima and maxima of the independent variables
	bounds = setBounds();

	// calculation of the discretization nets
	for (unsigned int j = 0; j < m_data[0].size() - 1; j++)
	{
		std::vector<double> discretj;
		for (int i = 0; i < nPoints; i++)
		{
			if (nPoints > 1)
			{
				discretj.push_back(bounds[j][0] + i * (bounds[j][1] - bounds[j][0]) / (nPoints - 1));
			}
			else
			{
				tlog += "The number of discretization points can not be less than two.\n";
				std::cout << "The number of discretization points can not be less than two." << std::endl;
				return;
			}
		}
		output.push_back(discretj);
	}

	tlog += std::string("Discretization completed!\n") + std::string("data matrix: ") + std::to_string(m_data.size()) + std::string(" x ") + std::to_string(m_data[0].size()) + std::string("\n");
	std::cout << "Discretization completed!" << std::endl;
	std::cout << "data matrix: " << m_data.size() << " x " << m_data[0].size() << std::endl;

	m_discret = output;;
}

// generate a variable-wise homogeneous discretization net from a vector of values
void Twinkle::genDiscret(std::vector<int> nPoints)
{
	std::vector<std::vector<double>> output;
	std::vector<std::vector<double>> bounds; // minima and maxima of the independent variables
	bounds = setBounds();

	// calculation of the discretization nets
	for (unsigned int j = 0; j < m_data[0].size() - 1; j++)
	{
		std::vector<double> discretj;
		// the user can set an exact discretization for variable j, by setting nPoints[j] = 0 in the vector
		if (nPoints[j] == 0)
		{
			output.push_back(genExactDiscret(j));
		}
		else
		{
			for (int i = 0; i < nPoints[j]; i++)
			{
				discretj.push_back(bounds[j][0] + i * (bounds[j][1] - bounds[j][0]) / (nPoints[j] - 1));
			}
			output.push_back(discretj);
		}
	}

	tlog += std::string("Discretization completed!\n") + std::string("data matrix: ") + std::to_string(m_data.size()) + std::string(" x ") + std::to_string(m_data[0].size()) + std::string("\n");
	std::cout << "Discretization completed!" << std::endl;
	std::cout << "data matrix: " << m_data.size() << " x " << m_data[0].size() << std::endl;

	m_discret = output;;
}

// generate a homogeneous discretization net from a file
void Twinkle::genDiscret(std::string filename)
{
	std::ifstream file;
	file.open(filename);

	try
	{
		std::vector<int> discret;
		std::string line;
		// create a vector of integers with the number of points of the discretization of each independent variable
		while(!file.eof())
		{
			line.clear();
			std::getline(file, line);
			if (line != "" && line != "\r" && line != "\n")
			{
				std::stringstream lineStream(line);
				std::string value;
				while(getline(lineStream, value))
				{
					discret.push_back(std::stoi(value));
				}
			}
		}
		file.close();

		tlog += "Discretization file: reading completed!\n";
		std::cout << "Discretization file: reading completed!" << std::endl;

		// call void Twinkle::genDiscret(std::vector<int> nPoints) with the previously created vector
		genDiscret(discret);
	}
	catch (...)
	{
		std::cout << "Please check discretization file format, it must be UTF-8." << std::endl;
		throw;
	}
}

// generate a non-homogeneous discretization net from a file
void Twinkle::genNonHomogDiscret(std::string filename)
{
	std::ifstream file;
	file.open(filename);

	try
	{
		std::string line;
		// fill the m_discret matrix with the double values of each independent variable
		m_discret.clear();
		while(!file.eof())
		{
			line.clear();
			std::getline(file, line);
			if (line != "" && line != "\r" && line != "\n")
			{
				std::stringstream lineStream(line);
				std::string value;
				std::vector<double> row;
				while(std::getline(lineStream, value, ';'))
				{
					row.push_back(std::stod(value));
				}
				m_discret.push_back(row);
			}
		}
		file.close();

		tlog += "Discretization file: reading completed!\n";
		std::cout << "Discretization file: reading completed!" << std::endl;
	}
	catch (...)
	{
		std::cout << "Please check discretization file format, it must be UTF-8." << std::endl;
		throw;
	}
}

// generate the exact discretization net for one variable
std::vector<double> Twinkle::genExactDiscret(int j)
{
	std::vector<double> output;
	// obtain all different values of the desired independent variable (j)
	for (unsigned int i = 0; i < m_data.size(); i++)
	{
		if (!containsValue(output, m_data[i][j]))
		{
			output.push_back(m_data[i][j]);
		}
	}
	// sort in increasing order all different values of the desired independent variable (j)
	for (unsigned int j = 0; j < m_data[0].size() - 1; j++)
	{
		std::sort(output.begin(), output.end());
	}

	return output;
}

// generate the exact discretization net for each variable
void Twinkle::genExactDiscret()
{
	std::vector<std::vector<double>> output;

	for (unsigned int j = 0; j < m_data[0].size() - 1; j++)
	{
		std::vector<double> discretj;
		output.push_back(discretj);
	}
	// obtain all different values of each independent variable
	for (unsigned int i = 0; i < m_data.size(); i++)
	{
		for (unsigned int j = 0; j < m_data[0].size() - 1; j++)
		{
			if (!containsValue(output[j], m_data[i][j]))
			{
				output[j].push_back(m_data[i][j]);
			}
		}
	}
	// sort in increasing order all different values of each independent variable
	for (unsigned int j = 0; j < m_data[0].size() - 1; j++)
	{
		std::sort(output[j].begin(), output[j].end());
	}

	tlog += std::string("Discretization completed!\n") + std::string("data matrix: ") + std::to_string(m_data.size()) + std::string(" x ") + std::to_string(m_data[0].size()) + std::string("\n");
	std::cout << "Discretization completed!" << std::endl;
	std::cout << "data matrix: " << m_data.size() << " x " << m_data[0].size() << std::endl;
	m_discret = output;
}

// set the maxima and minima from the dataset values to be used in the discretization net
std::vector<std::vector<double>> Twinkle::setBounds()
{
	std::vector<std::vector<double>> bounds; // minima and maxima of the independent variables
	for (unsigned int j = 0; j < m_data[0].size() - 1; j++)
	{
		std::vector<double> dims;
		dims.push_back(m_data[0][j]);
		dims.push_back(m_data[0][j]);
		bounds.push_back(dims);
	}

	for (unsigned int i = 0; i < m_data.size(); i++)
	{
		for (unsigned int j = 0; j < m_data[0].size() - 1; j++)
		{
			if (m_data[i][j] < bounds[j][0])
			{
				bounds[j][0] = m_data[i][j];
			}
			if (m_data[i][j] > bounds[j][1])
			{
				bounds[j][1] = m_data[i][j];
			}
		}
	}
	return bounds;
}

// check if value is contained in vector vec
bool Twinkle::containsValue(std::vector<double> vec, double value)
{
	for (unsigned int i = 0; i < vec.size(); i++)
	{
		if (vec[i] == value) return true;
	}
	return false;
}

// calculate the ROM
void Twinkle::calcROM(int maxTerm, int initMode)
{
	if (!m_data.empty() && !m_discret.empty() && isDataInDiscret()) // make sure there are data to work on and the discretization net is correct
	{
		tlog += "ROM calculation started.\n";
		std::cout << "ROM calculation started." << std::endl;

		m_bConverg = true; // overall convergence
		m_tmatrix.clear();
		std::vector<Term> tmatrix; // temporal list of Terms
		int currTerm = 0;
		double div = 0; // value used for the normalization of the error
		m_currEval.clear();
		// initial evaluation of the ROM and calculation of div
		for (int m = 0; m < m_rows; m++)
		{
			m_currEval.push_back(0.0);
			div += std::pow(m_data[m][m_cols], 2);
		}
		m_Tol = std::sqrt(div) / m_rows; // value used for the initializaton of the shape functions

		calcZs(tmatrix, currTerm, initMode); // calculation of the shape function of Term number currTerm

		// alpha calculation, special for first term
		double alpha = 1;
		for (int dim = 0; dim < m_cols; dim++)
		{
			double norm2 = norm(m_rcurr[dim], 2);
			alpha *= norm2;
			vec z_norm = m_rcurr[dim] / norm2;
			m_tmatrix[currTerm].m_fnD[dim].m_zs = z_norm;
		}
		m_tmatrix[currTerm].m_alpha = alpha;

		tlog += std::string("Term 1 alpha = ") + std::to_string(alpha) + std::string("\n");
		std::cout << "Term 1 alpha = " << alpha << std::endl;

		double globalError = fillCurrEval(div); // fill new evaluation of the ROM and return current error

		// calculate new terms until convergence
		currTerm = 1;
		int numRecTerm = 0;

		// loop for the calculation of the next terms
		while (globalError > m_globalTol && currTerm < maxTerm && numRecTerm < 100) // the loop stops when at least one of the following constraints is false: the globalError is higher than the required threshold, the number of terms is lower than the selected one, the number of times a term is recalculated (when globalError is not reduced) is lower than 100
		{
			m_errorMin = 1e10; // minimum error obtained in the ALS method with the different initializations

			calcZs(tmatrix, currTerm, initMode); // calculate current term's shape functions

			normZs(currTerm); // normalize the obtained shape functions

			std::vector<double> oldAlphas; // save previous alphas
			for (int i = 0; i < currTerm; i++)
			{
				oldAlphas.push_back(tmatrix[i].m_alpha);
			}

			calcAlphas(); // calculate new alphas

			double newGlobalError = fillCurrEval(div); // fill new evaluation of the ROM and return current error
			// check the error is reduced
			if (newGlobalError < globalError)
			{
				// if the error is reduced the term is accepted
				globalError = newGlobalError;
				currTerm++;
				numRecTerm = 0;
			}
			else
			{
				// if the error is not reduced the term is rejected and recalculated
				m_tmatrix.pop_back();

				tmatrix.pop_back();
				setAlphas(tmatrix, oldAlphas); // restore old alphas
				fillCurrEval(div, tmatrix); // recalculate the evaluation of the ROM at the input points

				numRecTerm++;

				if (numRecTerm == 100)
				{
					tlog += std::string("Term ") + std::to_string(currTerm + 1) + std::string(" was removed due to increased global error and the ROM calculation has been stopped.\n");
					std::cout << "Term " << currTerm + 1 << " was removed due to increased global error and the ROM calculation has been stopped." << std::endl;
				}
				else
				{
					tlog += std::string("Term ") + std::to_string(currTerm + 1) + std::string(" removed due to increased global error, a new attempt will be performed.\n");
					std::cout << "Term " << currTerm + 1 << " removed due to increased global error, a new attempt will be performed." << std::endl;
				}
			}
		}
		tlog += std::string("Overall z convergence achieved? ") + std::string(m_bConverg ? "true" : "false") + std::string("\n") + std::string("Number of terms: ") + std::to_string(m_tmatrix.size()) + std::string("\n");
		std::cout << "Overall z convergence achieved? " << std::boolalpha << m_bConverg << std::endl;
		std::cout << "Number of terms: " << m_tmatrix.size() << std::endl;
		for (unsigned int iTerm = 0; iTerm < m_tmatrix.size(); iTerm++)
		{
			tlog += std::string("Term ") + std::to_string(iTerm + 1) + std::string("\n") + std::string("Alpha = ") + std::to_string(m_tmatrix[iTerm].m_alpha) + std::string("\n");
			std::cout << "Term " << iTerm + 1 << std::endl;
			std::cout << "Alpha = " << m_tmatrix[iTerm].m_alpha << std::endl;
		}

		tlog += "ROM calculation finished.\n";
		std::cout << "ROM calculation finished." << std::endl;
	}
	else
	{
		if (!isDataInDiscret())
		{
			tlog += "The maxima and minima of the set discretization net does not match the data, please check it and try again.\n";
			std::cout << "The maxima and minima of the set discretization net does not match the data, please check it and try again." << std::endl;
		}
		else
		{
			tlog += "Can not execute calcROM, please read data first: setData function and set a discretization net.\n";
			std::cout << "Can not execute calcROM, please read data first: setData function and set a discretization net." << std::endl;
		}
		return;
	}
}

// check that both data and discretization net have the same minima and maxima
bool Twinkle::isDataInDiscret()
{
	std::vector<std::vector<double>> bounds; // minima and maxima of the independent variables
	bounds = setBounds();
	double eps = 0.00000000001;
	for (int i = 0; i < (int)m_discret.size(); i++)
	{
		if ((fabs(bounds[i][0] - m_discret[i][0]) > eps) || (fabs(bounds[i][1] - m_discret[i][m_discret[i].size() - 1]) > eps))
		{
			return false;
		}
	}
	return true;
}

// compute the Zs
void Twinkle::calcZs(std::vector<Term> &tmatrix, int currTerm, int initMode)
{
	bool bConverg = false; // control variable for the convergence of the shape functions of the current term
	int conv_count = 0;
	std::vector<vec> r_curr; // vector of the shape functions calculated in the current iteration of the ALS method
	std::vector<vec> r_prev; // vector of the shape functions calculated in the previous iteration of the ALS method
	while (bConverg == false && conv_count < m_maxZIter) // iterate until convergence is achieved for each term separately and repeat its calculation m_maxZIter times otherwise
	{
		Term termi(m_discret, pow(m_Tol, 1.0 / m_cols), initMode); // new instance of Term using a pow(m_Tol, 1.0 / m_cols) average value, an initMode initialization mode and a m_discret discretization net
		// add or substitute the new instance in the temporary list of terms
		if ((int)tmatrix.size() - 1 == currTerm)
		{
			tmatrix[currTerm] = termi;
		}
		else
		{
			tmatrix.push_back(termi);
			fillNtensor(currTerm, tmatrix); // required only if terms have different discretization
		}
		// set r_prev with the initializations of the shape functions
		r_prev.clear();
		for (int dim = 0; dim < m_cols; dim++)
		{
			r_prev.push_back(tmatrix[currTerm].m_fnD[dim].m_zs);
		}

		r_curr = getrcurr(tmatrix[currTerm]); // first approximation of the shape functions of all dimensions
		double error = getMaxrDiff(r_curr, r_prev); // calculate error between consecutive iterations of the ALS method
		int iter = 0;

		// iteration until convergence of all shape functions of all dimensions
		while (error > m_termTol && iter < m_maxIter)
		{
			r_prev = copyRs(r_curr);
			r_curr = getrcurr(tmatrix[currTerm]);
			error = getMaxrDiff(r_curr, r_prev);
			iter++;
		}

		bConverg = (error <= m_termTol); // check the convergence of the current term

		// check if the error obtained from the ALS method is lower than the desired threshold, if so copy the current term from the temporary matrix and store it into m_tmatrix. use the shape function calculation which gave the least error.
		if (error <= m_errorMin)
		{
			m_errorMin = error;
			if ((int)m_tmatrix.size() - 1 == currTerm)
			{
				m_tmatrix[currTerm] = tmatrix[currTerm];
			}
			else
			{
				m_tmatrix.push_back(tmatrix[currTerm]);
			}
			m_rcurr = copyRs(r_curr);
		}

		conv_count++;

		tlog += std::string("Convergence achieved in iteration ") + std::to_string(conv_count) +std::string("? ") + std::string(bConverg ? "true" : "false") + std::string(". Error = ") + std::to_string(error) + std::string("\n");
		std::cout << "Convergence achieved in iteration " << conv_count << "? " << std::boolalpha << bConverg << ". Error = " << error << std::endl;
	}
	m_bConverg = (m_bConverg && bConverg); // check the convergence of all terms up to the current term

	tlog += std::string("Final z error: ") + std::to_string(m_errorMin) + std::string("\n") + std::string("Tolerance value: ") + std::to_string(m_termTol) + std::string("\n");
	std::cout << "Final z error: " << m_errorMin << std::endl;
	std::cout << "Tolerance value: " << m_termTol << std::endl;
}

// fill a tensor with the data used for the interpolation of the shape functions, for each input point
void Twinkle::fillNtensor(int nTerm, std::vector<Term> matrix)
{
	if (m_data.size() == 0 || nTerm > (int)matrix.size() - 1)
	{
		return;
	}
	cube W(m_rows, m_cols, 3);
	W.zeros();
	m_Ntensor = W;
	for (int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_cols; j++)
		{
			vec iData(3);
			iData.zeros();
			iData = matrix[nTerm].m_fnD[j].getf1DinterData(m_data[i][j]);
			m_Ntensor(i, j, 0) = iData(0); // index of the point in the discretization net
			m_Ntensor(i, j, 1) = iData(1); // value of the shape function evaluated at the input point's independent variable (corresponding to the discretization net's point iData(0)-1)
			m_Ntensor(i, j, 2) = iData(2); // value of the shape function evaluated at the input point's independent variable (corresponding to the discretization net's point iData(0))
		}
	}
}

// compute the shape functions through the ALS method
std::vector<vec> Twinkle::getrcurr(Term& term)
{
	std::vector<vec> r;
	for (int dim = 0; dim < m_cols; dim++)
	{
		int Mk = term.m_fnD[dim].m_discretization.size(); // number of points in the discretization net of variable dim
		mat K1(Mk, Mk);
		K1.zeros();
		vec p1(Mk);
		p1.zeros();
		// calculation of the matrix K1 and constants vector p1 for the ALS method corresponding to variable dim
		for (int m = 0; m < m_rows; m++)
		{
			double bm = 1;
			for (int k = 0; k < m_cols; k++)
			{
				if (k != dim)
				{
					bm *= term.m_fnD[k].evalf1D(m_data[m][k]);
				}
			}
			unsigned int position = (unsigned int)m_Ntensor(m, dim, 0);
			double Ni =  m_Ntensor(m, dim, 1);
			double Nj =  m_Ntensor(m, dim, 2);
			if (Ni == 0)
			{
				K1(position, position) += std::pow(bm, 2) * Nj * Nj;
				p1(position) += bm * (m_data[m][m_cols] - m_currEval[m]) * Nj;
			}
			else if (Nj == 0)
			{
				K1(position - 1, position - 1) += std::pow(bm, 2) * Ni * Ni;
				p1(position - 1) += bm * (m_data[m][m_cols] - m_currEval[m]) * Ni;
			}
			else
			{
				K1(position - 1, position - 1) += std::pow(bm, 2) * Ni * Ni;
				K1(position - 1, position) += std::pow(bm, 2) * Ni * Nj;
				K1(position, position - 1) += std::pow(bm, 2) * Nj * Ni;
				K1(position, position) += std::pow(bm, 2) * Nj * Nj;
				p1(position - 1) += bm * (m_data[m][m_cols] - m_currEval[m]) * Ni;
				p1(position) += bm * (m_data[m][m_cols] - m_currEval[m]) * Nj;
			}
		}
		vec ri(Mk);
		ri.zeros();
		// solve the ALS method using Armadillo
		ri = solve(K1, p1, solve_opts::equilibrate);
		r.push_back(ri);
		term.m_fnD[dim].m_zs = ri;
	}
	return r;
}

// return a copy of the vector of shape functions (origin)
std::vector<vec> Twinkle::copyRs(std::vector<vec> origin)
{
	std::vector<vec> destiny;
	int originSize = origin.size();
	for (int i = 0; i < originSize; i++)
	{
		int Rsize = origin[i].size();
		vec Ri(Rsize);
		Ri.zeros();
		for (int j = 0; j < Rsize; j++)
		{
			Ri(j) = origin[i](j);
		}
		destiny.push_back(Ri);
	}
	return destiny;
}

// return the maximum of the normalization of the norm 2 of the differences of the vectors in r_curr and r_prev
double Twinkle::getMaxrDiff(std::vector<vec> r_curr, std::vector<vec> r_prev)
{
	double maxDiff =  norm(r_curr[0] - r_prev[0], 2) / norm(r_prev[0], 2);
	for (unsigned int i = 1; i < r_curr.size(); i++)
	{
		double diff = norm(r_curr[i] - r_prev[i], 2) / norm(r_prev[i], 2);
		if (diff > maxDiff)
		{
			maxDiff = diff;
		}
	}

	return maxDiff;
}

// normalize the shape functions of term currTerm
void Twinkle::normZs(int currTerm)
{
	for (int dim = 0; dim < m_cols; dim++)
	{
		m_tmatrix[currTerm].m_fnD[dim].m_zs = m_rcurr[dim] / norm(m_rcurr[dim], 2);
	}
}

// calculate the alphas
void Twinkle::calcAlphas()
{
	int T = m_tmatrix.size();
	mat K(T, T);
	K.zeros();
	vec f(T);
	f.zeros();
	// calculation of the matrix K and constants vector f for the LS method
	for (int m = 0; m < m_rows; m++)
	{
		for (int i = 0; i < T; i++)
		{
			for (int j = 0; j < T; j++)
			{
				K(i, j) += (m_tmatrix[i].evalUnweightTerm(m_data[m])) * (m_tmatrix[j].evalUnweightTerm(m_data[m]));
			}
			f(i) += (m_tmatrix[i].evalUnweightTerm(m_data[m])) * m_data[m][m_cols];
		}
	}
	// solve the LS method using Armadillo
	vec alphas = solve(K, f, solve_opts::equilibrate);

	// alphas assignment
	for (int i = 0; i < T; i++)
	{
		tlog += std::string("Term ") + std::to_string(i + 1) + std::string(" alpha = ") + std::to_string(alphas[i]) + std::string("\n");
		std::cout << "Term " << i + 1 << " alpha = " << alphas[i] << std::endl;

		m_tmatrix[i].m_alpha = alphas[i];
	}
}

// fill the vector m_currEval with the current evaluation of the ROM at the input points and return the global approximation error
double Twinkle::fillCurrEval(double div)
{
	m_currEval.clear();
	double E = 0;
	m_Tol = 0;
	double globalError = 0;
	for (int m = 0; m < m_rows; m++)
	{
		m_currEval.push_back(evalROM(m_data[m]));

		// error calculation
		double diff = m_data[m][m_cols] - m_currEval[m];
		E += std::pow(diff, 2);
		m_Tol += fabs(diff); // value used in the initialization of the shape functions of the next term
	}
	m_Tol /= m_rows;
	globalError = std::sqrt(E / div);

	tlog += std::string("Global approximation error: ") + std::to_string(globalError) + std::string("\n");
	std::cout << "Global approximation error: " << globalError << std::endl;

	return globalError;
}

// fill the vector m_currEval with the current evaluation of the ROM at the input points, using tmatrix
void Twinkle::fillCurrEval(double div, std::vector<Term> tmatrix)
{
	m_currEval.clear();
	double E = 0;
	m_Tol = 0;
	for (int m = 0; m < m_rows; m++)
	{
		double out = 0;
		for (unsigned int k = 0; k < tmatrix.size(); k++)
		{
			out += tmatrix[k].evalTerm(m_data[m]);
		}
		m_currEval.push_back(out);

		// error calculation
		double diffAvg = m_data[m][m_cols] - m_currEval[m];
		E += std::pow(diffAvg, 2);
		m_Tol += fabs(diffAvg); // value used in the initialization of the shape functions of the next term
	}
	m_Tol /= m_rows;
}

// restore alphas in tmatrix with the values in oldAlphas
void Twinkle::setAlphas(std::vector<Term> tmatrix, std::vector<double> oldAlphas)
{
	for (unsigned int i = 0; i < oldAlphas.size(); i++)
	{
		tmatrix[i].m_alpha = oldAlphas[i];
	}
}

// write a file with the ROM's parameters
void Twinkle::writeROM(std::string fileName)
{
	if(!m_data.empty())
	{
		std::ofstream file(fileName);
		file << "***********************************************************************" << "\n";
		file << "*** ROM summary" << "\n";
		file << "*** Number of terms:\t" <<  m_tmatrix.size() << "\n";
		file << "*** Number of dimensions:\t" <<  m_cols << "\n";
		for (int i = 0; i < m_cols; i++)
		{
			file << "*** Dimension " << i + 1 << " discretization:\t";
			for (unsigned int j = 0; j < m_tmatrix.size(); j++)
			{
				file << m_tmatrix[j].m_fnD[i].m_discretization.size() << "\t";
			}
			file << "\n";
		}
		file << "***********************************************************************" << "\n";
		file << "ROM Data:" << "\n";
		for (unsigned int i = 0; i < m_tmatrix.size(); i++)
		{
			file << "Term " << (i + 1) << "\n";
			std::stringstream alpha;
			alpha.precision(m_precision);
			alpha << m_tmatrix[i].m_alpha;
			file << "Alpha:\t" << alpha.str() << "\n";
			for (unsigned int j = 0; j < m_tmatrix[i].m_fnD.size(); j++)
			{
				file << "Dimension " << (j + 1) << "\n";
				for (unsigned int k = 0; k < m_tmatrix[i].m_fnD[j].m_discretization.size(); k++)
				{
					std::stringstream xdisc;
					xdisc.precision(m_precision);
					xdisc << m_tmatrix[i].m_fnD[j].m_discretization[k];
					std::stringstream ydisc;
					ydisc.precision(m_precision);
					ydisc << m_tmatrix[i].m_fnD[j].m_zs[k];
					file << xdisc.str() << "\t" << ydisc.str() << "\n";
				}
			}
		}
		file.close();

		tlog += "ROM file: writing completed!\n";
		std::cout << "ROM file: writing completed!" << std::endl;
	}
	else
	{
		tlog += "Can not execute writeROM, please read data first: setData function.\n";
		std::cout << "Can not execute writeROM, please read data first: setData function." << std::endl;
	}

	return;
}

// write a file (fileName) with the ROM's prediction evaluated at the input points' values
void Twinkle::writePrediction(std::string fileName)
{
	if ((!m_data.empty()) && (m_tmatrix.size() > 0))
	{
		std::ofstream fileP(fileName);
		fileP << "data" << "\t" << "prediction" << "\n";
		for (int i = 0; i < m_rows; i++)
		{
			std::stringstream odata;
			odata.precision(m_precision);
			odata << m_data[i][m_cols];
			std::stringstream pdata;
			pdata.precision(m_precision);
			pdata << evalROM(m_data[i]);
			fileP << odata.str() << "\t" << pdata.str() << "\n";
		}
		fileP.close();

		tlog += "Prediction file: writing completed!\n";
		std::cout << "Prediction file: writing completed!" << std::endl;
	}
	else
	{
		tlog += "Can not execute writePrediction, please read data first: setData function.\n";
		std::cout << "Can not execute writePrediction, please read data first: setData function." << std::endl;
	}

	return;
}

// write a file (evalname) with the ROM's prediction evaluated at the values in file dataname
void Twinkle::writePrediction(std::string dataname, std::string evalname)
{
	m_data.clear();
	setData(dataname);

	if ((!m_data.empty()) && (m_tmatrix.size() > 0) && m_cols == (int)m_tmatrix[0].m_fnD.size())
	{
		std::ofstream fileP(evalname);
		fileP << "data" << "\t" << "prediction" << "\n";
		for (int i = 0; i < m_rows; i++)
		{
			if (isInBounds(m_data[i]))
			{
				std::stringstream odata;
				odata.precision(m_precision);
				odata << m_data[i][m_cols];
				std::stringstream pdata;
				pdata.precision(m_precision);
				pdata << evalROM(m_data[i]);
				fileP << odata.str() << "\t" << pdata.str() << "\n";
			}
			else
			{
				tlog += "Warning, data out of ROM boundaries. The prediction will not be evaluated.\n";
				std::cout << "Warning, data out of ROM boundaries. The prediction will not be evaluated." << std::endl;
			}
		}
		fileP.close();

		tlog += "Prediction file: writing completed!\n";
		std::cout << "Prediction file: writing completed!" << std::endl;
	}
	else
	{
		if (m_cols != (int)m_tmatrix[0].m_fnD.size())
		{
			tlog += "The number of independent variables in the data file does not match the ROM file.\n";
			std::cout << "The number of independent variables in the data file does not match the ROM file." << std::endl;
		}
		else
		{
			tlog += "Can not execute writePrediction, please read data first: setData function.\n";
			std::cout << "Can not execute writePrediction, please read data first: setData function." << std::endl;
		}
	}

	return;
}

// write a file (evalname) with the ROM's prediction evaluated at the values in file dataname, using only the first nterms Terms
void Twinkle::writePrediction(std::string dataname, int nterms, std::string evalname)
{
	m_data.clear();
	setData(dataname);

	if ((!m_data.empty()) && (m_tmatrix.size() > 0))
	{
		std::ofstream fileP(evalname);
		fileP << "data" << "\t" << "prediction" << "\n";
		for (int i = 0; i < m_rows; i++)
		{
			if (isInBounds(m_data[i]))
			{
				std::stringstream odata;
				odata.precision(m_precision);
				odata << m_data[i][m_cols];
				std::stringstream pdata;
				pdata.precision(m_precision);
				pdata << evalROM(m_data[i], nterms);
				fileP << odata.str() << "\t" << pdata.str() << "\n";
			}
			else
			{
				tlog += "Warning, data out of ROM boundaries. The prediction will not be evaluated.\n";
				std::cout << "Warning, data out of ROM boundaries. The prediction will not be evaluated." << std::endl;
			}
		}
		fileP.close();

		tlog += "Prediction file: writing completed!\n";
		std::cout << "Prediction file: writing completed!" << std::endl;
	}
	else
	{
		tlog += "Can not execute writePrediction, please read data first: setData function.\n";
		std::cout << "Can not execute writePrediction, please read data first: setData function." << std::endl;
	}

	return;
}

// write a file with the ROM's log
void Twinkle::writeLog(std::string fileName)
{
	std::ofstream file(fileName);

	if(!m_tmatrix.empty())
	{
		tlog += "Log file: writing completed!\n";
		std::cout << "Log file: writing completed!" << std::endl;

		file << tlog;
	}
	else
	{
		tlog += "Can not execute writeLog, please calculate ROM first: calcROM function.\n";
		std::cout << "Can not execute writeLog, please calculate ROM first: calcROM function." << std::endl;

		file << tlog;
	}

	file.close();
	tlog = "";

	return;
}

// check if all values of the independent variable in x fall within the discretization net's boundaries
bool Twinkle::isInBounds(std::vector<double> x)
{
	bool output = true;

	for (int i = 0; i < m_cols; i++)
	{
		if (x[i] < m_tmatrix[0].m_fnD[i].m_discretization[0] || x[i] > m_tmatrix[0].m_fnD[i].m_discretization[m_tmatrix[0].m_fnD[i].m_discretization.size() - 1])
		{
			output = false;
		}
	}

	return output;
}

// evaluate the ROM using all Terms at the values of the independent variables in xs
double Twinkle::evalROM(std::vector<double> xs)
{
	double output = 0;

	for (unsigned int i = 0; i < m_tmatrix.size(); i++)
	{
		output += m_tmatrix[i].evalTerm(xs);
	}

	return output;
}

// evaluate the ROM using all Terms at the values of the independent variables in xs, using only the first nterms Terms
double Twinkle::evalROM(std::vector<double> xs, int nterms)
{
	double output = 0;

	for (unsigned int i = 0; i < std::fmin(nterms, m_tmatrix.size()); i++)
	{
		output += m_tmatrix[i].evalTerm(xs);
	}

	return output;
}

