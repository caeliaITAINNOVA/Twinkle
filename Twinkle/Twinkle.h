//*******************************************************************************************//
//* TWINKLE																				 	*//
//* Version: 0.2																		 	*//
//* copyright (c) 2018 Instituto Tecnológico de Aragón (ITAINNOVA) (España)				 	*//
//* Date: 4 December 2018																 	*//
//* Author: Instituto Tecnológico de Aragón (ita@itainnova.es) - R.Rodríguez and V.Zambrano *//
//* All rights reserved																	 	*//
//*******************************************************************************************//

#pragma once
#ifndef H__TWINKLE__H
#define H__TWINKLE__H

#include "Term.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

// class Twinkle description
class ROM::Twinkle
{
public:
	Twinkle();
	Twinkle(std::string);

	void calcROM(int, int);
	void setParams(double, double, int, int, int);
	void writeROM(std::string);
	void writePrediction(std::string);
	void writePrediction(std::string, std::string);
	void writePrediction(std::string, int, std::string);
	void writeLog(std::string);
	void setData(std::string);
	void setData(std::vector<std::vector<double>>);
	void genDiscret(int);
	void genDiscret(std::vector<int>);
	void genDiscret(std::string);
	void genNonHomogDiscret(std::string);
	void genExactDiscret();
	std::vector<double> genExactDiscret(int);
	double evalROM(std::vector<double>);
	double evalROM(std::vector<double>, int);

private:
	void fillNtensor(int, std::vector<Term>);
	void calcZs(std::vector<Term>&, int, int);
	void normZs(int);
	void calcAlphas();
	void setAlphas(std::vector<Term>, std::vector<double>);
	void fillCurrEval(double, std::vector<Term>);
	double fillCurrEval(double);
	double getMaxrDiff(std::vector<vec>, std::vector<vec>);
	std::vector<vec> getrcurr(Term&);
	std::vector<vec> copyRs(std::vector<vec>);
	std::vector<std::vector<double>> setBounds();
	bool isDataInDiscret();
	bool isInBounds(std::vector<double>);
	bool containsValue(std::vector<double>, double);

	bool m_bConverg;
	int m_maxIter;
	int m_maxZIter;
	int m_rows;
	int m_cols;
	int m_precision;
	double m_Tol;
	double m_globalTol;
	double m_termTol;
	double m_errorMin;
	std::vector<double> m_currEval;
	std::vector<vec> m_rcurr;
	std::vector<Term> m_tmatrix;
	std::vector<std::vector<double>> m_data;
	std::vector<std::vector<double>> m_discret;
	cube m_Ntensor;
};
#endif // !H__TWINKLE__H
