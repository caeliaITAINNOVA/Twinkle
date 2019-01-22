//*******************************************************************************************//
//* TWINKLE																				 	*//
//* Version: 0.2																		 	*//
//* copyright (c) 2018 Instituto Tecnológico de Aragón (ITAINNOVA) (España)				 	*//
//* Date: 4 December 2018																 	*//
//* Author: Instituto Tecnológico de Aragón (ita@itainnova.es) - R.Rodríguez and V.Zambrano *//
//* All rights reserved																	 	*//
//*******************************************************************************************//

#include "runTwinkle.h"

using namespace ROM;

// default constructor for the example class runTwinkle
runTwinkle::runTwinkle() {}

// main function
int main(int argc, char *argv[])
{
	try
	{
		// initialization of ROM variables with default values
		int nPoints = 5;
		std::vector<int> vPoints;
		double gTol = 1.0E-2;
		double tTol = 1.0E-2;
		int nIter = 50;
		int nTerms = 5;
		int nZIter = 5;
		int initMode = 0; // initialization mode. 0: random, 1: linear
		int precision = 6;
		std::string filename;
		std::string filesuffix;
		std::string romname;
		std::string evalname;
		std::string discfile;

		// processing of command line arguments and values
		bool isEvaluation = false;
		bool isVPoints = false;
		bool isMPoints = false;
		for (int i = 1; i < argc; i++)
		{
			std::string argtype = argv[i];

			if (argtype == "-file")
			{
				filename = argv[i + 1];
			}
			else if (argtype ==  "-init")
			{
				initMode = std::stoi(argv[i + 1]);
			}
			else if (argtype == "-precision")
			{
				precision = std::stoi(argv[i + 1]);
			}
			else if (argtype == "-ziter")
			{
				nZIter = std::stoi(argv[i + 1]);
			}
			else if (argtype ==  "-terms")
			{
				nTerms = std::stoi(argv[i + 1]);
			}
			else if (argtype ==  "-alsiter")
			{
				nIter = std::stoi(argv[i + 1]);
			}
			else if (argtype ==  "-ttol")
			{
				tTol = std::stod(argv[i + 1]);
			}
			else if (argtype == "-gtol")
			{
				gTol = std::stod(argv[i + 1]);
			}
			else if (argtype == "-npoints")
			{
				nPoints = std::stoi(argv[i + 1]);
			}
			else if (argtype == "-vpoints")
			{
				discfile = argv[i + 1];
				isVPoints = true;
			}
			else if (argtype == "-mpoints")
			{
				discfile = argv[i + 1];
				isMPoints = true;
			}
			else if (argtype ==  "-out")
			{
				filesuffix = argv[i + 1];
			}
			else if (argtype == "-rom")
			{
				romname = argv[i + 1];
			}
			else if (argtype == "-eval")
			{
				evalname = argv[i + 1];
				isEvaluation = true;
			}
			else if (argtype == "-help")
			{
				runTwinkle ex;
				ex.printHelp();
				return 0;
			}
			else
			{
				std::cout << "Invalid argument type, please type 'runTwinkle.exe -help' for help." << std::endl;
				return 0;
			}
			i++;
		}

		// create an instance of Twinkle
		Twinkle twin;
		runTwinkle ex;
		if (!isEvaluation)
		{
			if (filename != "")
			{
				if (nTerms == 0)
				{
					nTerms = 5;
				}
				twin.setData(filename);
				twin.setParams(gTol, tTol, nIter, nZIter, precision);
				if (!isVPoints && !isMPoints)
				{
					if (nPoints <= 0)
					{
						twin.genExactDiscret();
					}
					else
					{
						twin.genDiscret(nPoints);
					}
				}
				else if (isVPoints)
				{
					twin.genDiscret(discfile);
				}
				else if (isMPoints)
				{
					twin.genNonHomogDiscret(discfile);
				}
				else
				{
					std::cout << "Missing discretization net, please type 'runTwinkle.exe -help' for help." << std::endl;
					return 0;
				}
				twin.calcROM(nTerms, initMode);
				twin.writeROM("Results_" + filesuffix + ".txt");
				twin.writePrediction("Prediction_" + filesuffix + ".txt");
				twin.writeLog("Log_" + filesuffix + ".txt");
			}
			else
			{
				std::cout << "Missing data file name, please type 'runTwinkle.exe -help' for help." << std::endl;
				return 0;
			}
		}
		else
		{
			// create an instance of Twinkle for ROM evaluation
			if (romname != "" && evalname != "")
			{
				Twinkle twineval(romname);
				if (nTerms == 0)
				{
					twineval.writePrediction(evalname, "Prediction_" + filesuffix + ".txt");
					twin.writeLog("Log_" + filesuffix + ".txt");
				}
				else
				{
					twineval.writePrediction(evalname, nTerms, "Prediction_" + filesuffix + ".txt");
					twin.writeLog("Log_" + filesuffix + ".txt");
				}
			}
			else
			{
				std::cout << "Missing ROM and data files' names, please type 'runTwinkle.exe -help' for help." << std::endl;
				return 0;
			}

		}
		return 0;
	}
	catch (int e)
	{
		std::cout << "Invalid argument value, please type 'runTwinkle.exe -help' for help." << std::endl;
		return 0;
	}
	catch (...)
	{
		std::cout << "Exiting with an error." << std::endl;
		return 0;
	}
}

// print usage commands
void runTwinkle::printHelp()
{
	std::cout << "Arguments usage (if not specified default values will be applied): \n"
			"For ROM calculation: \n"
			"\t-file [file_path\\file_name.file_extension where the extension can be csv, txt or similar] - this field is mandatory \n"
			"\t-init [initialization type: 0 for random, 1 for linear] - default value 0 (random) \n"
			"\t-precision [output numerical precision] - default value 6 \n"
			"\t-ziter [number of initializations] - default value 5 \n"
			"\t-terms [number of terms] - default value 5 \n"
			"\t-alsiter [number of alternate least squares iterations over shape functions] - default value - 50 \n"
			"\t-ttol [term tolerance for convergence] - default value 1.0E-2 \n"
			"\t-gtol [global tolerance for convergence] - default value 1.0E-2 \n"
			"\t-npoints [number of interpolation points: there must be at least 2 interpolation points; if 0 is used the interpolation will be performed using all parameters' values] - not compatible with -vpoints and -mpoints flags - default value 5 \n"
			"\t-vpoints [file_path\\file_name.file_extension where the extension can be csv, txt or similar - file containing a column with the number of interpolation points, one different value for each variable - there must be at least 2 interpolation points; if 0 is used the interpolation will be performed using all parameters' values] - not compatible with -npoints and -mpoints flags \n"
			"\t-mpoints [file_path\\file_name.file_extension where the extension can be csv, txt or similar - file containing a number of rows equal to the number of independent variables, each of them containing the values to be used for the discretization net - there must be at least 2 interpolation points] - not compatible with -npoints and -vpoints flags \n"
			"\t-out [output files suffix, e.g. Prediction_<out>] \n"
			"For ROM evaluation: \n"
			"\t-rom [rom_file_path\\rom_file_name.txt] - this field is mandatory \n"
			"\t-eval [data_file_path\\data_file_name.file_extension where the extension can csv, txt or similar] - this field is mandatory \n"
			"\t-terms [number of terms] - default value all \n"
			"\t-out [output files suffix, i.e. Prediction_<out>] \n"
			<< std::endl;
}
