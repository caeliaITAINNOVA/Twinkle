# Twinkle
TWINKLE: A Digital-Twin-Building Kernel for Real-Time Computer-Aided Engineering

# About
TWINKLE is a library for building families of solvers to perform Canonical Polyadic Decomposition (CPD) of tensors. The common characteristic of these solvers is that the data structure supporting the tuneable solution strategies is based on a Galerkin projection of the phase space. This allows processing and recovering tensors described by highly sparse and unstructured data. For achieving high performance, TWINKLE is written in C++ and uses the Armadillo open source library for linear algebra and scientific computing, based on LAPACK (Linear Algebra PACKage) and BLAS (Basic Linear Algebra Subprograms) routines. The library has been implemented keeping in mind its future extensibility and adaptability to fulfil the different users' needs in academia and industry regarding Reduced Order Modelling (ROM) and data analysis by means of tensor decomposition. It is especially focused on post-processing data from Computer-Aided-Engineering (CAE) simulation tools.

# Licence
This product includes software developed at [ITAINNOVA](http://www.itainnova.es).

Twinkle has two kind of licenses: commercial and open-source.

#### * Commercial license
If you want to use Twinkle in a commercial way (any use out of open source applications under [GNU GPL license v3](https://www.gnu.org/licenses/gpl-3.0.html)), please contact with [ITAINNOVA](http://www.itainnova.es) (caelia@itainnova.es).

#### * Open source license
If you are creating an open source application under a license compatible with the [GNU GPL license v3](https://www.gnu.org/licenses/gpl-3.0.html) you may use Twinkle under its terms and conditions.

# Folders overview
Data: contains a data sample for code usage\
linux: contains the precompiled library and the executable file for linux OS – CentOS Linux 7 (Core)\
runTwinkle: contains the code of the executable file\
Twinkle: contains the code of the library\
win64: contains the precompiled library and the executable file for win64 OS together with a README.txt for library usage\
NOTICE.txt\
README.md

# Getting started
In order to compile Twinkle it is necessary to previously install Armadillo. Armadillo is an open source library for linear algebra and scientific computing, based on LAPACK (Linear Algebra PACKage) and BLAS (Basic Linear Algebra Subprograms) routines. It is therefore necessary to install all three packages following the instructions on Armadillo web page.\
If you prefer using precompiled files, please select a folder corresponding to your OS (available for linux and win64). Otherwise compile the library and the executable using the following commands respectively (see more under gcc online documentation).

_g++ -c -shared Twinkle.cpp -o Twinkle.so -I"Path_to/armadillo_bits" -I"Path_to/armadillo" -L"/usr/local/lib" -includearmadillo -O2 -g3 -Wall -fmessage-length=0_ 

_g++ -o runTwinkle runTwinkle.cpp Path_to_twinkle/*.cpp -I Path_to_twinkle/ -I Path_to/armadillo -L /usr/local/lib -llapack -lblas –lgfortran_ 

The executable file, that uses Twinkle library, provides a help guide for software usage. It can be visualized by typing the -help command.

# Running the code
TWINKLE library can be used directly through the runTwinkle.exe application. This is a command-line application whose usage is as follows:\
Arguments usage (if not specified default values will be applied):

For ROM calculation:
* _-file_ [file_path\file_name.file_extension where the extension can be csv, txt or similar] - this field is mandatory 
* _-init_ [initialization type: 0 for random, 1 for linear] - default value 0 (random)
* _-precision_ [output numerical precision] - default value 6
* _-ziter_ [number of initializations] - default value 5 
* _-terms_ [number of terms] - default value 5 
* _-alsiter_ [number of alternate least squares iterations over shape functions] - default value - 50 
* _-ttol_ [term tolerance for convergence] - default value 1.0E-2 
* _-gtol_ [global tolerance for convergence] - default value 1.0E-2 
* _-npoints_ [number of interpolation points: there must be at least 2 interpolation points; if 0 is used the interpolation will be performed using all parameters' values] - not compatible with -vpoints and -mpoints flags 
* _-vpoints_ [file_path\file_name.file_extension where the extension can be csv, txt or similar - file containing a column with the number of interpolation points, one different value for each variable - there must be at least 2 interpolation points; if 0 is used the interpolation will be performed using all parameters' values] - not compatible with -npoints and -mpoints flags 
* _-mpoints_ [file_path\file_name.file_extension where the extension can be csv, txt or similar - file containing a number of rows equal to the number of independent variables, each of them containing the values to be used for the discretization net - there must be at least 2 interpolation points] - not compatible with -npoints and -vpoints flags 
* _-out_ [output files suffix, e.g. Prediction_<out>] 
  
For ROM evaluation: 
* _-rom_ [rom_file_path\rom_file_name.txt] - this field is mandatory 
* _-eval_ [data_file_path\data_file_name.file_extension where the extension can csv, txt or similar] - this field is mandatory 
* _-terms_ [number of terms] - default value all 
* _-out_ [output files suffix, i.e. Prediction_<out>] 
  
For example:\
_runTwinkle.exe -file Rastrigin.csv -terms 10 -npoints 5 -out rastrigin_
