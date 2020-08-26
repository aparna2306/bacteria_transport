# bacteria_transport
Data including R codes from the transport of bacteria in columns with and without dissolved oxygen and nutrient


This README explains, the purpose and contents of the files within each folder titled E2-E3 (sub folders E3.1-E3.4). Each of the folders (and sub-folders) contain four files. They are the following: 
	- source.R (Source Code)
	- boundary.f95 (Fortran function file)
	- paper_data.xlsx or obs_data.xlsx (Experimental data)
	- parameters.xlsx (Parameter definitions)

By loading the source file into R console, (and ensuring that the other three files are present in the folder path), the “Run” icon can be clicked to get the result for each of the conditions. The contents, and purpose of each file in explained in the sections below.

A. Source code
This file contains the source code. In general, this code can be run, as long as the other 3 files (boundary.R, paper_data.xlsx, parameters.xlsx) are added to the folder path. The source file is broken into several subsections, they are as below:

	- All the packages needed to run the code are initialized. The system path has to be set to the folder which contains the other three files. (This is done using the library() function)

	- The guess values of the parameters are obtained from the “pars” sheet of parameters.xlsx. (usually the first two rows of the “pars” are the fitting parameters)

	- The model function with the ode method is defined within the function Model_fun(p), where p is the fitting parameter(s). The rodeo package generates a model (defined under “model” in source code). Using the function model$compile() the Fortran function file (containing the boundary condition) is implemented into the model. The variables and parameter values are separately defined for each grid, and the model$setVars() and model$setPars(), are used to initialize the variables and parameter values. Finally, the ode solver is implemented into model$dynamics() which is used to solve the ODE for the given time period defined in the function.  

	- The model cost or residual function is defined in the function Model_Cost, where the in-built FME function modCost(), generates the residual between the model and the experimental observations.

	- The “Fit” is the final result where we can get the fit of the model to the data using the modFit() function file. The Model_fun() and Model_cost() are used to calculate the residual for the model at each step of the parameter optimization. The modFit() contains options to choose various methods that can be used for parameter optimization. 

	- Finally, the plotting of the data using ggplot() is defined in the source code, for easy and clear visualization of the results. 

	B. Fortran function file
 The Fortran function file is defined in the source code using model$compile(). It contains the time-varying boundary conditions that are defined using linear interpolation. 

	C. Experimental Data
The processed experimental data with the time values data and the dimensionalized concentration values are shown. There is an additional column showing the standard deviation associated with the particular data point which arises from technical replicates. 

	D. Parameter definitions

This is an excel file with separate sheets for: 
	- variables: vars 
	- parameters: pars 
	- functions:funcs 
	- processes: pros
	- stoichiometric matrix: stoi
These are separately used in the model construction in the source file (“model”). For further details on the definitions refer to the help section for rodeo by typing “?rodeo” in R. 




