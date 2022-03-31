# PSI-CARN
Covid Arbitrary Reproductive number (CARN) code applied to the U.S. and the world during the first wave of the COVID-19 pandemic.
##
## Structure of Repository and Compilation Instructions
##

### Code
This repository includes fortran code. 
After cloning the code navigate to the 'code/src' directory

%cd code/src

and use the python script to compile

%./compile.py

This should create four dynamic libraries: mcmc.so detcovid.so stochcovid.so tanhrt.so


R scripts for the U.S and the world can be found in the two sub-directories:

code/usa 

code/world


In each sub-directory there are R scripts for a 'calendar time' and 'pandemic time' calculations:

code/usa: usa_mcmc_fit_calendar_time.R	usa_mcmc_fit_pandemic_time.R

code/world: world_mcmc_fit_calendar_time.R	world_mcmc_fit_pandemic_time.R

In each sub-directory there is also a python script that is set to launch multiple R threads so that each one calculates a subset of the locations. The user can change the number of threads and should also set the script to call the desired R script, i.e. calendar or pandemic time.

### Data
The data sub-directory includes two 'rds' files: World_data.rds State_data.rds. The first is used by the 'world' calculations and the second by the U.S. calculations.
Please note that the R scripts assume that the directory structure is as this repository.   If you would like to change this you will need to modify the lines of code that load the fortran libraries and the data.

### Required Packages
squire, countrycode, EpiEstim, incidence, coda, ggplot2, grid, gridExtra
