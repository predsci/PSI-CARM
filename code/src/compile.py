#!/usr/bin/python
import os
import string
os.system("make")
os.system("R CMD SHLIB detcovid.o ignbin.o toms343.o utils.o")
os.system("R CMD SHLIB stochcovid.o ignbin.o toms343.o utils.o")
os.system("R CMD SHLIB mcmc.o detcovid.o ignbin.o toms343.o utils.o")
os.system("R CMD SHLIB tanhrt.o")
os.system("/bin/rm -rf *.o")
